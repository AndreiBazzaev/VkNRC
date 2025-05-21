#include "Dependency.hpp"
#include <algorithm>
#include <unordered_set>
#include <iostream> // For debug output
#include <cassert>  // For debug assertions

namespace myvk_rg_executor {

	struct Dependency::PassVisitor {
		const Args& args;
		Dependency& self;

		void operator()(const GraphicsPassBase* p_pass) const { handle_pass(p_pass); }
		void operator()(const ComputePassBase* p_pass) const { handle_pass(p_pass); }
		void operator()(const TransferPassBase* p_pass) const { handle_pass(p_pass); } // example

		// Catch-all for unhandled types
		void operator()(const auto* unused) const {
			assert(false && "Unhandled pass type in PassVisitor");
		}

	private:
		template <typename PassType>
		void handle_pass(const PassType* p_pass) const {
			//std::cout << "[PassVisitor] Visiting pass: " << p_pass->GetName() << '\n';

			for (const auto& it : p_pass->GetInputPoolData()) {
				const InputBase* p_input = it.second.template Get<InputBase>();
				self.get_dep_info(p_input).p_pass = p_pass;
				self.get_dep_info(p_pass).inputs.push_back(p_input);

				InputVisitor visitor{ args, self, p_pass, p_input };
				p_input->GetInputAlias().Visit(visitor);
			}
		}
	};

	struct Dependency::InputVisitor {
		const Args& args;
		Dependency& self;
		const PassBase* p_pass;
		const InputBase* p_input;

		void operator()(const OutputAlias auto* p_output_alias) const {
			std::cout << "[InputVisitor] Traversing OutputAlias\n";

			const InputBase* p_src_input = self.traverse_output_alias(args, *p_output_alias);
			if (!p_src_input) Throw(error::PassNotDAG{});

			const PassBase* p_src_pass = self.get_dep_info(p_src_input).p_pass;
			const ResourceBase* p_resource = self.get_dep_info(p_src_input).p_resource;

			if (!p_resource || p_src_pass == p_pass) {
				Throw(error::PassNotDAG{});
			}

			self.get_dep_info(p_input).p_resource = p_resource;
			self.m_pass_graph.AddEdge(p_src_pass, p_pass, PassEdge{ p_src_input, p_input, p_resource });
		}

		void operator()(const RawAlias auto* p_raw_alias) const {
			std::cout << "[InputVisitor] Handling RawAlias\n";

			const ResourceBase* p_resource = args.collection.FindResource(p_raw_alias->GetSourceKey());
			if (!p_resource) Throw(error::PassNotDAG{});

			self.m_resource_graph.AddVertex(p_resource);
			self.get_dep_info(p_input).p_resource = p_resource;

			ResourceVisitor visitor{ args, self, p_pass, p_input, p_resource };
			p_resource->Visit(visitor);
		}

		void operator()(const auto* unused) const {
			assert(false && "Unhandled input alias type in InputVisitor");
		}
	};

	struct Dependency::ResourceVisitor {
		const Args& args;
		Dependency& self;
		const PassBase* p_pass;
		const InputBase* p_input;
		const ResourceBase* p_resource;

		void operator()(const CombinedResource auto* p_combined_resource) const {
			std::cout << "[ResourceVisitor] Handling CombinedResource\n";

			for (const auto& src_alias : p_combined_resource->GetSubAliases()) {
				const InputBase* p_src_input = self.traverse_output_alias(args, src_alias);
				const PassBase* p_src_pass = self.get_dep_info(p_src_input).p_pass;
				const ResourceBase* p_sub_resource = self.get_dep_info(p_src_input).p_resource;

				if (!p_sub_resource || p_src_pass == p_pass) {
					Throw(error::PassNotDAG{});
				}

				self.m_pass_graph.AddEdge(p_src_pass, p_pass, PassEdge{ p_src_input, p_input, p_sub_resource });
				self.m_resource_graph.AddEdge(p_resource, p_sub_resource, {});
			}
		}

		void operator()(const auto*) const {
			std::cout << "[ResourceVisitor] Handling default resource\n";

			self.m_pass_graph.AddEdge(nullptr, p_pass, PassEdge{ nullptr, p_input, p_resource });
		}
	};

	Dependency Dependency::Create(const Args& args) {
		args.collection.ClearInfo(&PassInfo::dependency, &InputInfo::dependency, &ResourceInfo::dependency);

		Dependency g = {};
		for (const auto& it : args.render_graph.GetResultPoolData()) {
			it.second.Visit([&](const auto* p_alias) {
				g.traverse_output_alias(args, *p_alias);
				});
		}
		g.tag_resources(args);
		g.add_war_edges();
		g.sort_passes();
		g.get_pass_relation();
		g.get_resource_relation();
		g.add_image_read_edges();
		return g;
	}

	const InputBase* Dependency::traverse_output_alias(const Args& args, const OutputAlias auto& output_alias) {
		const PassBase* p_src_pass = args.collection.FindPass(output_alias.GetSourcePassKey());
		const InputBase* p_src_input = args.collection.FindInput(output_alias.GetSourceKey());

		if (!p_src_pass || !p_src_input) {
			std::cerr << "[Error] Invalid output alias: missing pass or input\n";
			return nullptr;
		}

		//std::cout << "[traverse_output_alias] From pass: " << p_src_pass->GetName() << '\n';
		traverse_pass(args, p_src_pass);
		return p_src_input;
	}

	void Dependency::traverse_pass(const Args& args, const PassBase* p_pass) {
		if (m_pass_graph.HasVertex(p_pass)) return;

		static thread_local std::unordered_set<const PassBase*> pass_stack;
		if (pass_stack.contains(p_pass)) {
			Throw(error::PassNotDAG{});
		}

		//std::cout << "[traverse_pass] Visiting pass: " << p_pass->GetName() << '\n';

		m_pass_graph.AddVertex(p_pass);
		pass_stack.insert(p_pass);

		try {
			p_pass->Visit(PassVisitor{ args, *this });
			pass_stack.erase(p_pass);
		}
		catch (...) {
			pass_stack.erase(p_pass);
			throw;
		}
	}
struct AccessEdgeInfo {
	std::vector<std::size_t> reads;
	std::optional<std::size_t> opt_write;
};
void Dependency::add_war_edges() {
	auto view = m_pass_graph.MakeView(kAnyFilter, kPassEdgeFilter<PassEdgeType::kBarrier>);

	for (const PassBase *p_pass : view.GetVertices()) {
		std::unordered_map<const ResourceBase *, AccessEdgeInfo> access_edges;

		for (auto [_, e, edge_id] : view.GetOutEdges(p_pass)) {
			auto &info = access_edges[e.p_resource];
			if (UsageIsReadOnly(e.p_dst_input->GetUsage())) {
				info.reads.push_back(edge_id);
			} else {
				// Forbid multiple writes
				if (info.opt_write)
					Throw(error::MultipleWrite{.alias = e.p_dst_input->GetInputAlias()});

				info.opt_write = edge_id;
			}
		}

		for (const auto &[p_resource, info] : access_edges) {
			// If no writes or no reads, skip
			if (!info.opt_write || info.reads.empty())
				continue;

			std::size_t write_id = *info.opt_write;
			const InputBase *p_write_dst = m_pass_graph.GetEdge(write_id).p_dst_input;

			// Add edges from read to write
			for (std::size_t read_id : info.reads)
				m_pass_graph.AddEdge(m_pass_graph.GetToVertex(read_id), m_pass_graph.GetToVertex(write_id),
				                     PassEdge{.opt_p_src_input = m_pass_graph.GetEdge(read_id).p_dst_input,
				                              .p_dst_input = p_write_dst,
				                              .p_resource = p_resource});

			// Turn to Indirect WAW Edge
			m_pass_graph.GetEdge(write_id).type = PassEdgeType::kIndirectWAW;
		}
	}
}

void Dependency::sort_passes() {
	// Exclude nullptr Pass, use Barrier edges only
	auto view = m_pass_graph.MakeView([](const PassBase *p_pass) -> bool { return p_pass; },
	                                  kPassEdgeFilter<PassEdgeType::kBarrier>);

	auto kahn_result = view.KahnTopologicalSort();

	if (!kahn_result.is_dag)
		Throw(error::PassNotDAG{});

	// Assign topo-id to passes
	m_passes = std::move(kahn_result.sorted);
	for (std::size_t topo_id = 0; const PassBase *p_pass : m_passes)
		get_dep_info(p_pass).topo_id = topo_id++;
}

void Dependency::tag_resources(const Args &args) {
	// Validate and Tag Resources
	for (const ResourceBase *p_resource : m_resource_graph.GetVertices()) {
		p_resource->Visit(overloaded(
		    [&](const ExternalResource auto *p_resource) {
			    // External Resource should not have parent resource
			    if (!m_resource_graph.GetInEdges(p_resource).empty())
				    Throw(error::ResourceExtParent{.key = p_resource->GetGlobalKey()});
		    },
		    [](auto &&) {}));

		m_resources.push_back(p_resource);
	}

	// Resolve Resource Tree
	auto find_trees_result = m_resource_graph.FindTrees(
	    [](const ResourceBase *p_root, const ResourceBase *p_sub) { get_dep_info(p_sub).p_root_resource = p_root; });

	if (!find_trees_result.is_forest)
		Throw(error::ResourceNotTree{});

	m_root_resources = std::move(find_trees_result.roots);

	// Assign Root ID to resources
	for (std::size_t root_id = 0; const ResourceBase *p_root_resource : m_root_resources)
		get_dep_info(p_root_resource).root_id = root_id++;
	for (const ResourceBase *p_resource : m_resources)
		if (!IsRootResource(p_resource))
			get_dep_info(p_resource).root_id = get_dep_info(GetRootResource(p_resource)).root_id;
}

void Dependency::get_pass_relation() {
	// Exclude nullptr Pass, use Barrier edges only
	auto view = m_pass_graph.MakeView([](const PassBase *p_pass) -> bool { return p_pass; },
	                                  kPassEdgeFilter<PassEdgeType::kBarrier>);

	m_pass_relation =
	    view.TransitiveClosure(GetPassTopoID, [this](std::size_t topo_id) { return GetTopoIDPass(topo_id); });
}

void Dependency::get_resource_relation() {
	Relation resource_pass_access{GetRootResourceCount(), GetPassCount()};
	// Tag access bits of root resources
	for (std::size_t topo_id = 0; const PassBase *p_pass : m_passes) {
		for (const InputBase *p_input : GetPassInputs(p_pass))
			resource_pass_access.Add(GetResourceRootID(GetInputResource(p_input)), topo_id);
		++topo_id;
	}

	// Pass < Resource
	Relation pass_resource_relation{GetPassCount(), GetRootResourceCount()};
	for (std::size_t pass_topo_id = 0; pass_topo_id < GetPassCount(); ++pass_topo_id)
		for (std::size_t root_id = 0; root_id < GetRootResourceCount(); ++root_id) {
			// Pass < Resource <==> forall x in Resource Access Passes, Pass < x
			if (m_pass_relation.All(pass_topo_id, resource_pass_access.GetRowData(root_id)))
				pass_resource_relation.Add(pass_topo_id, root_id);
		}

	// Resource > Pass
	Relation resource_pass_relation = pass_resource_relation.GetInversed();

	// Resource < Resource
	m_resource_relation.Reset(GetRootResourceCount(), GetRootResourceCount());
	for (std::size_t l_root_id = 0; l_root_id < GetRootResourceCount(); ++l_root_id)
		for (std::size_t r_root_id = 0; r_root_id < GetRootResourceCount(); ++r_root_id) {
			// Resource_L < Resource_R <==> forall x in Resource_L Access Passes, Resource_R > x
			if (resource_pass_relation.All(r_root_id, resource_pass_access.GetRowData(l_root_id)))
				m_resource_relation.Add(l_root_id, r_root_id);
		}
}

void Dependency::add_image_read_edges() {
	// Add extra image read edges, since multiple reads to the same image can break merging if, for example the first
	// (topo_id) reads as Input attachment and the second reads as Sampler

	// Only process Barrier & Image access Edges
	auto view = m_pass_graph.MakeView(Dependency::kAnyFilter, [](const Dependency::PassEdge &e) -> bool {
		return e.type == PassEdgeType::kBarrier && e.p_resource->GetType() == ResourceType::kImage;
	});

	for (const PassBase *p_pass : view.GetVertices()) {
		std::vector<std::size_t> out_edge_ids;
		for (auto [_, _1, edge_id] : view.GetOutEdges(p_pass))
			out_edge_ids.push_back(edge_id);

		// Sort Output Edges with Topological Order of its 'To' Vertex
		std::ranges::sort(out_edge_ids, [&](std::size_t l_edge_id, std::size_t r_edge_id) {
			const PassBase *p_l_pass = m_pass_graph.GetToVertex(l_edge_id),
			               *p_r_pass = m_pass_graph.GetToVertex(r_edge_id);
			return Dependency::GetPassTopoID(p_l_pass) < Dependency::GetPassTopoID(p_r_pass);
		});

		std::unordered_map<const ResourceBase *, std::size_t> access_edges;

		for (std::size_t edge_id : out_edge_ids) {
			const auto &e = m_pass_graph.GetEdge(edge_id);
			const ResourceBase *p_resource = e.p_resource;

			auto it = access_edges.find(p_resource);
			if (it == access_edges.end())
				access_edges[p_resource] = edge_id;
			else {
				std::size_t prev_edge_id = it->second;
				it->second = edge_id;

				const auto &prev_e = m_pass_graph.GetEdge(prev_edge_id);
				m_pass_graph.AddEdge(m_pass_graph.GetToVertex(prev_edge_id), m_pass_graph.GetToVertex(edge_id),
				                     {.opt_p_src_input = prev_e.p_dst_input,
				                      .p_dst_input = e.p_dst_input,
				                      .p_resource = p_resource,
				                      .type = PassEdgeType::kImageRead});
			}
		}
	}
}

} // namespace myvk_rg_executor
