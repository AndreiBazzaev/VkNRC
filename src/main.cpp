#include <myvk/FrameManager.hpp>
#include <myvk/GLFWHelper.hpp>
#include <myvk/ImGuiHelper.hpp>

#include <spdlog/spdlog.h>

#include "rg/NRCRenderGraph.hpp"

constexpr uint32_t kFrameCount = 3, kWidth = 1280, kHeight = 720;

int main(int argc, char **argv) {
	--argc, ++argv;
	// Scene scene = Scene::LoadOBJShapeInstanceSAH(argv[0], 7); // at most 128 instances
	Scene scene;
	if (argc == 0) {
		spdlog::error("No OBJ file, loading default");
		// Put your scene name here for faster testing
		scene = Scene::LoadOBJSingleInstance("living_room/living_room.obj");
	}
	else
		scene = Scene::LoadOBJSingleInstance(argv[0]);
	GLFWwindow *window = myvk::GLFWCreateWindow("VkNRC", kWidth, kHeight, true);
	if (scene.Empty())
		return EXIT_FAILURE;
	spdlog::info("Loaded {} Vertices, {} Texcoords, {} Materials, {} Instances", scene.GetVertices().size(),
	             scene.GetTexcoords().size(), scene.GetMaterials().size(), scene.GetInstances().size());

	auto instance = myvk::Instance::CreateWithGlfwExtensions();
	myvk::Ptr<myvk::Queue> generic_queue, compute_queue;
	myvk::Ptr<myvk::PresentQueue> present_queue;
	myvk::Ptr<myvk::PhysicalDevice> physical_device;
	for (const auto &candidate_physical_device : myvk::PhysicalDevice::Fetch(instance)) {
		if (candidate_physical_device->GetExtensionSupport(VK_NV_COOPERATIVE_MATRIX_EXTENSION_NAME)) {
			physical_device = candidate_physical_device;
			break;
		}
	}
	if (physical_device == nullptr) {
		spdlog::error("No Device Supporting VK_NV_cooperative_matrix Found");
		return -1;
	}
	spdlog::info("Physical Device: {}", physical_device->GetProperties().vk10.deviceName);
	auto features = physical_device->GetDefaultFeatures();
	features.vk11.storageBuffer16BitAccess = VK_TRUE;
	features.vk12.bufferDeviceAddress = VK_TRUE;
	features.vk12.hostQueryReset = VK_TRUE;
	features.vk12.shaderFloat16 = VK_TRUE;
	features.vk12.vulkanMemoryModel = VK_TRUE;
	features.vk12.vulkanMemoryModelDeviceScope = VK_TRUE;
	features.vk13.computeFullSubgroups = VK_TRUE;
	features.vk13.subgroupSizeControl = VK_TRUE;
	VkPhysicalDeviceShaderAtomicFloatFeaturesEXT atomic_float_features = {
	    .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_ATOMIC_FLOAT_FEATURES_EXT,
	    .shaderBufferFloat32AtomicAdd = VK_TRUE};
	VkPhysicalDeviceCooperativeMatrixFeaturesNV cooperative_matrix_features = {
	    .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_COOPERATIVE_MATRIX_FEATURES_NV,
	    .pNext = &atomic_float_features,
	    .cooperativeMatrix = VK_TRUE,
	    .cooperativeMatrixRobustBufferAccess = VK_FALSE};
	VkPhysicalDeviceRayQueryFeaturesKHR ray_query_features = {
	    .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_RAY_QUERY_FEATURES_KHR,
	    .pNext = &cooperative_matrix_features,
	    .rayQuery = VK_TRUE};
	VkPhysicalDeviceAccelerationStructureFeaturesKHR accel_features = {
	    .sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_ACCELERATION_STRUCTURE_FEATURES_KHR,
	    .pNext = &ray_query_features,
	    .accelerationStructure = VK_TRUE};
	features.vk13.pNext = &accel_features;
	auto device = myvk::Device::Create(
	    physical_device,
	    myvk::GenericPresentQueueSelector{&generic_queue, myvk::Surface::Create(instance, window), &present_queue},
	    features,
	    {VK_KHR_SWAPCHAIN_EXTENSION_NAME, VK_KHR_DEFERRED_HOST_OPERATIONS_EXTENSION_NAME,
	     VK_KHR_ACCELERATION_STRUCTURE_EXTENSION_NAME, VK_KHR_RAY_QUERY_EXTENSION_NAME,
	     VK_NV_COOPERATIVE_MATRIX_EXTENSION_NAME, VK_EXT_SHADER_ATOMIC_FLOAT_EXTENSION_NAME});
	myvk::ImGuiInit(window, myvk::CommandPool::Create(generic_queue));

	auto camera = myvk::MakePtr<Camera>();
	Camera::Control cam_control{.sensitivity = 0.005f, .speed = 0.5f, .prev_cursor_pos = {}};

	auto vk_scene = myvk::MakePtr<VkScene>(generic_queue, scene);
	auto vk_scene_blas = myvk::MakePtr<VkSceneBLAS>(vk_scene);
	auto vk_scene_tlas = myvk::MakePtr<VkSceneTLAS>(vk_scene_blas);
	auto vk_nrc_state = myvk::MakePtr<VkNRCState>(generic_queue, VkExtent2D{kWidth, kHeight});

	auto frame_manager = myvk::FrameManager::Create(generic_queue, present_queue, false, kFrameCount);
	frame_manager->SetResizeFunc([&](VkExtent2D extent) {
		vk_nrc_state->ResetAccumulateImage(extent);
		vk_nrc_state->ResetAccumulateCount();
	});
	std::array<myvk::Ptr<rg::NRCRenderGraph>, kFrameCount> render_graphs;
	for (auto &rg : render_graphs)
		rg = myvk::MakePtr<rg::NRCRenderGraph>(frame_manager, vk_scene_tlas, vk_nrc_state, camera);

	bool view_accumulate = vk_nrc_state->IsAccumulate();
	int view_left_method = static_cast<int>(vk_nrc_state->GetLeftMethod());
	int view_right_method = static_cast<int>(vk_nrc_state->GetRightMethod());
	bool nrc_use_ema = vk_nrc_state->IsUseEMAWeights(), nrc_lock = false, nrc_train_one_frame = false;

	double prev_time = glfwGetTime();
	while (!glfwWindowShouldClose(window)) {
		double delta;
		{
			double cur_time = glfwGetTime();
			delta = cur_time - prev_time;
			prev_time = cur_time;
		}
		glfwPollEvents();

		nrc_train_one_frame = false;

		myvk::ImGuiNewFrame();
		ImGui::Begin("Panel");
		ImGui::Text("FPS %.1f", ImGui::GetIO().Framerate);
		if (ImGui::CollapsingHeader("View")) {
			if (ImGui::Checkbox("Accumulate", &view_accumulate))
				vk_nrc_state->SetAccumulate(view_accumulate);
			if (vk_nrc_state->IsAccumulate()) {
				ImGui::SameLine();
				ImGui::Text("SPP %d", vk_nrc_state->GetAccumulateCount());
			}
			constexpr const char *kViewTypeComboStr = "None\0NRC\0Cache\0";
			if (ImGui::Combo("Left", &view_left_method, kViewTypeComboStr)) {
				vk_nrc_state->SetLeftMethod(static_cast<VkNRCState::Method>(view_left_method));
				vk_nrc_state->ResetAccumulateCount();
			}
			if (ImGui::Combo("Right", &view_right_method, kViewTypeComboStr)) {
				vk_nrc_state->SetRightMethod(static_cast<VkNRCState::Method>(view_right_method));
				vk_nrc_state->ResetAccumulateCount();
			}
		}
		if (ImGui::CollapsingHeader("NRC")) {
			if (ImGui::Checkbox("Use EMA", &nrc_use_ema)) {
				vk_nrc_state->SetUseEMAWeights(nrc_use_ema);
				vk_nrc_state->ResetAccumulateCount();
			}
			if (ImGui::Checkbox("Lock", &nrc_lock))
				vk_nrc_state->ResetAccumulateCount();
			if (nrc_lock) {
				ImGui::SameLine();
				if (ImGui::Button("Train 1-Frame")) {
					vk_nrc_state->ResetAccumulateCount();
					nrc_train_one_frame = true;
				}
			}
			if (ImGui::Button("Re-Train")) {
				vk_nrc_state->ResetAccumulateCount();
				vk_nrc_state->ResetMLPBuffers();
			}
		}
		ImGui::End();
		ImGui::Render();

		if (camera->DragControl(window, &cam_control, delta))
			vk_nrc_state->ResetAccumulateCount();

		if (nrc_lock && !nrc_train_one_frame)
			vk_nrc_state->SetTrainProbability(0.0f);
		else
			vk_nrc_state->SetTrainProbability(VkNRCState::GetDefaultTrainProbability());

		if (frame_manager->NewFrame()) {
			const auto &command_buffer = frame_manager->GetCurrentCommandBuffer();
			auto &render_graph = render_graphs[frame_manager->GetCurrentFrame()];

			command_buffer->Begin(VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT);
			render_graph->SetCanvasSize(frame_manager->GetExtent());
			render_graph->CmdExecute(command_buffer);
			command_buffer->End();

			frame_manager->Render();
		}

		vk_nrc_state->NextFrame();
	}

	frame_manager->WaitIdle();
	glfwTerminate();
	return 0;
}
