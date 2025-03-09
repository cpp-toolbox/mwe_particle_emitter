#include "graphics/batcher_x_particle_emitter_logic/batcher_x_particle_emitter_logic.hpp"
#include "utility/unique_id_generator/unique_id_generator.hpp"
#include <random>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image.h>
#include <stb_image_write.h>

#include "graphics/particle_emitter/particle_emitter.hpp"
#include "graphics/vertex_geometry/vertex_geometry.hpp"
#include "graphics/texture_packer/texture_packer.hpp"
#include "graphics/batcher/generated/batcher.hpp"
#include "graphics/shader_cache/shader_cache.hpp"
#include "graphics/fps_camera/fps_camera.hpp"
#include "graphics/draw_info/draw_info.hpp"
#include "graphics/window/window.hpp"

#include "input/glfw_lambda_callback_manager/glfw_lambda_callback_manager.hpp"
#include "input/input_state/input_state.hpp"

#include "utility/texture_loader/texture_loader.hpp"
#include "utility/resource_path/resource_path.hpp"
#include "utility/print_utils/print_utils.hpp"

#include <filesystem>

class SmokeParticleEmitter {
  public:
    ParticleEmitter particle_emitter;

    SmokeParticleEmitter(std::function<void(int)> on_particle_spawn_callback,
                         std::function<void(int)> on_particle_death_callback)
        : particle_emitter(life_span_lambda(), initial_velocity_lambda(), velocity_change_lambda(), scaling_lambda(),
                           rotation_lambda(), spawn_delay_lambda(), on_particle_spawn_callback,
                           on_particle_death_callback) {}

  private:
    static std::function<float()> life_span_lambda() {
        return []() -> float {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> dist(1.0f, 3.0f);
            return dist(rng);
        };
    }

    static std::function<glm::vec3()> initial_velocity_lambda() {
        return []() -> glm::vec3 {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> horizontal_dist(-0.5f, 0.5f); // Small lateral variance
            std::uniform_real_distribution<float> upward_dist(1.0f, 2.0f);      // Strong upward push

            // Initial upward push with slight lateral drift
            float dx = horizontal_dist(rng);
            float dy = upward_dist(rng); // Main upward velocity
            float dz = horizontal_dist(rng);

            return glm::vec3(dx, dy, dz);
            /*return glm::vec3(0, 0, 0);*/
        };
    }

    static std::function<glm::vec3(float, float)> velocity_change_lambda() {
        return [](float life_percentage, float delta_time) -> glm::vec3 {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> horizontal_dist(-0.010f, 0.010f); // Small lateral variance
            std::uniform_real_distribution<float> vertical_dist(0.2f, 0.3f);

            float accel_x = horizontal_dist(rng);
            float accel_y = vertical_dist(rng);
            float accel_z = horizontal_dist(rng);

            glm::vec3 smoke_push_down = -glm::vec3(accel_x, accel_y, accel_z) * delta_time;
            return smoke_push_down;
        };
    }

    static std::function<float(float)> scaling_lambda() {
        return [](float life_percentage) -> float { return std::max(life_percentage * 2.0f, 0.0f); };
    }

    static std::function<float(float)> rotation_lambda() {
        return [](float life_percentage) -> float { return life_percentage / 5.0f; };
    }

    static std::function<float()> spawn_delay_lambda() {
        return []() -> float {
            return 0.01f; // Spawn a new particle every 0.05 seconds
        };
    }
};

int main() {
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::debug);
    std::vector<spdlog::sink_ptr> sinks = {console_sink};

    unsigned int width_px = 800;
    unsigned int height_px = 800;
    bool fullscreen = false;

    Window window;
    window.initialize_glfw_glad_and_return_window(width_px, height_px, "mwe particle emitter", fullscreen, false,
                                                  false);

    InputState input_state;

    FPSCamera fps_camera;

    std::function<void(unsigned int)> char_callback = [](unsigned int codepoint) {};
    std::function<void(int, int, int, int)> key_callback = [&](int key, int scancode, int action, int mods) {
        input_state.glfw_key_callback(key, scancode, action, mods);
    };
    std::function<void(double, double)> mouse_pos_callback = [&](double xpos, double ypos) {
        fps_camera.mouse_callback(xpos, ypos);
    };
    std::function<void(int, int, int)> mouse_button_callback = [](int button, int action, int mods) {};
    std::function<void(int, int)> frame_buffer_size_callback = [](int width, int height) {};
    GLFWLambdaCallbackManager glcm(window.glfw_window, char_callback, key_callback, mouse_pos_callback,
                                   mouse_button_callback, frame_buffer_size_callback);

    glfwSetInputMode(window.glfw_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED); // Hide and capture the mouse

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    auto requested_shaders = {ShaderType::TEXTURE_PACKER_CWL_V_TRANSFORMATION_UBOS_1024};
    ShaderCache shader_cache(requested_shaders, sinks);
    Batcher batcher(shader_cache);

    ResourcePath rp(false);

    const auto textures_directory = std::filesystem::path("assets");
    std::filesystem::path output_dir = std::filesystem::path("assets") / "packed_textures";
    int container_side_length = 1024;

    TexturePacker texture_packer(textures_directory, output_dir, container_side_length);
    shader_cache.set_uniform(ShaderType::TEXTURE_PACKER_CWL_V_TRANSFORMATION_UBOS_1024,
                             ShaderUniformVariable::PACKED_TEXTURE_BOUNDING_BOXES, 1);

    BatcherXParticleEmitterLogic batcher_x_particle_emitter_logic(
        texture_packer, batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.ltw_object_id_generator,
        batcher);

    std::function<void(int)> on_particle_spawn = [&](int particle_id) {
        batcher_x_particle_emitter_logic.on_particle_spawn(particle_id);
    };

    std::function<void(int)> on_particle_death = [&](int particle_id) {
        batcher_x_particle_emitter_logic.on_particle_death(particle_id);
    };

    SmokeParticleEmitter spe(on_particle_spawn, on_particle_death);

    double previous_time = glfwGetTime();
    while (not glfwWindowShouldClose(window.glfw_window)) {
        double current_time = glfwGetTime();
        double delta_time = current_time - previous_time;
        previous_time = current_time;

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        fps_camera.process_input(input_state.is_pressed(EKey::LEFT_CONTROL), input_state.is_pressed(EKey::TAB),
                                 input_state.is_pressed(EKey::w), input_state.is_pressed(EKey::a),
                                 input_state.is_pressed(EKey::s), input_state.is_pressed(EKey::d),
                                 input_state.is_pressed(EKey::SPACE), input_state.is_pressed(EKey::LEFT_SHIFT),
                                 delta_time);

        glm::mat4 projection = fps_camera.get_projection_matrix(width_px, height_px);
        glm::mat4 view = fps_camera.get_view_matrix();

        shader_cache.set_uniform(ShaderType::TEXTURE_PACKER_CWL_V_TRANSFORMATION_UBOS_1024,
                                 ShaderUniformVariable::CAMERA_TO_CLIP, projection);
        shader_cache.set_uniform(ShaderType::TEXTURE_PACKER_CWL_V_TRANSFORMATION_UBOS_1024,
                                 ShaderUniformVariable::WORLD_TO_CAMERA, view);

        batcher_x_particle_emitter_logic.draw_particles(spe.particle_emitter, fps_camera, delta_time, projection, view);

        // load in the matrices
        batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.upload_ltw_matrices();
        batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.draw_everything();

        glfwSwapBuffers(window.glfw_window);
        glfwPollEvents();
        TemporalBinarySignal::process_all();
    }
}
