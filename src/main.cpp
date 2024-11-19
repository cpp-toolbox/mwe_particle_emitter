#include <random>
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#include "graphics/batcher/generated/batcher.hpp"
#include "graphics/fps_camera/fps_camera.hpp"
#include "graphics/shader_cache/shader_cache.hpp"
#include "graphics/window/window.hpp"
#include "graphics/vertex_geometry/vertex_geometry.hpp"
#include "input/glfw_lambda_callback_manager/glfw_lambda_callback_manager.hpp"
#include "utility/texture_loader/texture_loader.hpp"
#include "graphics/particle_emitter/particle_emitter.hpp"
#include "graphics/texture_packer/texture_packer.hpp"

class SmokeParticleEmitter {
  public:
    ParticleEmitter particle_emitter;

    SmokeParticleEmitter(unsigned int max_particles)
        : particle_emitter(life_span_lambda(), initial_velocity_lambda(), velocity_change_lambda(), scaling_lambda(),
                           rotation_lambda(), spawn_delay_lambda(), max_particles) {}

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

GLuint create_texture_from_data(const TextureData &texture_data) {
    GLuint texture_id;
    glGenTextures(1, &texture_id);
    glBindTexture(GL_TEXTURE_2D, texture_id);

    // Load the texture data into OpenGL
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texture_data.width, texture_data.height, 0,
                 texture_data.num_components == 4 ? GL_RGBA : GL_RGB, GL_UNSIGNED_BYTE, texture_data.image_data.data());

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glBindTexture(GL_TEXTURE_2D, 0); // Unbind the texture

    return texture_id;
}

// Wrapper that automatically creates a lambda for member functions
template <typename T, typename R, typename... Args> auto wrap_member_function(T &obj, R (T::*f)(Args...)) {
    // Return a std::function that wraps the member function in a lambda
    return std::function<R(Args...)>{[&obj, f](Args &&...args) { return (obj.*f)(std::forward<Args>(args)...); }};
}

int main() {
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::debug);
    std::vector<spdlog::sink_ptr> sinks = {console_sink};

    unsigned int width_px = 800;
    unsigned int height_px = 800;
    bool fullscreen = false;

    GLFWwindow *window =
        initialize_glfw_glad_and_return_window(width_px, height_px, "mwe particle emitter", fullscreen, false, false);

    FPSCamera camera(glm::vec3(0, 0, 3), 50, width_px, height_px, 90, 0.1, 50);
    std::function<void(unsigned int)> char_callback = [](unsigned int _) {};
    std::function<void(int, int, int, int)> key_callback = [](int _, int _1, int _2, int _3) {};
    std::function<void(double, double)> mouse_pos_callback = wrap_member_function(camera, &FPSCamera::mouse_callback);
    std::function<void(int, int, int)> mouse_button_callback = [](int _, int _1, int _2) {};
    GLFWLambdaCallbackManager glcm(window, char_callback, key_callback, mouse_pos_callback, mouse_button_callback);

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED); // Hide and capture the mouse

    TextureData texture_data = load_texture_from_file("assets/smoke_64px.png");
    GLuint texture_id = create_texture_from_data(texture_data);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    /*glBlendFunc(GL_SRC_ALPHA, GL_ONE);*/
    glBindTexture(GL_TEXTURE_2D, texture_id);

    auto requested_shaders = {ShaderType::TEXTURE_PACKER_CWL_V_TRANSFORMATION_UBOS_1024};
    ShaderCache shader_cache(requested_shaders, sinks);
    Batcher batcher(shader_cache);

    TexturePacker texture_packer("assets/packed_textures/packed_texture.json",
                                 {"assets/packed_textures/packed_texture_0.png"});

    auto vertices = generate_square_vertices(0, 0, 0.5);
    auto indices = generate_rectangle_indices();
    std::vector<glm::vec2> local_uvs = generate_rectangle_texture_coordinates();
    auto texture_coordinates = texture_packer.get_packed_texture_coordinates("assets/smoke_64px.png", local_uvs);
    auto pt_idx = texture_packer.get_packed_texture_index_of_texture("assets/smoke_64px.png");
    std::vector<int> pt_idxs(4, pt_idx); // 4 because square

    GLuint ltw_matrices_gl_name;
    glm::mat4 ltw_matrices[1024];

    // initialize all matrices to identity matrices
    for (int i = 0; i < 1024; ++i) {
        ltw_matrices[i] = glm::mat4(1.0f);
    }

    glGenBuffers(1, &ltw_matrices_gl_name);
    glBindBuffer(GL_UNIFORM_BUFFER, ltw_matrices_gl_name);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(ltw_matrices), ltw_matrices, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ltw_matrices_gl_name);

    SmokeParticleEmitter spe(1000);

    double previous_time = glfwGetTime();
    while (not glfwWindowShouldClose(window)) {
        double current_time = glfwGetTime();
        double delta_time = current_time - previous_time;
        previous_time = current_time;

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        camera.process_input(window, delta_time);

        glm::mat4 projection = camera.get_projection_matrix();
        glm::mat4 view = camera.get_view_matrix();

        shader_cache.set_uniform(ShaderType::TEXTURE_PACKER_CWL_V_TRANSFORMATION_UBOS_1024,
                                 ShaderUniformVariable::CAMERA_TO_CLIP, projection);
        shader_cache.set_uniform(ShaderType::TEXTURE_PACKER_CWL_V_TRANSFORMATION_UBOS_1024,
                                 ShaderUniformVariable::WORLD_TO_CAMERA, view);

        spe.particle_emitter.update(delta_time, projection * view);
        auto particles = spe.particle_emitter.get_particles_sorted_by_distance();

        for (size_t i = 0; i < particles.size(); ++i) {
            auto &curr_particle = particles[i];

            //  compute the up vector (assuming we want it to be along the y-axis)
            glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
            glm::vec3 forward = camera.transform.compute_forward_vector();

            glm::vec3 right = glm::normalize(glm::cross(up, forward));

            up = glm::normalize(glm::cross(forward, right));

            // this makes it billboarded
            glm::mat4 rotation_matrix = glm::mat4(1.0f);
            rotation_matrix[0] = glm::vec4(right, 0.0f);
            rotation_matrix[1] = glm::vec4(up, 0.0f);
            rotation_matrix[2] = glm::vec4(-forward, 0.0f); // We negate the direction for correct facing

            glm::mat4 transform = glm::translate(glm::mat4(1.0f), curr_particle.transform.position);
            transform *= rotation_matrix;
            transform = glm::scale(transform, curr_particle.transform.scale);

            ltw_matrices[i] = transform;

            if (curr_particle.is_alive()) {

                /*auto nv =*/
                /*    generate_rectangle_vertices_3d(particle.transform.position,
                 * camera.transform.compute_right_vector(),*/
                /*                                   camera.transform.compute_up_vector(), 1, 1);*/

                std::vector<unsigned int> ltw_mat_idxs(4, i);
                batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.queue_draw(
                    indices, vertices, ltw_mat_idxs, pt_idxs, texture_coordinates);
            }
        }

        // load in the matrices
        glBindBuffer(GL_UNIFORM_BUFFER, ltw_matrices_gl_name);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(ltw_matrices), ltw_matrices);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.draw_everything();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}
