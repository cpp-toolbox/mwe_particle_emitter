#include "graphics/window/window.hpp"

int main() {
    unsigned int width_px = 800;
    unsigned int height_px = 800;
    bool fullscreen = false;

    GLFWwindow *window =
        initialize_glfw_glad_and_return_window(width_px, height_px, "mwe particle emitter", fullscreen, false, false);

    while (not glfwWindowShouldClose(window)) {
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}
