#ifndef WINDOW_HPP
#define WINDOW_HPP

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#pragma GCC diagnostic push
//needed to gcc
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wpedantic"
//needed to clang
#pragma GCC diagnostic ignored "-Wgnu-anonymous-struct"
#pragma GCC diagnostic ignored "-Wnested-anon-types"
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#pragma GCC diagnostic pop

namespace GLFW {
    class Session {
    private:
        bool initialized_glad;
        friend class Window;
    public:
        Session();
        ~Session();
        void pollEvents();
    };

    class Window {
    private:
        Session& session_m;
        GLFWwindow* window;
        
        GLFWcursor* arrow { nullptr };
        GLFWcursor* hand { nullptr };
    public:
        struct options {
        public:
            int window_width { 800 };
            int window_height { 600 };
            int version_major { 4 };
            int version_minor { 5 };
            int decorated { 0 };
            int resizable { 0 };
            int profile { GLFW_OPENGL_CORE_PROFILE };
            int samples { 0 };
        };
        Window(Session& session, const char* title, options opts);
        ~Window();

        GLFWwindow* justGimmeTheWindow() { return window; }

        // bool shouldClose() { return glfwWindowShouldClose(this->window); }
        // void makeContextCurrent() { glfwMakeContextCurrent(this->window); }
        // void processInput();

        // void swapBuffers();
    };
}

#endif // WINDOW_HPP