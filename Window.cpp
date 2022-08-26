#include "Window.hpp"
// #include <glad/glad.h> included by Window.hpp
// #include <GLFW/glfw3.h>
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wgnu-anonymous-struct"
// #pragma GCC diagnostic ignored "-Wnested-anon-types"
// #include <glm/vec4.hpp>
// #include <glm/mat4x4.hpp>
// #include <glm/gtc/matrix_transform.hpp>
// #include <glm/gtc/type_ptr.hpp>
// #pragma GCC diagnostic pop
#include <exception>
#include <stdexcept>

#include <iostream>

namespace GLFW {
    // Classe Session
    Session::Session() : initialized_glad {false} {
        if (!glfwInit()) {
            throw std::runtime_error("Error initializing GLFW session\n");
        }
        static_assert(GLFW_TRUE == 1, "GLFW_TRUE diferente de 1");
        static_assert(GLFW_FALSE == 0, "GLFW_FALSE diferente de 0");
    }

    Session::~Session() {
        glfwTerminate();
    }

    void Session::pollEvents() {
        glfwPollEvents();
    }
    // Fim da classe Session
    //----------------------

    void framebuffer_size_callback(GLFWwindow* window, int width, int height);

    Window::Window(Session& session, const char* title, options opts) :
    session_m{session} {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, opts.version_major);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, opts.version_minor);
        glfwWindowHint(GLFW_DECORATED, opts.decorated);
        glfwWindowHint(GLFW_RESIZABLE, opts.resizable);
        glfwWindowHint(GLFW_OPENGL_PROFILE, opts.profile);
        glfwWindowHint(GLFW_SAMPLES, opts.samples);

        int width = opts.window_width;
        int height = opts.window_height;
        this->window = {
            glfwCreateWindow(width, height, title, NULL, NULL)
            };
        if (!window) {
            throw std::runtime_error { "Error creating window\n" };
        }
        glfwMakeContextCurrent(window);

        if (!session_m.initialized_glad) {
            if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
                throw std::runtime_error { "Error loading GLAD\n"};
            }
            session_m.initialized_glad = true;
        }

        glViewport(0, 0, width, height);
        
        this->arrow = glfwCreateStandardCursor(GLFW_ARROW_CURSOR);
        this->hand = glfwCreateStandardCursor(GLFW_HAND_CURSOR);

        glfwSetWindowUserPointer(this->window, this);
        glfwSetFramebufferSizeCallback(this->window, framebuffer_size_callback);
        glfwMakeContextCurrent(nullptr);
    }

    /*STUDY THAT::::::::::::

    class MyGlWindow {
    public:
        std::function<void(MyGlWindow*)> onClose;
        std::function<void(MyGlWindow*, int, int, int)> onMouseClick = [](auto self, int, int, int) { };
    };

    void makeWindow() {
        GLFWwindow* glfwWindow;
        MyGlWindow* myWindow;

        // ... Initialize everything here ...

        glfwSetWindowUserPointer(glfwWindow, myWindow);

        #define genericCallback(functionName)
            [](GLFWwindow* window, auto... args) {
                auto pointer = static_cast<MyGlWindow*>(glfwGetWindowUserPointer(window));
                if (pointer->functionName) pointer->functionName(pointer, args...);
            }

        glfwSetWindowCloseCallback(glfwWindow, genericCallback(onClose));
        glfwSetMouseButtonCallback(glfwWindow, genericCallback(onMouseClick));
        #undef genericCallback

        myWindow->onMouseClick = [](auto self, int, int, int) {
            std::cout << "I'm such a rebel" << std::endl;
            self->onClose = [](auto self) {
                std::cout << "I'm such a rebellion" << std::endl;
            };
        };
    }*/

    Window::~Window() {
        // não é necessário chamar outras funções, já que os recursos serão
        // liberados quando glfwTerminate for chamada
        glfwDestroyWindow(this->window);
    }

    ////////////////////////////////////////////////////////////////////////////////

    // void Window::swapBuffers() {
    //     glfwSwapBuffers(this->window);
    // }

    void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
        glViewport(0, 0, width, height);
    }
}