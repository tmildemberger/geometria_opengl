#define LER 0
#if LER==1

#include <iostream>
#include <array>
#include <cstdlib>
#include <cmath>
#include <glad/glad.h>
// #include <GLFW/glfw3.h>
#include "shader.hpp"
#include "Window.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wgnu-anonymous-struct"
#pragma GCC diagnostic ignored "-Wnested-anon-types"
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#pragma GCC diagnostic pop

extern "C" {
    #define DLL_EXPORT [[gnu::dllexport]]
    DLL_EXPORT unsigned long NvOptimusEnablement = 0x00000001;
    DLL_EXPORT int AmdPowerXpressRequestHighPerformance = 1;
}

// struct vec2 {
//     double x;
//     double y;
// };

// void framebuffer_size_callback(GLFWwindow* window, int width, int height);
// void cursor_enter_callback(GLFWwindow* window, int entered);
// void mouse_button_callback(GLFWwindow* window, int button,
//                            int action, int mods);
// void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
// double area(glm::dvec2 pA, glm::dvec2 pB, glm::dvec2 pC);
// bool is_inside(GLFWwindow* window, glm::dvec2 points[3]);

const unsigned int screen_width = 800;
const unsigned int screen_height = 600;

// double area(glm::dvec2 pA, glm::dvec2 pB, glm::dvec2 pC) {
//     return std::abs((pA.y-pB.y)*pC.x + (pB.x-pA.x)*pC.y + pA.x*pB.y - pA.y*pB.x)
//     ;
// }

// bool is_inside(GLFWwindow* window, glm::dvec2 points[3]) {
//     glm::dvec2 cursor;
//     glfwGetCursorPos(window, &cursor.x, &cursor.y);
//     double possibleAreas { 0.0 };
//     // for (int i { 0 }; i < 3; ++i) {
//     //     possibleAreas += area(points[i], points[(i+1)%3], cursor);
//     // }
//     glm::dmat3 pointsMat { glm::dvec3 { points[0], 1.0 },
//                            glm::dvec3 { points[1], 1.0 },
//                            glm::dvec3 { points[2], 1.0 } };
//     for (int i { 0 }; i < 3; ++i) {
//         glm::dmat3 possibleMat { pointsMat };
//         possibleMat[i] = glm::dvec3 { cursor, 1.0 };
//         possibleAreas += std::abs(glm::determinant(possibleMat));
//     }
//     // return possibleAreas <= area(points[0], points[1], points[2]);
//     return possibleAreas <= std::abs(glm::determinant(pointsMat));
// }

#include "stb_img.hpp"
#include "texture.hpp"
#include "drawing.hpp"
#include "renderer.hpp"

void message_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, GLchar const* message, void const* user_param);

int mwain() {

    //desired:
    //GLFW::SESSION
    //GLFW::WINDOW??
    // if (!glfwInit()) {
    //     return 1;
    // }
    GLFW::Session session {};

    GLFW::Window::options opts {};
    opts.samples = 128;
    GLFW::Window win {"LearnOpenGL", opts};

    win.makeContextCurrent();

    // GLFW::Window::options opts {}; //all but samples to default
    // //OpenGL version 4.5, no title bar, core profile is the default
    // opts.samples = 128;
    // GLFW::Window win {screen_width, screen_height, "LearnOpenGL", opts};

    // opts.decorated = 1;
    // GLFW::Window winK {screen_width /4, screen_height /4, "LearnClasses", opts};

    // opts.decorated = 0;

    // glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    // glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    // glfwWindowHint(GLFW_DECORATED, 0);
    // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // glfwWindowHint(GLFW_SAMPLES, 128);
    // GLFWwindow* window {
    //     glfwCreateWindow(screen_width, screen_height, "LearnOpenGL", NULL, NULL)
    //     };
    // if (!window) {
    //     std::cout << "Failed to create GLFW window\n";
    //     glfwTerminate();
    //     return 2; //change this all today
    // }
    // glfwMakeContextCurrent(window);

    // if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
    //     std::cout << "Failed to initialize GLAD\n";
    //     return 3;
    // }
    // glViewport(0, 0, screen_width, screen_height);
    // glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // glfwSetCursorEnterCallback(window, cursor_enter_callback);
    // glfwSetMouseButtonCallback(window, mouse_button_callback);
    // glfwSetCursorPosCallback(window, cursor_position_callback);
    // glfwSetScrollCallback(window, scroll_callback);

    Renderer renderer;

    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    glDebugMessageCallback(message_callback, nullptr);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);

    Shader program {"shaders/vertex.glsl", "shaders/fragment.glsl"};

    // float vertices[] {
    //     -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,
    //      0.5f, -0.5f, -0.5f,  1.0f, 0.0f,
    //      0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
    //      0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
    //     -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
    //     -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,

    //     -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
    //      0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
    //      0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
    //      0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
    //     -0.5f,  0.5f,  0.5f,  0.0f, 1.0f,
    //     -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,

    //     -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
    //     -0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
    //     -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
    //     -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
    //     -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
    //     -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

    //      0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
    //      0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
    //      0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
    //      0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
    //      0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
    //      0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

    //     -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
    //      0.5f, -0.5f, -0.5f,  1.0f, 1.0f,
    //      0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
    //      0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
    //     -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
    //     -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,

    //     -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
    //      0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
    //      0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
    //      0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
    //     -0.5f,  0.5f,  0.5f,  0.0f, 0.0f,
    //     -0.5f,  0.5f, -0.5f,  0.0f, 1.0f
    // };
    // float vertices[] {
    //     //  0.5f,  0.5f, 0.0f, // top right
    //     //  0.5f, -0.5f, 0.0f, // bottom right
    //     // -0.5f, -0.5f, 0.0f, // bottom left
    //     // -0.5f,  0.5f, 0.0f  // top left
    //     // -0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f,
    //     //  0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f,
    //     //  0.0f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
    //     // -0.8f,  0.8f, 0.0f, 1.0f, 0.5f, 0.2f,
    //     // -0.5f,  0.2f, 0.0f, 0.8f, 0.8f, 0.0f,
    //     // -0.2f,  0.5f, 0.0f, 1.0f, 0.5f, 0.2f,
    //     // -0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f,
    //     //  0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f,
    //     //  0.0f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
    //     // -0.8f,  0.8f, 0.0f, 1.0f, 0.0f, 0.0f,
    //     // -0.5f,  0.2f, 0.0f, 0.0f, 0.0f, 1.0f,
    //     // -0.2f,  0.5f, 0.0f, 0.0f, 1.0f, 0.0f,
    //     // // positions         // colors          // texture coords
    //     //  0.5f,  0.5f, 0.0f,  1.0f, 0.0f, 0.0f,  1.0f, 1.0f, // top right
    //     //  0.5f, -0.5f, 0.0f,  0.0f, 1.0f, 0.0f,  1.0f, 0.0f, // bottom right
    //     // -0.5f, -0.5f, 0.0f,  0.0f, 0.0f, 1.0f,  0.0f, 0.0f, // bottom left
    //     // -0.5f,  0.5f, 0.0f,  1.0f, 1.0f, 0.0f,  0.0f, 1.0f, // top left
    //      0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
    //      0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
    //     -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
    //     -0.5f,  0.5f,  0.5f,  0.0f, 1.0f,

    //      0.5f,  0.5f, -0.5f,  1.0f, 0.0f,
    //      0.5f, -0.5f, -0.5f,  1.0f, 1.0f,
    //     -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
    //     -0.5f,  0.5f, -0.5f,  0.0f, 0.0f,
    // };

    // unsigned int indices[] {
    //     // 0, 1, 2,    // first triangle
    //     // 0, 2, 3,    // second triangle
    // //     0, 1, 2,    // first triangle
    // //     3, 4, 5,    // second triangle
    // //     0, 4, 5,
    // //     5, 0, 2,
    //     0, 1, 2,
    //     0, 2, 3,
    //     4, 5, 1,
    //     4, 1, 0,
    //     7, 6, 5,
    //     7, 5, 4,
    //     3, 2, 6,
    //     3, 6, 7,
    //     4, 0, 3,
    //     4, 3, 7,
    //     1, 5, 6,
    //     1, 6, 2,
    // };

    // unsigned int EBO { 0 };
    // glGenBuffers(1, &EBO);

    // unsigned int VBO { 0 };
    // glCreateBuffers(1, &VBO);
    // glNamedBufferStorage(VBO, sizeof vertices, vertices, GL_DYNAMIC_STORAGE_BIT)
    // ;
    std::array vertices {
        -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,
         0.5f, -0.5f, -0.5f,  1.0f, 0.0f,
         0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
         0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, 0.0f,

        -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 1.0f,
        -0.5f,  0.5f,  0.5f,  0.0f, 1.0f,
        -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,

        -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
        -0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
        -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
        -0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

         0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
         0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
         0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
         0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
         0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 0.0f,

        -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,
         0.5f, -0.5f, -0.5f,  1.0f, 1.0f,
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f,
        -0.5f, -0.5f,  0.5f,  0.0f, 0.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, 1.0f,

        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f,
         0.5f,  0.5f, -0.5f,  1.0f, 1.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 0.0f,
        -0.5f,  0.5f,  0.5f,  0.0f, 0.0f,
        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f
    };

    //desired:
    //vertices and texcoords array definition
    //layout definition? (you won cherno)
    //Drawing cube {attributes, layout};
    // cube.render(renderer);
    BufferLayout layout;
    layout.push<float>(3);
    layout.push<float>(2);
    Drawing cube_model {vertices, layout};
    std::array<Drawing, 10> cubes {cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,};

    // unsigned int VAO { 0 };
    // glCreateVertexArrays(1, &VAO);

    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof indices, indices,
    //              GL_STATIC_DRAW);
    // glVertexArrayVertexBuffer(VAO, 0, VBO, 0, 5 * sizeof (float));

    // glEnableVertexArrayAttrib(VAO, 0);
    // glEnableVertexArrayAttrib(VAO, 1);

    // glVertexArrayAttribFormat(VAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
    // glVertexArrayAttribFormat(VAO, 1, 2, GL_FLOAT, GL_FALSE, 3 * sizeof (float));

    // glVertexArrayAttribBinding(VAO, 0, 0);
    // glVertexArrayAttribBinding(VAO, 1, 0);

    Texture texture0 {"wooden_container.png", GL_BGR};
    Texture texture1 {"awesomeface.png", GL_RGBA};

    program.setUniform("texture0", 0);
    program.setUniform("texture1", 1);
    texture0.useAsUnit(0);
    texture1.useAsUnit(1);
    glm::vec3 cubePositions[] {
        glm::vec3( 0.0f,  0.0f,  0.0f),
        glm::vec3( 2.0f,  5.0f, -15.0f),
        glm::vec3(-1.5f, -2.2f, -2.5f),
        glm::vec3(-3.8f, -2.0f, -12.3f),
        glm::vec3( 2.4f, -0.4f, -3.5f),
        glm::vec3(-1.7f,  3.0f, -7.5f),
        glm::vec3( 1.3f, -2.0f, -2.5f),
        glm::vec3( 1.5f,  2.0f, -2.5f),
        glm::vec3( 1.5f,  0.2f, -1.5f),
        glm::vec3(-1.3f,  1.0f, -1.5f),
    };

    float mixing { 0.3 };
    float red { 0.7 };
    while (!win.shouldClose()) {
        win.processInput();

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        program.use();
        mixing = std::sin(static_cast<float>(glfwGetTime()) * 1.8 + 0.3);
        program.setUniform("mixing", mixing);
        // program.setFloat("mixing", mixing);
        red = std::sin(static_cast<float>(glfwGetTime()) * 0.4f) * 0.7f + 0.3f;
        program.setUniform("otherColorRed", red);

        glm::mat4 view { glm::mat4(1.0f) };
        view = glm::lookAt(
            win.cam.cameraPos,
            win.cam.cameraPos + win.cam.cameraFront,
            win.cam.cameraUp);

        program.setViewTransform(glm::value_ptr(view));

        glm::mat4 projection { glm::mat4(1.0f) };
        projection = glm::perspective(
            glm::radians(win.cam.fov),
            static_cast<float>(800) / static_cast<float>(600),
            0.1f,
            100.0f);
        program.setProjectionTransform(glm::value_ptr(projection));

        // glBindVertexArray(VAO);

        for (int i { 0 }; i < 10; ++i) {
            glm::mat4 model { glm::mat4(1.0f) };
            model = glm::translate(model, cubePositions[i]);
            float angle {
                20.f * i +
                2 * static_cast<float>(glfwGetTime()) * (i % 2 + 1)};
            model = glm::rotate(
                model,
                glm::radians(angle),
                glm::vec3 {1.0f, 0.3f, 0.5f});
            program.setModelTransform(glm::value_ptr(model));
            cubes[static_cast<unsigned>(i)].render(renderer);
        }

        // glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        // glfwSwapBuffers(window);
        win.swapBuffers();
        // glfwPollEvents();
        session.pollEvents();
    }

    // glDeleteVertexArrays(1, &VAO);
    // glDeleteBuffers(1, &EBO);
    // glDeleteBuffers(1, &VBO);

    //win is destroyed
    //session is destroyed
    // glfwDestroyWindow(window);
    // glfwTerminate();
    return 0;
}

void message_callback(GLenum source, GLenum type, GLuint id, GLenum severity,
                      GLsizei length, GLchar const* message, void const*) {
    const auto src_str = [source] () {
        switch (source)	{
            case GL_DEBUG_SOURCE_API: return "API";
            case GL_DEBUG_SOURCE_WINDOW_SYSTEM: return "WINDOW SYSTEM";
            case GL_DEBUG_SOURCE_SHADER_COMPILER: return "SHADER COMPILER";
            case GL_DEBUG_SOURCE_THIRD_PARTY: return "THIRD PARTY";
            case GL_DEBUG_SOURCE_APPLICATION: return "APPLICATION";
            case GL_DEBUG_SOURCE_OTHER: return "OTHER";
            default: return "IMPOSSIBLE";
        }
    }();

    const auto type_str = [type] () {
        switch (type) {
            case GL_DEBUG_TYPE_ERROR: return "ERROR";
            case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: return "DEPRECATED_BEHAVIOR"
            ;
            case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: return "UNDEFINED_BEHAVIOR";
            case GL_DEBUG_TYPE_PORTABILITY: return "PORTABILITY";
            case GL_DEBUG_TYPE_PERFORMANCE: return "PERFORMANCE";
            case GL_DEBUG_TYPE_MARKER: return "MARKER";
            case GL_DEBUG_TYPE_OTHER: return "OTHER";
            default: return "IMPOSSIBLE";
        }
    }();

    const auto severity_str = [severity] () {
        switch (severity) {
            case GL_DEBUG_SEVERITY_NOTIFICATION: return "NOTIFICATION";
            case GL_DEBUG_SEVERITY_LOW: return "LOW";
            case GL_DEBUG_SEVERITY_MEDIUM: return "MEDIUM";
            case GL_DEBUG_SEVERITY_HIGH: return "HIGH";
            default: return "IMPOSSIBLE";
        }
    }();
    if (severity_str == std::string{"HIGH"}) {
        std::cout << src_str << ", " << type_str << ", " <<
                    severity_str << ", " << id << ": ";
        throw std::runtime_error { message };
    }
    std::cout << src_str << ", " << type_str << ", " <<
                severity_str << ", " << id << ": " << message << '\n';
}

#endif