#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <numbers>
#include <limits>
#include <glad/glad.h>
// #include <GLFW/glfw3.h>
#include "shader.hpp"
#include "Window.hpp"
#include "stb_img.hpp"
#include "texture.hpp"
#include "drawing.hpp"
#include "renderer.hpp"

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

using Ponto = std::array<double, 2>;
double area_orientada(Ponto p1, Ponto p2, Ponto p3);
bool left(Ponto p1, Ponto p2, Ponto p3);
std::vector<Ponto> fecho_convexo(std::vector<Ponto> pontos);
double dist(Ponto p1, Ponto p2);
double angulo_interno(Ponto p1, Ponto p2, Ponto p3);

// área do paralelogramo, não do triângulo;
// produto vetorial do vetor p1p2 com o p1p3;
double area_orientada(Ponto p1, Ponto p2, Ponto p3) {
    // (x2 - x1)(y3 - y1) - (x3 - x1)(y2 - y1)
    return (p2[0] - p1[0])*(p3[1] - p1[1]) - (p3[0] - p1[0])*(p2[1] - p1[1]);
}

bool left(Ponto p1, Ponto p2, Ponto p3) {
    return area_orientada(p1, p2, p3) > 0.;
}

double dist(Ponto p1, Ponto p2) {
    return std::sqrt((p2[0] - p1[0])*(p2[0] - p1[0]) + (p2[1] - p1[1])*(p2[1] - p1[1]));
}

// ângulo entre p1p2 e p1p3
double angulo_interno(Ponto p1, Ponto p2, Ponto p3) {
    double seno = area_orientada(p1, p2, p3)/(dist(p1, p2)*dist(p1, p3));
    double angulo = std::asin(seno);
    if (angulo < 0.) {
        angulo += std::numbers::pi;
    }
    return angulo;
}

std::vector<Ponto> fecho_convexo(std::vector<Ponto> pontos) {
    // recebe 'pontos' como cópia mesmo, já que teremos que mudar a ordem dele
    std::vector<Ponto> fecho;
    std::size_t n = pontos.size();
    fecho.reserve(n);

    std::sort(pontos.begin(), pontos.end(), [](Ponto p1, Ponto p2) {
        // ordena pela coordenada x, do maior até o menor
        return p1[0] > p2[0];
    });
    
    // std::cout << "terminada ordenação, veja:" << std::endl;
    // for (std::size_t i = 0; i < pontos.size(); ++i) {
    //     std::cout << pontos[i][0] << ' ' << pontos[i][1] << std::endl;
    // }
    // std::cout << "----- Começo do calculo do fecho superior" << std::endl;

    // adiciona os dois pontos mais da esquerda para o fecho
    fecho.push_back(pontos[pontos.size() - 1]);
    fecho.push_back(pontos[pontos.size() - 2]);
    
    for (std::size_t i = pontos.size() - 2; i > 0; --i) {
        std::size_t j = fecho.size();
        while (j >= 2 && !left(pontos[i-1], fecho[j-1], fecho[j-2])) {
            --j;
            fecho.pop_back();
        }
        fecho.push_back(pontos[i-1]);
    }

    // std::cout << "fim do calculo do fecho superior, veja:" << std::endl;
    // for (std::size_t i = 0; i < fecho.size(); ++i) {
    //     std::cout << fecho[i][0] << ' ' << fecho[i][1] << std::endl;
    // }
    // std::cout << "-----começo do calculo do resto" << std::endl;


    std::vector<Ponto> fecho_inferior;
    fecho_inferior.reserve(n);

    // adiciona os dois pontos mais da direita para o fecho
    fecho_inferior.push_back(pontos[0]);
    fecho_inferior.push_back(pontos[1]);

    for (std::size_t i = 2; i < n; ++i) {
        std::size_t j = fecho_inferior.size();
        while (j >= 2 && !left(pontos[i], fecho_inferior[j-1], fecho_inferior[j-2])) {
            --j;
            fecho_inferior.pop_back();
        }
        fecho_inferior.push_back(pontos[i]);
    }

    // retira o último do fecho superior já que ele é o mesmo do início do
    // fecho inferior
    fecho.pop_back();

    fecho.insert(
        fecho.end(),
        std::make_move_iterator(fecho_inferior.begin()),
        std::make_move_iterator(fecho_inferior.end())
    );

    // retira o último do fecho já que ele é o mesmo do início
    fecho.pop_back();

    // vira ao contrário já que ele está na ordem horária por algum motivo
    std::reverse(fecho.begin(), fecho.end());

    return fecho;
}

using Reta = std::array<Ponto, 2>;

struct RetornoAlg {
    Ponto p;
    Reta r;
    double distancia;
};

double dist(Ponto p, Reta r);
RetornoAlg algoritmo(std::vector<Ponto> poligono);

double dist(Ponto p, Reta r) {
    double x3_x1 = p[0] - r[0][0];
    double x2_x1 = r[1][0] - r[0][0];
    double y3_y1 = p[1] - r[0][1];
    double y2_y1 = r[1][1] - r[0][1];
    double c = (x2_x1*x3_x1 + y2_y1*y3_y1) / (x2_x1*x2_x1 + y2_y1*y2_y1);
    double dist_x = x3_x1 - x2_x1*c;
    double dist_y = y3_y1 - y2_y1*c;
    return std::sqrt( dist_x*dist_x + dist_y*dist_y );
}

RetornoAlg algoritmo(std::vector<Ponto> poligono) {
    std::vector<Ponto> fecho = fecho_convexo(poligono);
    
    // n é o número de pontos no fecho convexo
    std::size_t n = fecho.size();
    std::vector<double> angulos(n, 0.);

    // adiciona manualmente primeiro e último ângulo
    // angulos[0] = angulo_interno(fecho[n-1], fecho[0], fecho[1]);
    // angulos[n-1] = angulo_interno(fecho[n-2], fecho[n-1], fecho[0]);
    angulos[0] = angulo_interno(fecho[0], fecho[n-1], fecho[1]);
    angulos[n-1] = angulo_interno(fecho[n-1], fecho[n-2], fecho[0]);

    for (std::size_t i = 1; i < n-1; ++i) {
        angulos[i] = angulo_interno(fecho[i], fecho[i-1], fecho[i+1]);
    }

    // pensar como ir somando os ângulos para poder descobrir
    // o somatório entre quaisquer pontos com uma só subtração
    std::vector<double> angulos_acumulados(n+1, 0.);

    for (std::size_t i = 1; i <= n; ++i) {
        angulos_acumulados[i] = angulos_acumulados[i-1] + angulos[i-1];
    }

    for (std::size_t i = 0; i < n; ++i) {
        std::cout << angulos[i] << ' ';
    }
    std::cout << std::endl;

    for (std::size_t i = 0; i <= n; ++i) {
        std::cout << angulos_acumulados[i] << ' ';
    }
    std::cout << std::endl;

    Ponto menor_ponto {};
    Reta menor_reta {};
    double menor_distancia { std::numeric_limits<double>::max() };
    // agora para cada ponto, encontrar o ponto/linha oposto
    double metade = (std::numbers::pi * (n - 2)) / 2;
    for (std::size_t i = 0; i < n; ++i) {
        // busca binária:
        std::size_t l = 0;
        std::size_t r = n;
        Reta encontrada {};
        while (l <= r) {
            std::size_t m = (l + r) / 2;
            std::size_t atual = i;
            std::size_t meio = (i+m >= n) ? (i+m-n) : (i+m);
            std::size_t prox = (i+m+1 >= n) ? (i+m+1-n) : (i+m+1);
            double phi_meio_atual = 0.;
            if (meio > atual) {
                phi_meio_atual = angulos_acumulados[meio+1] - angulos_acumulados[atual] - angulos[meio]/2. - angulos[atual]/2.;
            } else {
                // meio <= atual
                phi_meio_atual = angulos_acumulados[n] + angulos_acumulados[meio+1] - angulos_acumulados[atual] - angulos[meio]/2. - angulos[atual]/2.;
            }
            double phi_prox_atual = phi_meio_atual + angulos[meio]/2. + angulos[prox]/2.;
            std::cout << phi_meio_atual << ' ' << phi_prox_atual << ' ' << metade << std::endl; exit(0);
            if (phi_prox_atual >= metade) {
                if (phi_meio_atual < metade) {
                    // encontrado
                    encontrada = std::array<Ponto, 2>{fecho[meio], fecho[prox]};
                    break;
                } else {
                    r = m - 1;
                }
            } else {
                l = m + 1;
            }
        }
        double distancia = dist(fecho[i], encontrada);
        if (distancia < menor_distancia) {
            menor_distancia = distancia;
            menor_ponto = fecho[i];
            menor_reta = encontrada;
        }
    }

    return RetornoAlg {menor_ponto, menor_reta, menor_distancia};
}

void message_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, GLchar const* message, void const* user_param);

int main() {

    // std::vector<Ponto> pontos = {{15, 15}, {8, -9}, {-2, 3}, {-13, -5}, {-17, 20}};
    std::vector<Ponto> pontos = {{15, 15}, {8, -9}, {-2, 3}, {-13, -5}, {-17, 20}, {-20, -10}, {-10, 10}, {25, -10}, {-10, -15}, {10, -20}};
    for (std::size_t i = 0; i < pontos.size(); ++i) {
        std::cout << pontos[i][0] << ' ' << pontos[i][1] << std::endl;
    }
    std::cout << "Começando calculo do fecho convexo:" << std::endl;
    std::vector<Ponto> fecho = fecho_convexo(pontos);

    for (std::size_t i = 0; i < fecho.size(); ++i) {
        std::cout << fecho[i][0] << ' ' << fecho[i][1] << std::endl;
    }
    
    std::cout << "calculo da menor distancia:" << std::endl;

    RetornoAlg coiso = algoritmo(pontos);

    std::cout << "ponto: " << coiso.p[0] << ' ' << coiso.p[1] << std::endl;
    std::cout << "reta: " << coiso.r[0][0] << ' ' << coiso.r[0][1] << std::endl;
    std::cout << "    : " << coiso.r[1][0] << ' ' << coiso.r[1][1] << std::endl;
    std::cout << "distancia: " << coiso.distancia << std::endl;

    return 0;

    GLFW::Session session {};

    // tamanho padrão de tela: 800x600
    GLFW::Window::options opts { .decorated {1}, .samples {128} };
    GLFW::Window win {session, "Geometria", opts};

    GLFWwindow* window = win.justGimmeTheWindow();

    glfwMakeContextCurrent(window);

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

    /*std::array vertices {
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

    BufferLayout layout;
    layout.push<float>(3);
    layout.push<float>(2);
    Drawing cube_model {vertices, layout};
    std::array<Drawing, 10> cubes {cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,cube_model,};

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
    */

    // float mixing { 0.3 };
    // float red { 0.7 };
    while (!glfwWindowShouldClose(window)) {
        // win.processInput();
        // processar entradas??

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // program.use();
        // mixing = std::sin(static_cast<float>(glfwGetTime()) * 1.8 + 0.3);
        // program.setUniform("mixing", mixing);
        // // program.setFloat("mixing", mixing);
        // red = std::sin(static_cast<float>(glfwGetTime()) * 0.4f) * 0.7f + 0.3f;
        // program.setUniform("otherColorRed", red);

        // glm::mat4 view { glm::mat4(1.0f) };
        // // view = glm::lookAt(
        // //     win.cam.cameraPos,
        // //     win.cam.cameraPos + win.cam.cameraFront,
        // //     win.cam.cameraUp);
        // // tempor´rrio:::
        // view = glm::lookAt(glm::vec3{}, glm::vec3{}, glm::vec3{});

        // program.setViewTransform(glm::value_ptr(view));

        // glm::mat4 projection { glm::mat4(1.0f) };
        // projection = glm::perspective(
        //     glm::radians(45.f),
        //     static_cast<float>(800) / static_cast<float>(600),
        //     0.1f,
        //     100.0f);
        // program.setProjectionTransform(glm::value_ptr(projection));

        glfwSwapBuffers(window);
        session.pollEvents();
    }
    
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