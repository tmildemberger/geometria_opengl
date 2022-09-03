#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <map>
//#include <numbers>
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

#ifdef _WIN32
extern "C" {
    #define DLL_EXPORT [[gnu::dllexport]]
    DLL_EXPORT unsigned long NvOptimusEnablement = 0x00000001;
    DLL_EXPORT int AmdPowerXpressRequestHighPerformance = 1;
}
#endif

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

// double angulo_interno(Ponto p1, Ponto p2, Ponto p3) {
//     double seno = area_orientada(p1, p2, p3)/(dist(p1, p2)*dist(p1, p3));
//     double angulo = std::asin(seno);
//     if (angulo < 0.) {
//         angulo += std::numbers::pi;
//     }
//     return angulo;
// }

// ângulo entre p1p2 e p1p3
double angulo_interno(Ponto p1, Ponto p2, Ponto p3) {
    double angulo_p2 = std::atan2(p2[1] - p1[1], p2[0] - p1[0]);
    double angulo_p3 = std::atan2(p3[1] - p1[1], p3[0] - p1[0]);
    double angulo = angulo_p2 - angulo_p3;
    if (angulo < 0.) {
        angulo += 2 * 3.14159265358979323846;
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
    Ponto intersecao_encontrada;
};

std::array<double, 2> vetor_reta_ponto(Ponto p, Reta r);
double dist(Ponto p, Reta r);
RetornoAlg algoritmo(std::vector<Ponto> poligono);

std::array<double, 2> vetor_reta_ponto(Ponto p, Reta r) {
    double x3_x1 = p[0] - r[0][0];
    double x2_x1 = r[1][0] - r[0][0];
    double y3_y1 = p[1] - r[0][1];
    double y2_y1 = r[1][1] - r[0][1];
    // double produto_escalar = x2_x1*x3_x1 + y2_y1*y3_y1;
    // double tamanho_ao_quadrado = x2_x1*x2_x1 + y2_y1*y2_y1;
    // double c = (produto_escalar) / (tamanho_ao_quadrado);
    double c = (x2_x1*x3_x1 + y2_y1*y3_y1) / (x2_x1*x2_x1 + y2_y1*y2_y1);
    double dist_x = x3_x1 - x2_x1*c;
    double dist_y = y3_y1 - y2_y1*c;
    // double dist_x = p[0] - x2_x1*c - r[0][0];
    // double dist_y = p[1] - y2_y1*c - r[0][1];
    return {dist_x, dist_y};
}

double dist(Ponto p, Reta r) {
    auto [dist_x, dist_y] = vetor_reta_ponto(p, r);
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
    double metade = (3.14159265358979323846 * (n - 2)) / 2;
    for (std::size_t i = 0; i < n; ++i) {
        bool debug = false;
        if (i == 0) debug = true;
        else debug = false;
        // busca binária:
        std::size_t l = 0;
        std::size_t r = n;
        Reta encontrada {};
        Ponto a {};
        std::size_t indice_a {};
        while (l <= r) {
            std::size_t m = (l + r) / 2;
            std::size_t atual = i;
            std::size_t meio = (i+m >= n) ? (i+m-n) : (i+m);
            std::size_t prox = (i+m+1 >= n) ? (i+m+1-n) : (i+m+1);
            double phi_meio_atual = 0.;
            if (meio > atual) {
                phi_meio_atual = angulos_acumulados[meio] - angulos_acumulados[atual+1];
            } else {
                // meio <= atual
                phi_meio_atual = angulos_acumulados[n] + angulos_acumulados[meio] - angulos_acumulados[atual+1];
            }
            double phi_prox_atual = phi_meio_atual + angulos[meio] + angulos[atual] / 2.;
            if (debug) {
                std::cout << atual << ' ' << meio << ' ' << phi_meio_atual << ' ' << phi_prox_atual << ' ' << metade << std::endl;
            }
            if (phi_prox_atual >= metade) {
                if (phi_meio_atual < metade) {
                    // encontrado
                    // encontrada = std::array<Ponto, 2>{fecho[meio], fecho[prox]};
                    a = fecho[meio];
                    indice_a = meio;
                    break;
                } else {
                    r = m - 1;
                }
            } else {
                l = m + 1;
            }
        }
        // supostamente 'a' é o ponto oposto
        std::size_t indice_a_prox = (indice_a + 1 >= n) ? (indice_a + 1 - n) : (indice_a + 1);
        std::size_t indice_a_prev = (indice_a == 0) ? (n - 1) : (indice_a - 1);
        Reta r1 {fecho[indice_a], fecho[indice_a_prox]};
        Reta r2 {fecho[indice_a], fecho[indice_a_prev]};
        
        std::size_t prox = (i+1 >= n) ? (i+1-n) : (i+1);
        // std::size_t prev = (i == 0) ? (n-1) : (i-1);
        double dx = fecho[prox][0] - fecho[i][0];
        double dy = fecho[prox][1] - fecho[i][1];
        double rotacao = angulos[i] / 2.;
        double new_dx = dx * std::cos(rotacao) - dy * std::sin(rotacao);
        double new_dy = dx * std::sin(rotacao) + dy * std::cos(rotacao);
        Ponto p_bissetriz {fecho[i][0] + new_dx, fecho[i][1] + new_dy};
        if (debug) {
            std::cout << "bissetriz: " << p_bissetriz[0] << ' ' << p_bissetriz[1] << std::endl;
        }
        
        double distancia { std::numeric_limits<double>::max() };
        if (left(fecho[i], p_bissetriz, r1[0]) != left(fecho[i], p_bissetriz, r1[1])) {
            distancia = dist(fecho[i], r1);
            encontrada = r1;
            if (debug) {
                std::cout << "ala: " << std::endl;
                std::cout << fecho[i][0] << ' ' << fecho[i][1] << std::endl;
                std::cout << p_bissetriz[0] << ' ' << p_bissetriz[1] << std::endl;
                std::cout << r1[0][0] << ' ' << r1[0][1] << std::endl;
                std::cout << r1[1][0] << ' ' << r1[1][1] << std::endl;
                std::cout << left(fecho[i], p_bissetriz, r1[0]) << ' ' << left(fecho[i], p_bissetriz, r1[1]) << std::endl;
                std::cout << distancia << std::endl;
            }
        } else if (left(fecho[i], p_bissetriz, r2[0]) != left(fecho[i], p_bissetriz, r2[1])) {
            distancia = dist(fecho[i], r2);
            encontrada = r2;
            if (debug) {
                std::cout << "ala2: " << std::endl;
                std::cout << fecho[i][0] << ' ' << fecho[i][1] << std::endl;
                std::cout << p_bissetriz[0] << ' ' << p_bissetriz[1] << std::endl;
                std::cout << r1[0][0] << ' ' << r1[0][1] << std::endl;
                std::cout << r1[1][0] << ' ' << r1[1][1] << std::endl;
                std::cout << left(fecho[i], p_bissetriz, r1[0]) << ' ' << left(fecho[i], p_bissetriz, r1[1]) << std::endl;
                std::cout << distancia << std::endl;
            }
        } else {
            std::cout << "estranho: " << i << std::endl;
        }
        double dist2 {};
        {
            std::size_t indice_prox = (i + 1 >= n) ? (i + 1 - n) : (i + 1);
            std::size_t indice_prev = (i == 0) ? (n - 1) : (i - 1);
            Reta r1_a {fecho[i], fecho[indice_prox]};
            Reta r2_a {fecho[i], fecho[indice_prev]};
            
            // usar indice_a_prox
            // std::size_t prox = (i+1 >= n) ? (i+1-n) : (i+1);
            // std::size_t prev = (i == 0) ? (n-1) : (i-1);
            double dx_a = fecho[indice_a_prox][0] - fecho[indice_a][0];
            double dy_a = fecho[indice_a_prox][1] - fecho[indice_a][1];
            double rotacao_a = angulos[indice_a] / 2.;
            double new_dx_a = dx_a * std::cos(rotacao_a) - dy_a * std::sin(rotacao_a);
            double new_dy_a = dx_a * std::sin(rotacao_a) + dy_a * std::cos(rotacao_a);
            Ponto p_a_bissetriz {fecho[indice_a][0] + new_dx_a, fecho[indice_a][1] + new_dy_a};
            if (debug) {
                std::cout << "bissetriz_a: " << p_a_bissetriz[0] << ' ' << p_a_bissetriz[1] << std::endl;
            }
            if (left(fecho[indice_a], p_a_bissetriz, r1_a[0]) != left(fecho[indice_a], p_a_bissetriz, r1_a[1])) {
                dist2 = dist(fecho[indice_a], r1_a);
                if (debug) {
                    std::cout << "aaaa: " << std::endl;
                    std::cout << fecho[indice_a][0] << ' ' << fecho[indice_a][1] << std::endl;
                    std::cout << p_a_bissetriz[0] << ' ' << p_a_bissetriz[1] << std::endl;
                    std::cout << r1_a[0][0] << ' ' << r1_a[0][1] << std::endl;
                    std::cout << r1_a[1][0] << ' ' << r1_a[1][1] << std::endl;
                    std::cout << left(fecho[indice_a], p_a_bissetriz, r1_a[0]) << ' ' << left(fecho[indice_a], p_a_bissetriz, r1_a[1]) << std::endl;
                    std::cout << dist2 << std::endl;
                }
            } else if (left(fecho[indice_a], p_a_bissetriz, r2_a[0]) != left(fecho[indice_a], p_a_bissetriz, r2_a[1])) {
                dist2 = dist(fecho[indice_a], r2_a);
                if (debug) {
                    std::cout << "aaaaaaaaaaaa: " << std::endl;
                    std::cout << fecho[indice_a][0] << ' ' << fecho[indice_a][1] << std::endl;
                    std::cout << p_a_bissetriz[0] << ' ' << p_a_bissetriz[1] << std::endl;
                    std::cout << r2_a[0][0] << ' ' << r2_a[0][1] << std::endl;
                    std::cout << r2_a[1][0] << ' ' << r2_a[1][1] << std::endl;
                    std::cout << left(fecho[indice_a], p_a_bissetriz, r2_a[0]) << ' ' << left(fecho[indice_a], p_a_bissetriz, r2_a[1]) << std::endl;
                    std::cout << dist2 << std::endl;
                }
            } else {
                std::cout << "hmm: " << i << ' ' << indice_a << std::endl;
                if (debug) {
                    std::cout << "qqqqqqqqqqqqqqqq: " << std::endl;
                    std::cout << fecho[indice_a][0] << ' ' << fecho[indice_a][1] << std::endl;
                    std::cout << p_a_bissetriz[0] << ' ' << p_a_bissetriz[1] << std::endl;
                    std::cout << r1_a[0][0] << ' ' << r1_a[0][1] << std::endl;
                    std::cout << r1_a[1][0] << ' ' << r1_a[1][1] << std::endl;
                    std::cout << r2_a[0][0] << ' ' << r2_a[0][1] << std::endl;
                    std::cout << r2_a[1][0] << ' ' << r2_a[1][1] << std::endl;
                    std::cout << left(fecho[indice_a], p_a_bissetriz, r1_a[0]) << ' ' << left(fecho[indice_a], p_a_bissetriz, r1_a[1]) << std::endl;
                    std::cout << left(fecho[indice_a], p_a_bissetriz, r2_a[0]) << ' ' << left(fecho[indice_a], p_a_bissetriz, r2_a[1]) << std::endl;
                    std::cout << dist2 << std::endl;
                }
            }
        }
        if (distancia >= dist2 && distancia < menor_distancia) {
            menor_distancia = distancia;
            menor_ponto = fecho[i];
            menor_reta = encontrada;
        }
    }

    return RetornoAlg {menor_ponto, menor_reta, menor_distancia, {}};
}

double area_poligono(std::vector<Ponto> poligono);

double area_poligono(std::vector<Ponto> poligono) {
    std::size_t n = poligono.size();
    double area {};
    for (std::size_t i = 1; i < n - 1; ++i) {
        area += area_orientada(poligono[0], poligono[i], poligono[i + 1]);
    }
    // todas as áreas calculadas eram a do paralelogramo ao invés do triângulo
    // a divisão por 2 foi deixada para o final
    return area / 2.;
}

void message_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, GLchar const* message, void const* user_param);

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods);

enum class Cor {
    DESCONHECIDO,
    FORA,
    DENTRO,
};

std::tuple<double,double> intersecao(Ponto p1, Ponto p2, Ponto p3, Ponto p4);

// calcula 's' e 't' com as equações paramétricas das retas
std::tuple<double,double> intersecao(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    double x4_x3 = p4[0] - p3[0];
    double x3_x1 = p3[0] - p1[0];
    double x2_x1 = p2[0] - p1[0];
    double y4_y3 = p4[1] - p3[1];
    double y3_y1 = p3[1] - p1[1];
    double y2_y1 = p2[1] - p1[1];
    double det = x4_x3 * y2_y1 - x2_x1 * y4_y3;
    double s = (x4_x3 * y3_y1 - x3_x1 * y4_y3) / det;
    double t = (x2_x1 * y3_y1 - x3_x1 * y2_y1) / det;
    // if (std::abs(p1[0]) < 0.1 && std::abs(p1[1]) < 0.1) {
    //     std::cout << p1[0] << ' ' << p1[1] << std::endl;
    //     std::cout << p2[0] << ' ' << p2[1] << std::endl;
    //     std::cout << p3[0] << ' ' << p3[1] << std::endl;
    //     std::cout << p4[0] << ' ' << p4[1] << std::endl;
    //     std::cout << s << ' ' << t << std::endl;
    // }
    return std::make_tuple(s, t);
}

enum class Intersecao {
    PROPRIA,
    IMPROPRIA,
    NAO,
};

Intersecao intersecao_semireta_segmento(Ponto p1, Ponto p2, Ponto p3, Ponto p4);

// Intersecao intersecao_semireta_segmento(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
//     if (std::abs(p1[0]) < 0.1 && std::abs(p1[1]) < 0.1) {
//         std::cout << p1[0] << ' ' << p1[1] << std::endl;
//         std::cout << p2[0] << ' ' << p2[1] << std::endl;
//         std::cout << p3[0] << ' ' << p3[1] << std::endl;
//         std::cout << p4[0] << ' ' << p4[1] << std::endl;
//         std::cout << area_orientada(p1, p2, p3) << ' ' << area_orientada(p1, p2, p4) << std::endl;
//     }
//     if (area_orientada(p1, p2, p3) == 0. || area_orientada(p1, p2, p4) == 0.) {
//         return Intersecao::IMPROPRIA;
//     }
//     if (left(p1, p2, p3) != left(p1, p2, p4)) {
//         return Intersecao::PROPRIA;
//     }
//     return Intersecao::NAO;
// }
// checa se há interseção entre a semireta p1p2 e o segmento p3p4
Intersecao intersecao_semireta_segmento(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    auto [s, t] = intersecao(p1, p2, p3, p4);
    if (s > 0 && t > 0 && t < 1) {
        return Intersecao::PROPRIA;
    } else if (s >= 0 && t >= 0 && t <= 1) {
        return Intersecao::IMPROPRIA;
    } else {
        return Intersecao::NAO;
    }
}

Intersecao intersecao_com_left(Ponto p1, Ponto p2, Ponto p3, Ponto p4);

// checa se há interseção entre a semireta p1p2 e o segmento p3p4
Intersecao intersecao_com_left(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    double p1_p2_p3 = area_orientada(p1, p2, p3);
    double p1_p2_p4 = area_orientada(p1, p2, p4);
    double p3_p4_p1 = area_orientada(p3, p4, p1);
    double p3_p4_p2 = area_orientada(p3, p4, p2);
    if (p1_p2_p3 == 0. || p1_p2_p4 == 0. || p3_p4_p1 == 0. || p3_p4_p2 == 0.) {
        return Intersecao::IMPROPRIA;
    }
    bool left1 = p1_p2_p3 > 0.;
    bool left2 = p1_p2_p4 > 0.;
    bool left3 = p3_p4_p1 > 0.;
    bool left4 = p3_p4_p2 > 0.;
    if (left1 != left2 && left3 != left4) {
        return Intersecao::PROPRIA;
    } else {
        return Intersecao::NAO;
    }
}

Ponto ponto_intersecao(Ponto p1, Ponto p2, Ponto p3, Ponto p4);

Ponto ponto_intersecao(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    auto [s, t] = intersecao(p1, p2, p3, p4);
    return {p3[0] * (1. - t) + p4[0] * t, p3[1] * (1. - t) + p4[1] * t};
}

Cor point_in_polygon(Ponto ponto, std::vector<Ponto> poligono);

Cor point_in_polygon(Ponto ponto, std::vector<Ponto> poligono) {
    Ponto auxiliar {ponto[0] + 1.0, ponto[1]};
    std::size_t n = poligono.size();
    poligono.push_back(poligono[0]);
    std::size_t count = 0;
    for (std::size_t i = 0; i < n; ++i) {
        if (poligono[i][0] < ponto[0] && poligono[i+1][0] < ponto[0]) {
            continue;
        }
        Intersecao tipo = intersecao_semireta_segmento(
            ponto, auxiliar, poligono[i], poligono[i+1]
        );
        if (tipo == Intersecao::PROPRIA) {
            ++count;
        } else if (tipo == Intersecao::IMPROPRIA) {
            std::cout << "será" << std::endl;
            if (poligono[i][1] > ponto[1]) ++count;
            if (poligono[i+1][1] > ponto[1]) ++count;
        }
    }
    // std::cout << count << std::endl;
    if (count % 2 == 1) {
        // ímpar => dentro
        return Cor::DENTRO;
    } else {
        return Cor::FORA;
    }
}

double produto_vetorial(Ponto p1, Ponto p2, Ponto p3, Ponto p4);

double produto_vetorial(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    double v1_x = p2[0] - p1[0];
    double v1_y = p2[1] - p1[1];
    double v2_x = p4[0] - p3[0];
    double v2_y = p4[1] - p3[1];
    return v1_x * v2_y - v1_y * v2_x;
}

double produto_escalar(Ponto p1, Ponto p2, Ponto p3, Ponto p4);

double produto_escalar(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    double v1_x = p2[0] - p1[0];
    double v1_y = p2[1] - p1[1];
    double v2_x = p4[0] - p3[0];
    double v2_y = p4[1] - p3[1];
    return v1_x * v2_x + v1_y * v2_y;
}

double produto_escalar_com_ortogonal(Ponto p1, Ponto p2, Ponto p3, Ponto p4);

double produto_escalar_com_ortogonal(Ponto p1, Ponto p2, Ponto p3, Ponto p4) {
    double v1_x = p2[0] - p1[0];
    double v1_y = p2[1] - p1[1];
    double v2_x = p4[1] - p3[1]; // y e x trocados
    double v2_y = p3[0] - p4[0]; // y é x invertido
    return v1_x * v2_x + v1_y * v2_y;
}

Cor convexidade_do_vertice(std::vector<Ponto> poligono, std::size_t i);

Cor convexidade_do_vertice(std::vector<Ponto> poligono, std::size_t i) {
    auto& p = poligono;
    std::size_t prev = (i == 0) ? (p.size()-1) : (i-1);
    std::size_t prox = (i+1 >= p.size()) ? (0) : (i+1);
    auto& p1 = p[prev];
    auto& p2 = p[i];
    auto& p3 = p[prox];
    Cor nova_cor {};
    if (area_orientada(p1, p2, p3) >= 0.) {
        nova_cor = Cor::DENTRO;
    } else {
        nova_cor = Cor::FORA;
    }
    return nova_cor;
}

bool in_cone_convexo(Ponto p_i, Ponto p_j, Ponto p_i_menos, Ponto p_i_mais);

bool in_cone_convexo(Ponto p_i, Ponto p_j, Ponto p_i_menos, Ponto p_i_mais) {
    return left(p_i, p_j, p_i_menos) && left(p_j, p_i, p_i_mais);
}

bool in_cone_reflexo(Ponto p_i, Ponto p_j, Ponto p_i_menos, Ponto p_i_mais);

bool in_cone_reflexo(Ponto p_i, Ponto p_j, Ponto p_i_menos, Ponto p_i_mais) {
    return left(p_i_menos, p_i, p_j) || left(p_i, p_i_mais, p_j);
}

bool diagonal(const std::vector<Ponto>& poligono, std::size_t i, std::size_t j);

bool diagonal(const std::vector<Ponto>& poligono, std::size_t i, std::size_t j) {
    auto& p = poligono;
    std::size_t prev = (i == 0) ? (p.size()-1) : (i-1);
    std::size_t prox = (i+1 >= p.size()) ? (0) : (i+1);

    if (convexidade_do_vertice(p, i) == Cor::DENTRO) {
        // isso significa convexo (por enquanto)
        if (!in_cone_convexo(p[i], p[j], p[prev], p[prox])) {
            return false;
        }
    } else {
        if (!in_cone_reflexo(p[i], p[j], p[prev], p[prox])) {
            return false;
        }
    }

    for (std::size_t k = 0; k < p.size(); ++k) {
        auto k_prox = (k+1 >= p.size()) ? (0) : k+1;
        if (k == i || k == j || k_prox == i || k_prox == j) {
            continue;
        }
        if (intersecao_com_left(p[i], p[j], p[k], p[k_prox]) != Intersecao::NAO) {
            return false;
        }
    }
    return true;
}

bool orelha(const std::vector<Ponto>& poligono, std::size_t i);

bool orelha(const std::vector<Ponto>& poligono, std::size_t i) {
    auto& p = poligono;
    std::size_t prev = (i == 0) ? (p.size()-1) : (i-1);
    std::size_t prox = (i+1 >= p.size()) ? (0) : (i+1);
    if (convexidade_do_vertice(p, i) == Cor::DENTRO) {
        // isso significa convexo (por enquanto)
        return diagonal(poligono, prev, prox);
    } else {
        return false;
    }
}

using PoligonoComFuros = std::vector<std::vector<Ponto>>;

// PoligonoComFuros preparacao(const PoligonoComFuros& poly);

// PoligonoComFuros preparacao(const PoligonoComFuros& poly) {
//     PoligonoComFuros preparado {};
//     for (auto& comp : poly) {
//         std::vector<Ponto> nova_componente;
//         nova_componente.reserve(comp.size() * 2);
//         preparado.push_back(nova_componente);
//     }
//     return preparado;
// }

PoligonoComFuros intersecao_poligonos(PoligonoComFuros poly1, PoligonoComFuros poly2);

PoligonoComFuros intersecao_poligonos(PoligonoComFuros poly1, PoligonoComFuros poly2) {
    // considerando que cada componente já tem o primeiro e último ponto iguais
    // para poder iterar por todas as arestas dentro do loop

    // considerando que cada poligono é composto por um vetor de sequências de pontos,
    // onde a primeira é a única sequência anti-horária, e as seguintes são os buracos,
    // que devem estar inteiramente dentro do primeiro

    // std::vector<std::map<std::size_t, std::tuple<std::size_t, std::size_t, bool>>> idas(poly1.size());
    // std::vector<std::map<std::size_t, std::tuple<std::size_t, std::size_t, bool>>> voltas(poly2.size());
    // std::size_t num_intersecoes = 0;

    // // std::vector<std::vector<std::pair<Ponto, double>>> intersecoes;
    // // std::map<std::pair<std::size_t, std::size_t>, std::vector<std::tuple<Ponto, double, std::size_t, std::size_t, std::size_t>>> intersecoes;

    // for (std::size_t p1_idx = 0; p1_idx < poly1.size(); ++p1_idx) {
    //     auto& comp1 = poly1[p1_idx];
    //     for (std::size_t i = 0; i < comp1.size() - 1; ++i) {
            
    //         for (std::size_t p2_idx = 0; p2_idx < poly2.size(); ++p2_idx) {
    //             auto& comp2 = poly2[p2_idx];

    //             for (std::size_t j = 0; j < comp2.size() - 1; ++j) {
    //                 Ponto& p1 = comp1[i];
    //                 Ponto& p2 = comp1[i+1];
    //                 Ponto& p3 = comp2[j];
    //                 Ponto& p4 = comp2[j+1];
    //                 auto [s, t] = intersecao(p1, p2, p3, p4);
    //                 if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
    //                     Ponto inter {p3[0] * (1. - t) + p4[0] * t, p3[1] * (1. - t) + p4[1] * t};
    //                     bool entrando = produto_escalar_com_ortogonal(p1, p2, p3, p4) < 0;
                        
    //                 }
    //             }
    //         }
    //     }
    // }
    // for (auto& [key, vec] : intersecoes) {
    //     std::sort(vec.begin(), vec.end(), [](auto a, auto b) {
    //         return std::get<1>(a) < std::get<1>(b);
    //     });
    // }
    // std::vector<std::map<std::size_t, std::tuple<std::size_t, std::size_t, bool>>> idas(poly1.size());
    // std::vector<std::map<std::size_t, std::tuple<std::size_t, std::size_t, bool>>> voltas(poly2.size());
    // std::size_t num_intersecoes = 0;
    // for (std::size_t p2_idx = 0; p2_idx < poly2.size(); ++p2_idx) {
    //     auto& comp2 = poly2[p2_idx];
    //     novo2.push_back({});
    //     for (std::size_t j = 0; j < comp2.size() - 1; ++j) {
    //         novo2[p2_idx].push_back(comp2[j]);
    //         for (auto [ponto, posicao, p1_idx, i, pos] : intersecoes[std::make_pair(p2_idx, j)]) {
    //             novo2[p2_idx].push_back(ponto);
    //             auto poly2_pos = novo2[p2_idx].size() - 1;
    //             bool entrando = produto_escalar_com_ortogonal(
    //                 poly1[p1_idx][i],
    //                 poly1[p1_idx][i+1],
    //                 comp2[j],
    //                 comp2[j+1]
    //             ) < 0;
    //             idas[p1_idx][pos] = {p2_idx, poly2_pos, entrando};
    //             voltas[p2_idx][poly2_pos] = {p1_idx, pos, !entrando};
    //             ++num_intersecoes;
    //         }
    //     }
    // }
    // if (num_intersecoes == 0) {
    //     // nesse caso não houve interseções
    //     // preciso adicionar muitos testes depois para ver se alguma componente
    //     // de algum está dentro outra, ou algo assim
    //     return {};
    // }
    // // começa a percorrer
    // std::size_t p_idx = 0;
    // for (; p_idx < poly1.size(); ++p_idx) {

    // }

    return {};

}

/*
PoligonoComFuros intersecao_poligonos(PoligonoComFuros poly1, PoligonoComFuros poly2) {
    // considerando que cada componente já tem o primeiro e último ponto iguais
    // para poder iterar por todas as arestas dentro do loop
    PoligonoComFuros novo1 = preparacao(poly1);
    PoligonoComFuros novo2 = preparacao(poly2);
    // std::map<std::pair<std::size_t, std::size_t>, Ponto> intersecoes;
    // std::map<std::size_t, std::pair<Ponto, double>> intersecoes;
    
    // std::vector<std::vector<std::pair<Ponto, double>>> intersecoes;
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::tuple<Ponto, double, std::size_t, std::size_t, std::size_t>>> intersecoes;

    for (std::size_t p1_idx = 0; p1_idx < poly1.size(); ++p1_idx) {
        auto& comp1 = poly1[p1_idx];
        novo1.push_back({});
        for (std::size_t i = 0; i < comp1.size() - 1; ++i) {
            novo1[p1_idx].push_back(comp1[i]);
            
            for (std::size_t p2_idx = 0; p2_idx < poly2.size(); ++p2_idx) {
                auto& comp2 = poly2[p2_idx];
                // intersecoes.push_back({});
                // intersecoes[p2_idx].resize(comp2.size());
                for (std::size_t j = 0; j < comp2.size() - 1; ++j) {
                    Ponto& p1 = comp1[i];
                    Ponto& p2 = comp1[i+1];
                    Ponto& p3 = comp2[j];
                    Ponto& p4 = comp2[j+1];
                    auto [s, t] = intersecao(p1, p2, p3, p4);
                    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
                        Ponto inter {p3[0] * (1. - t) + p4[0] * t, p3[1] * (1. - t) + p4[1] * t};
                        novo1[p1_idx].push_back(inter);
                        // intersecoes[std::make_pair(i, j)] = inter;
                        // intersecoes[j] = std::make_pair(inter, t);
                        intersecoes[std::make_pair(p2_idx, j)].push_back(std::make_tuple(inter, t, p1_idx, i, novo1[p1_idx].size() - 1));
                    }
                }
            }
        }
    }
    for (auto& [key, vec] : intersecoes) {
        std::sort(vec.begin(), vec.end(), [](auto a, auto b) {
            return std::get<1>(a) < std::get<1>(b);
        });
    }
    std::vector<std::map<std::size_t, std::tuple<std::size_t, std::size_t, bool>>> idas(poly1.size());
    std::vector<std::map<std::size_t, std::tuple<std::size_t, std::size_t, bool>>> voltas(poly2.size());
    std::size_t num_intersecoes = 0;
    for (std::size_t p2_idx = 0; p2_idx < poly2.size(); ++p2_idx) {
        auto& comp2 = poly2[p2_idx];
        novo2.push_back({});
        for (std::size_t j = 0; j < comp2.size() - 1; ++j) {
            novo2[p2_idx].push_back(comp2[j]);
            for (auto [ponto, posicao, p1_idx, i, pos] : intersecoes[std::make_pair(p2_idx, j)]) {
                novo2[p2_idx].push_back(ponto);
                auto poly2_pos = novo2[p2_idx].size() - 1;
                bool entrando = produto_escalar_com_ortogonal(
                    poly1[p1_idx][i],
                    poly1[p1_idx][i+1],
                    comp2[j],
                    comp2[j+1]
                ) < 0;
                idas[p1_idx][pos] = {p2_idx, poly2_pos, entrando};
                voltas[p2_idx][poly2_pos] = {p1_idx, pos, !entrando};
                ++num_intersecoes;
            }
        }
    }
    if (num_intersecoes == 0) {
        // nesse caso não houve interseções
        // preciso adicionar muitos testes depois para ver se alguma componente
        // de algum está dentro outra, ou algo assim
        return {};
    }
    // começa a percorrer
    std::size_t p_idx = 0;
    for (; p_idx < poly1.size(); ++p_idx) {

    }

    return {};

}
*/

enum class Tela {
    ORIGINAL,
    OPERACOES_BOOLEANAS,
    TRIANGULACAO,
    ATIVIDADE,
};

struct State {
    std::vector<Ponto> cliques;
    std::vector<std::tuple<Ponto,Cor>> outros;
    float pointSize;
    bool should_recalculate_convex_hull;
    bool should_recalculate_area;
    bool should_recalculate_point_in_polygon;
    bool passo_a_passo_em_andamento;
    bool proximo_passo;
    bool passo_a_passo_acabou_de_acabar;
    bool mostrar_resultado_passo_a_passo;
    bool mostrando_resultado_passo_a_passo;
    bool comecar_passo_a_passo;

    Tela tela;

    // parte da atividade de hoje
    std::vector<Ponto> entrada;
    std::vector<Cor> cores_entrada;
    bool resetar_pontos;
    bool recebendo_pontos;
    bool recalcular_orientacao;
    bool recalcular_convexidade_dos_vertices;
    bool recalcular_orelhas;
};

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
    State& estado = *(static_cast<State*> (glfwGetWindowUserPointer(window)));
    if (estado.tela == Tela::ORIGINAL) {
        if (estado.passo_a_passo_em_andamento || estado.passo_a_passo_acabou_de_acabar || estado.mostrando_resultado_passo_a_passo) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE && mods == (GLFW_MOD_CONTROL | GLFW_MOD_SHIFT)) {
                if (!estado.passo_a_passo_em_andamento) {
                    if (estado.passo_a_passo_acabou_de_acabar) {
                        estado.mostrar_resultado_passo_a_passo = true;
                        estado.passo_a_passo_acabou_de_acabar = false;
                    } else if (estado.mostrando_resultado_passo_a_passo) {
                        estado.mostrar_resultado_passo_a_passo = false;
                    } else {
                        std::cout << "isso aqui nunca roda" << std::endl;
                    }
                } else {
                    estado.proximo_passo = true;
                }
            }
            return;
        }
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
            if (!mods) {
                double xpos {};
                double ypos {};
                glfwGetCursorPos(window, &xpos, &ypos);
                int width {};
                int height {};
                glfwGetWindowSize(window, &width, &height);
                double x {xpos / static_cast<double> (width) * 2. - 1.};
                double y {1. - ypos / static_cast<double> (height) * 2.};
                estado.cliques.push_back({x, y});
            } else if (mods == GLFW_MOD_CONTROL) {
                double xpos {};
                double ypos {};
                glfwGetCursorPos(window, &xpos, &ypos);
                int width {};
                int height {};
                glfwGetWindowSize(window, &width, &height);
                double x {xpos / static_cast<double> (width) * 2. - 1.};
                double y {1. - ypos / static_cast<double> (height) * 2.};
                estado.outros.push_back({{x, y}, Cor::DESCONHECIDO});
            } else if (mods == (GLFW_MOD_CONTROL | GLFW_MOD_SHIFT)) {
                estado.comecar_passo_a_passo = true;
            } else if (mods == GLFW_MOD_SHIFT) {
                std::cout << estado.pointSize << std::endl;
            }
        } else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
            if (!mods) {
                estado.should_recalculate_convex_hull = true;
            } else if (mods == GLFW_MOD_SHIFT) {
                estado.should_recalculate_area = true;
            } else if (mods == GLFW_MOD_CONTROL) {
                estado.should_recalculate_point_in_polygon = true;
            }
        }
    } else if (estado.tela == Tela::OPERACOES_BOOLEANAS) {

    } else if (estado.tela == Tela::ATIVIDADE) {
        if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE && estado.recebendo_pontos) {
            if (!mods) {
                double xpos {};
                double ypos {};
                glfwGetCursorPos(window, &xpos, &ypos);
                int width {};
                int height {};
                glfwGetWindowSize(window, &width, &height);
                double x {xpos / static_cast<double> (width) * 2. - 1.};
                double y {1. - ypos / static_cast<double> (height) * 2.};
                estado.entrada.push_back({x, y});
                estado.cores_entrada.push_back(Cor::DESCONHECIDO);
            } else if (mods == GLFW_MOD_CONTROL) {
                estado.recebendo_pontos = false;
            }
        } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_RELEASE && !estado.recebendo_pontos) {
            if (mods == GLFW_MOD_CONTROL) {
                estado.resetar_pontos = true;
                estado.recebendo_pontos = true;
            } else if (mods == GLFW_MOD_SHIFT) {
                estado.recalcular_orientacao = true;
            } else if (!mods) {
                estado.recalcular_convexidade_dos_vertices = true;
            } else if (mods == (GLFW_MOD_CONTROL | GLFW_MOD_SHIFT)) {
                estado.recalcular_orelhas = true;
            }
        }
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    State& estado = *(static_cast<State*> (glfwGetWindowUserPointer(window)));
    if (estado.pointSize < 60.0f) {
        estado.pointSize += static_cast<float>(yoffset);
    } else {
        estado.pointSize += static_cast<float>(yoffset) * 3.0f;
    }
    if (estado.pointSize > 320.0f) {
        estado.pointSize = 320.0f;
    } else if (estado.pointSize < 2.0f) {
        estado.pointSize = 2.0f;
    }
    glLineWidth(std::max(estado.pointSize / 2.0f, 1.0f));
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    State& estado = *(static_cast<State*> (glfwGetWindowUserPointer(window)));
    if (action == GLFW_RELEASE && !mods) {
        switch (key) {
            case GLFW_KEY_1:
                estado.tela = Tela::ORIGINAL;
                break;
            case GLFW_KEY_2:
                estado.tela = Tela::OPERACOES_BOOLEANAS;
                break;
            case GLFW_KEY_3:
                estado.tela = Tela::TRIANGULACAO;
                break;
            case GLFW_KEY_4:
                estado.tela = Tela::ATIVIDADE;
                break;
            default:
                break;
        }
    }
}

enum class Etapa {
    ETAPA_1,
    ETAPA_2,
};

struct PassoAPasso {
    Etapa etapa_do_passo_executado;
    Reta desenhar_essa;
    bool desenhar_a_outra;
    Reta tambem_desenhar_essa;
    Ponto colorir_esse;
    Ponto esse_tambem;
    bool acabou;
    RetornoAlg resultado_ate_agora;
};

PassoAPasso algoritmo_v1_passo_a_passo(std::vector<Ponto> fecho);

// recebe o fecho de início já
PassoAPasso algoritmo_v1_passo_a_passo(std::vector<Ponto> fecho) {
    
    // essa parte toda precisa -- ou só é no momento -- recalculada todo passo

    // n é o número de pontos no fecho convexo
    std::size_t n = fecho.size();
    std::vector<double> angulos(n, 0.);

    // adiciona manualmente primeiro e último ângulo
    angulos[0] = angulo_interno(fecho[0], fecho[n-1], fecho[1]);
    angulos[n-1] = angulo_interno(fecho[n-1], fecho[n-2], fecho[0]);

    for (std::size_t i = 1; i < n-1; ++i) {
        angulos[i] = angulo_interno(fecho[i], fecho[i-1], fecho[i+1]);
    }
    std::vector<double> angulos_acumulados(n+1, 0.);
    for (std::size_t i = 1; i <= n; ++i) {
        angulos_acumulados[i] = angulos_acumulados[i-1] + angulos[i-1];
    }

    // guarda o passo atual; é resetado quando uma execução se completa
    static std::size_t passo = 0;
    // outras informações necessárias
    static std::size_t ponto_atual = 0;
    static std::size_t l_atual = 0;
    static std::size_t r_atual = n;
    static bool encontrado = false;
    static std::size_t indice_encontrado = 0;

    double metade = (3.14159265358979323846 * (n - 2)) / 2;

    std::size_t i = ponto_atual;
    std::size_t m = (l_atual + r_atual) / 2;
    std::cout << l_atual << ' ' << r_atual << ' ' << i << ' ' << m << std::endl;
    // corajosamente vou supor que r_atual nunca será menor que l_atual
    std::size_t atual = i;
    std::size_t meio = (i+m >= n) ? (i+m-n) : (i+m);
    std::size_t prox = (i+m+1 >= n) ? (i+m+1-n) : (i+m+1);
    if (!encontrado && l_atual <= r_atual) {
        double phi_meio_atual = 0.;
        if (meio > atual) {
            phi_meio_atual = angulos_acumulados[meio] - angulos_acumulados[atual+1];
        } else {
            // meio <= atual
            phi_meio_atual = angulos_acumulados[n] + angulos_acumulados[meio] - angulos_acumulados[atual+1];
        }
        double phi_prox_atual = phi_meio_atual + angulos[meio] + angulos[atual] / 2.;
        if (phi_prox_atual >= metade) {
            if (phi_meio_atual < metade) {
                // encontrado
                indice_encontrado = meio;
                encontrado = true;
            } else {
                r_atual = m - 1;
            }
        } else {
            l_atual = m + 1;
        }
    }
    if (!encontrado) {
        // retorna, pois esse passo acabou sem encontrar o meio
        
        //atualiza coisas
        ++passo;

        return {Etapa::ETAPA_1, {fecho[i], fecho[meio]}, false, {}, fecho[i], fecho[meio], false, {}};
    }

    // começo da etapa 2

    // segmentos
    std::size_t indice_encontrado_prox = (indice_encontrado + 1 >= n) ? (indice_encontrado + 1 - n) : (indice_encontrado + 1);
    std::size_t indice_encontrado_prev = (indice_encontrado == 0) ? (n - 1) : (indice_encontrado - 1);
    std::size_t indice_prox = (i + 1 >= n) ? (i + 1 - n) : (i + 1);
    std::size_t indice_prev = (i == 0) ? (n - 1) : (i - 1);
    Reta op1 {fecho[indice_encontrado], fecho[indice_encontrado_prox]};
    Reta op2 {fecho[indice_encontrado], fecho[indice_encontrado_prev]};
    Reta sg1 {fecho[i], fecho[indice_prox]};
    Reta sg2 {fecho[i], fecho[indice_prev]};
    
    // primeira bissetriz
    double dx = fecho[indice_prox][0] - fecho[i][0];
    double dy = fecho[indice_prox][1] - fecho[i][1];
    double rotacao = angulos[i] / 2.;
    double new_dx = dx * std::cos(rotacao) - dy * std::sin(rotacao);
    double new_dy = dx * std::sin(rotacao) + dy * std::cos(rotacao);
    Ponto p_bissetriz {fecho[i][0] + new_dx, fecho[i][1] + new_dy};

    // segunda bissetriz
    double dx_a = fecho[indice_encontrado_prox][0] - fecho[indice_encontrado][0];
    double dy_a = fecho[indice_encontrado_prox][1] - fecho[indice_encontrado][1];
    double rotacao_a = angulos[indice_encontrado] / 2.;
    double new_dx_a = dx_a * std::cos(rotacao_a) - dy_a * std::sin(rotacao_a);
    double new_dy_a = dx_a * std::sin(rotacao_a) + dy_a * std::cos(rotacao_a);
    Ponto p_a_bissetriz {fecho[indice_encontrado][0] + new_dx_a, fecho[indice_encontrado][1] + new_dy_a};

    Ponto p = fecho[i];
    Ponto p_oposto = fecho[indice_encontrado];
    // Ponto prox_oposto = fecho[indice_encontrado_prox];
    // Ponto prev_oposto = fecho[indice_encontrado_prev];
    double distancia {};
    Reta encontrada {};
    Ponto intersecao_encontrada {};
    if (intersecao_semireta_segmento(p, p_bissetriz, op1[0], op1[1]) != Intersecao::NAO) {
        distancia = dist(p, op1);
        encontrada = op1;
        intersecao_encontrada = ponto_intersecao(p, p_bissetriz, op1[0], op1[1]);
    } else if (intersecao_semireta_segmento(p, p_bissetriz, op2[0], op2[1]) != Intersecao::NAO) {
        distancia = dist(p, op2);
        encontrada = op2;
        intersecao_encontrada = ponto_intersecao(p, p_bissetriz, op2[0], op2[1]);
    } else {
        std::cout << "estranho" << std::endl;
        distancia = std::numeric_limits<double>::max();
        intersecao_encontrada = p_bissetriz;
    }

    double outra_distancia {};
    Ponto outra_intersecao_encontrada {};
    if (intersecao_semireta_segmento(p_oposto, p_a_bissetriz, sg1[0], sg1[1]) != Intersecao::NAO) {
        outra_distancia = dist(p_oposto, sg1);
        outra_intersecao_encontrada = ponto_intersecao(p_oposto, p_a_bissetriz, sg1[0], sg1[1]);
    } else if (intersecao_semireta_segmento(p_oposto, p_a_bissetriz, sg2[0], sg2[1]) != Intersecao::NAO) {
        outra_distancia = dist(p_oposto, sg2);
        outra_intersecao_encontrada = ponto_intersecao(p_oposto, p_a_bissetriz, sg2[0], sg2[1]);
    } else {
        std::cout << "estranho2" << std::endl;
        outra_distancia = distancia + 1.;
        outra_intersecao_encontrada = p_a_bissetriz;
        // std::cout << p_oposto << std::endl;
        // dx_a
        // dy_a
        // rotacao_a
        // new_dx_a
        // new_dy_a
        // p_a_bissetriz
    }

    static double menor_distancia { std::numeric_limits<double>::max() };
    static Ponto menor_ponto {};
    static Reta menor_segmento {};
    static Ponto menor_intersecao {};
    if (outra_distancia > distancia) {
        // a coisa aconteceu (essa distância não é a procurada)
    } else if (distancia < menor_distancia) {
        menor_distancia = distancia;
        menor_ponto = p;
        menor_segmento = encontrada;
        menor_intersecao = intersecao_encontrada;
    }
    // de qualquer forma, aqui acaba outro passo, então temos que retornar
    // e fazer as atualizações necessárias
    ++ponto_atual;
    if (ponto_atual >= n) {
        // esse foi o último ponto
        // precisamos resetar tudo e retornar que acabou

        // antes criar o objeto de retorno
        PassoAPasso retorno {
            Etapa::ETAPA_2,
            {p, intersecao_encontrada},
            true,
            {p_oposto, outra_intersecao_encontrada},
            p,
            p_oposto,
            true,
            {menor_ponto, menor_segmento, menor_distancia, menor_intersecao}
        };
        passo = 0;
        ponto_atual = 0;
        l_atual = 0;
        r_atual = n;
        encontrado = false;
        indice_encontrado = 0;
        
        menor_distancia = std::numeric_limits<double>::max();
        menor_ponto = Ponto{};
        menor_segmento = Reta{};
        menor_intersecao = Ponto{};

        return retorno;
    } else {
        // o algoritmo continuará
        ++passo;
        // algumas partes são resetadas
        l_atual = 0;
        r_atual = n;
        encontrado = false;
        indice_encontrado = 0;
        return {
            Etapa::ETAPA_2,
            {p, intersecao_encontrada},
            true,
            {p_oposto, outra_intersecao_encontrada},
            p,
            p_oposto,
            false,
            {menor_ponto, menor_segmento, menor_distancia, menor_intersecao}
        };
    }
}

PassoAPasso algoritmo_guedes_v1_passo_a_passo(std::vector<Ponto> fecho);

// recebe o fecho de início já
PassoAPasso algoritmo_guedes_v1_passo_a_passo(std::vector<Ponto> fecho) {

    // n é o número de pontos no fecho convexo
    std::size_t n = fecho.size();

    static bool resetar_a = false;
    static bool resetar_b = false;


    // guarda o passo atual; é resetado quando uma execução se completa
    static std::size_t passo = 0;
    // outras informações necessárias
    static std::size_t segmento_atual = 0;
    static std::size_t l_atual = 2;
    static std::size_t r_atual = n - 1;
    static bool encontrado = false;
    static std::size_t indice_encontrado = 0;

    // usados na segunda parte
    static double menor_distancia { std::numeric_limits<double>::max() };
    static Ponto menor_ponto {};
    static Reta menor_segmento {};
    static Ponto menor_intersecao {};

    double metade = (3.14159265358979323846 * (n - 2)) / 2;

    if (resetar_a) {
        l_atual = 2;
        r_atual = n - 1;
        encontrado = false;
        indice_encontrado = 0;
        resetar_a = false;
    }
    if (resetar_b) {
        passo = 0;
        segmento_atual = 0;
        l_atual = 2;
        r_atual = n - 1;
        encontrado = false;
        indice_encontrado = 0;
        
        menor_distancia = std::numeric_limits<double>::max();
        menor_ponto = Ponto{};
        menor_segmento = Reta{};
        menor_intersecao = Ponto{};
        resetar_b = false;
    }

    std::size_t i = segmento_atual;
    std::size_t outro = (i + 1 >= n) ? (0) : (i + 1);
    std::size_t m = (l_atual + r_atual) / 2;
    // std::cout << l_atual << ' ' << r_atual << ' ' << i << ' ' << m << std::endl;
    // corajosamente vou supor que r_atual nunca será menor que l_atual
    std::size_t meio = (i+m >= n) ? (i+m-n) : (i+m);
    std::size_t prev = (meio == 0) ? (n-1) : (meio-1);
    std::size_t prox = (meio+1 >= n) ? (meio+1-n) : (meio+1);
    if (l_atual > r_atual) {
        encontrado = true;
        std::cout << "aviso" << std::endl;
    }
    // checagem redundante
    if (!encontrado && l_atual <= r_atual) {
        double primeiro = produto_vetorial(fecho[i], fecho[outro], fecho[prev], fecho[meio]);
        double segundo = produto_vetorial(fecho[i], fecho[outro], fecho[meio], fecho[prox]);
        if (primeiro >= 0. && segundo <= 0.) {
            indice_encontrado = meio;
            encontrado = true;
        } else if (primeiro < 0.) {
            r_atual = m - 1;
        } else {
            l_atual = m + 1;
        }
    }
    if (!encontrado) {
        // retorna, pois esse passo acabou sem encontrar o meio
        
        //atualiza coisas
        ++passo;

        return {Etapa::ETAPA_1, {fecho[i], fecho[outro]}, false, {}, fecho[i], fecho[meio], false, {}};
    }

    // começo da etapa 2

    Ponto p = fecho[i];
    Ponto p_o = fecho[outro];
    Ponto p_oposto = fecho[indice_encontrado];
    // Ponto prox_oposto = fecho[indice_encontrado_prox];
    // Ponto prev_oposto = fecho[indice_encontrado_prev];
    // double x3_x1 = p_oposto[0] - p[0];
    // double x2_x1 = p_o[0] - p[0];
    // double y3_y1 = p_oposto[1] - p[1];
    // double y2_y1 = p_o[1] - p[1];
    // double d_p = x2_x1*x3_x1 + y2_y1*y3_y1;
    // double tam = x2_x1*x2_x1 + y2_y1*y2_y1;
    // double c = (d_p) / (tam);
    // double dist_x = x3_x1 - x2_x1*c;
    // double dist_y = y3_y1 - y2_y1*c;
    // double i_x = x2_x1*c + p[0];
    // double i_y = y2_y1*c + p[1];
    // if (i == 2) {
    //     std::cout << p[0] << ' ' << p[1] << std::endl;
    //     std::cout << p_o[0] << ' ' << p_o[1] << std::endl;
    //     std::cout << p_oposto[0] << ' ' << p_oposto[1] << std::endl;
    //     std::cout << d_p << std::endl;
    //     std::cout << c << std::endl;
    //     std::cout << tam << std::endl;
    //     std::cout << x2_x1 << ' ' << y2_y1 << std::endl;
    //     std::cout << x2_x1*c << ' ' << y2_y1*c << std::endl;
    //     std::cout << i_x << ' ' << i_y << std::endl;
    // }
    auto diff = vetor_reta_ponto(p_oposto, {p, p_o});
    // double distancia {dist(p_oposto, {p, p_o})};
    Ponto intersecao_encontrada {p_oposto[0] - diff[0], p_oposto[1] - diff[1]};
    // Ponto intersecao_encontrada {i_x, i_y};
    double distancia {dist(p_oposto, intersecao_encontrada)};
    Reta encontrada {p, p_o};

    if (distancia < menor_distancia) {
        menor_distancia = distancia;
        menor_ponto = p_oposto;
        menor_segmento = encontrada;
        menor_intersecao = intersecao_encontrada;
    }
    // de qualquer forma, aqui acaba outro passo, então temos que retornar
    // e fazer as atualizações necessárias
    ++segmento_atual;
    if (segmento_atual >= n) {
        // esse foi o último segmento
        // precisamos resetar tudo e retornar que acabou

        // antes criar o objeto de retorno
        PassoAPasso retorno {
            Etapa::ETAPA_2,
            {p, p_o},
            true,
            {p_oposto, intersecao_encontrada},
            p,
            p_oposto,
            true,
            {menor_ponto, menor_segmento, menor_distancia, menor_intersecao}
        };

        // isso é o resetar_b:
        // passo = 0;
        // segmento_atual = 0;
        // l_atual = 2;
        // r_atual = n - 1;
        // encontrado = false;
        // indice_encontrado = 0;
        
        // menor_distancia = std::numeric_limits<double>::max();
        // menor_ponto = Ponto{};
        // menor_segmento = Reta{};
        // menor_intersecao = Ponto{};
        resetar_b = true;

        return retorno;
    } else {
        // o algoritmo continuará
        ++passo;
        // algumas partes são resetadas
        // isso é o resetar_a:
        // l_atual = 2;
        // r_atual = n - 1;
        // encontrado = false;
        // indice_encontrado = 0;
        resetar_a = true;
        return {
            Etapa::ETAPA_2,
            {p, p_o},
            true,
            {p_oposto, intersecao_encontrada},
            p,
            p_oposto,
            false,
            {menor_ponto, menor_segmento, menor_distancia, menor_intersecao}
        };
    }
}

class AlgoritmoPassoAPasso {
public:
    AlgoritmoPassoAPasso(State& arg_estado, const Shader& arg_point_program, const Shader& arg_line_program) :
        estado {arg_estado},
        point_program {arg_point_program},
        line_program {arg_line_program},
        resultado_ate_agora {},
        resultado_arrumado_para_renderizacao {false},
        quantos {0} {
        
        glGenBuffers(1, &passo_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
        glBufferData(GL_ARRAY_BUFFER, 2*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
        
        glGenVertexArrays(1, &passo_vao);
        glBindVertexArray(passo_vao);
        
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
        glEnableVertexAttribArray(1);

        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    void reset() {
        resultado_ate_agora = {};
        quantos = 0;
    }
    void proximo_passo(const std::vector<Ponto>& fecho) {
        ultimo_retorno = algoritmo_guedes_v1_passo_a_passo(fecho);
        if (ultimo_retorno.etapa_do_passo_executado == Etapa::ETAPA_2) {
            resultado_ate_agora = ultimo_retorno.resultado_ate_agora;
        }
        if (ultimo_retorno.acabou) {
            estado.passo_a_passo_em_andamento = false;
            estado.passo_a_passo_acabou_de_acabar = true;
            resultado_arrumado_para_renderizacao = false;
        }
    }
    void arruma_renderizacao() {
        std::array<Ponto, 6> pontos {};
        std::size_t num = 4;
        pontos[0] = ultimo_retorno.colorir_esse;
        pontos[1] = ultimo_retorno.esse_tambem;
        pontos[2] = ultimo_retorno.desenhar_essa[0];
        pontos[3] = ultimo_retorno.desenhar_essa[1];
        if (ultimo_retorno.desenhar_a_outra) {
            num += 2;
            pontos[4] = ultimo_retorno.tambem_desenhar_essa[0];
            pontos[5] = ultimo_retorno.tambem_desenhar_essa[1];
        }
        quantos = num;

        std::array<std::array<float, 3>, 6> cores {};
        if (ultimo_retorno.etapa_do_passo_executado == Etapa::ETAPA_1) {
            // usado como base: #26a6c9
            cores[0] = {27, 181, 224};
            cores[1] = {101, 197, 224};
            cores[2] = {38, 166, 201};
            cores[3] = {38, 166, 201};
            if (ultimo_retorno.desenhar_a_outra) {
                cores[4] = {111, 182, 201};
                cores[5] = {111, 182, 201};
            }
        } else {
            // usado como base: #c9262b
            cores[0] = {230, 32, 39};
            cores[1] = {230, 78, 83};
            cores[2] = {201, 38, 43};
            cores[3] = {201, 38, 43};
            if (ultimo_retorno.desenhar_a_outra) {
                cores[4] = {201, 71, 75};
                cores[5] = {201, 71, 75};
            }
        }
        std::vector<float> ps {};
        ps.reserve(num * 5 * sizeof (float));
        for (std::size_t i = 0; i < num; ++i) {
            ps.push_back(pontos[i][0]);
            ps.push_back(pontos[i][1]);
            ps.push_back(cores[i][0] / 255.0f);
            ps.push_back(cores[i][1] / 255.0f);
            ps.push_back(cores[i][2] / 255.0f);
        }
        glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
        glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(num * 5 * sizeof (float)), ps.data());

    }
    void renderiza_passo() {
        // isso vem junto, acho
        // glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
        glBindVertexArray(passo_vao);
        line_program.use();
        glDrawArrays(GL_LINES, 2, quantos - 2);
        point_program.use();
        glDrawArrays(GL_POINTS, 0, 2);
    }
    void renderiza_resultado() {
        if (!resultado_arrumado_para_renderizacao) {
            std::array<Ponto, 4> pontos {};
            std::size_t num = 4;
            pontos[0] = resultado_ate_agora.p;
            pontos[1] = resultado_ate_agora.intersecao_encontrada;
            pontos[2] = resultado_ate_agora.r[0];
            pontos[3] = resultado_ate_agora.r[1];
            std::vector<float> ps {};
            ps.reserve(num * 5 * sizeof (float));
            for (std::size_t i = 0; i < num; ++i) {
                ps.push_back(pontos[i][0]);
                ps.push_back(pontos[i][1]);
                ps.push_back(0.149f); // 38
                ps.push_back(0.788f); // 201
                ps.push_back(0.682f); // 174
            }
            glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
            glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(num * 5 * sizeof (float)), ps.data());

            resultado_arrumado_para_renderizacao = true;
        }

        glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
        glBindVertexArray(passo_vao);
        line_program.use();
        glDrawArrays(GL_LINES, 0, 4);
        point_program.use();
        glDrawArrays(GL_POINTS, 0, 4);
    }
private:
    State& estado;
    const Shader& point_program;
    const Shader& line_program;
    RetornoAlg resultado_ate_agora;
    PassoAPasso ultimo_retorno;
    bool resultado_arrumado_para_renderizacao = false;
    std::size_t quantos;
    unsigned int passo_vbo;
    unsigned int passo_vao;
};

int main() {

    // {
    //     // std::vector<Ponto> pontos = {{15, 15}, {8, -9}, {-2, 3}, {-13, -5}, {-17, 20}};
    //     std::vector<Ponto> pontos = {{15, 15}, {8, -9}, {-2, 3}, {-13, -5}, {-17, 20}, {-20, -10}, {-10, 10}, {25, -10}, {-10, -15}, {10, -20}};
    //     // std::vector<Ponto> pontos = {{2, 3}, {2, -3}, {-2, 3}, {-2, -3}};
    //     for (std::size_t i = 0; i < pontos.size(); ++i) {
    //         std::cout << pontos[i][0] << ' ' << pontos[i][1] << std::endl;
    //     }
    //     std::cout << "Começando calculo do fecho convexo:" << std::endl;
    //     std::vector<Ponto> fecho = fecho_convexo(pontos);

    //     for (std::size_t i = 0; i < fecho.size(); ++i) {
    //         std::cout << fecho[i][0] << ' ' << fecho[i][1] << std::endl;
    //     }
        
    //     std::cout << "calculo da menor distancia:" << std::endl;

    //     RetornoAlg coiso = algoritmo(pontos);

    //     std::cout << "ponto: " << coiso.p[0] << ' ' << coiso.p[1] << std::endl;
    //     std::cout << "reta: " << coiso.r[0][0] << ' ' << coiso.r[0][1] << std::endl;
    //     std::cout << "    : " << coiso.r[1][0] << ' ' << coiso.r[1][1] << std::endl;
    //     std::cout << "distancia: " << coiso.distancia << std::endl;
    // }

    //return 0;

    GLFW::Session session {};

    // tamanho padrão de tela: 800x600
    GLFW::Window::options opts;
    opts.window_width = 600;
    opts.version_major = 3;
    opts.version_minor = 3;
    opts.decorated = 1;
    opts.samples = 128;
    GLFW::Window win {session, "Geometria", opts};

    GLFWwindow* window = win.justGimmeTheWindow();

    glfwMakeContextCurrent(window);
    
    State estado {};
    estado.recebendo_pontos = true;
    glfwSetWindowUserPointer(window, &estado);
    
    // glfwSetCursorEnterCallback(window, cursor_enter_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    // glfwSetCursorPosCallback(window, cursor_position_callback);
    // glfwSetScrollCallback(window, scroll_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetKeyCallback(window, key_callback);

    Renderer renderer;

    /*
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    glDebugMessageCallback(message_callback, nullptr);
    */

    //glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    Shader program {"shaders/vertex.glsl", "shaders/fragment.glsl"};
    Shader point_program {"shaders/point_vertex.glsl", "shaders/point_fragment.glsl"};
    Shader color_line_program {"shaders/color_line_vertex.glsl", "shaders/color_line_fragment.glsl"};
    
    color_line_program.setFloat("alpha", 1.0f);
    /*
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
    
    unsigned vbo {};
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, 5*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    unsigned vao {};
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    
    /////////////////////////////
    unsigned outros_vbo {};
    glGenBuffers(1, &outros_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, outros_vbo);
    glBufferData(GL_ARRAY_BUFFER, 5*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    unsigned outros_vao {};
    glGenVertexArrays(1, &outros_vao);
    glBindVertexArray(outros_vao);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    
    /////////////////////////////
    unsigned fecho_vbo {};
    glGenBuffers(1, &fecho_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, fecho_vbo);
    glBufferData(GL_ARRAY_BUFFER, 2*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    unsigned fecho_vao {};
    glGenVertexArrays(1, &fecho_vao);
    glBindVertexArray(fecho_vao);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
    
    /////////////////////////////
    unsigned atividade_vbo {};
    glGenBuffers(1, &atividade_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
    glBufferData(GL_ARRAY_BUFFER, 2*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    unsigned atividade_vao {};
    glGenVertexArrays(1, &atividade_vao);
    glBindVertexArray(atividade_vao);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    
    std::vector<Ponto> fecho_calculado {};
    std::size_t last_size = 0;
    std::size_t outros_control = 0;
    estado.pointSize = 20.0f;
    glLineWidth(estado.pointSize / 2.0f);
    // estado.cliques.push_back({-.53726, -.48185});
    // estado.cliques.push_back({.37386, .09127});
    // estado.cliques.push_back({.46278, .52544});

    // RetornoAlg resultado_ate_agora {};
    // bool resultado_arrumado_para_renderizacao = false;
    // std::size_t passo_quantos {};
    std::size_t atividade_size = 0;
    AlgoritmoPassoAPasso passo_a_passo_manager {estado, point_program, color_line_program};
    while (!glfwWindowShouldClose(window)) {
        // win.processInput();
        // processar entradas??

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (estado.tela == Tela::ORIGINAL) {
            
            if (estado.cliques.size() > last_size) {
                std::size_t diff = estado.cliques.size() - last_size;
                std::vector<float> ps {};
                ps.reserve(diff * 5 * sizeof (float));
                for (std::size_t i = last_size; i < estado.cliques.size(); ++i) {
                    ps.push_back(estado.cliques[i][0]);
                    ps.push_back(estado.cliques[i][1]);
                    ps.push_back(1.0f);
                    ps.push_back(1.0f);
                    ps.push_back(1.0f);
                }
                glBindBuffer(GL_ARRAY_BUFFER, vbo);
                glBufferSubData(GL_ARRAY_BUFFER, static_cast<GLintptr>(last_size * 5 * sizeof (float)), static_cast<GLintptr>(diff * 5 * sizeof (float)), ps.data());
                last_size = estado.cliques.size();
            }

            if (estado.outros.size() > outros_control) {
                std::size_t diff = estado.outros.size() - outros_control;
                std::vector<float> ps {};
                ps.reserve(diff * 5 * sizeof (float));
                for (std::size_t i = outros_control; i < estado.outros.size(); ++i) {
                    auto [ponto, cor] = estado.outros[i];
                    ps.push_back(ponto[0]);
                    ps.push_back(ponto[1]);
                    // sempre começa com amarelo
                    ps.push_back(0.788f);
                    ps.push_back(0.682f);
                    ps.push_back(0.078f);
                }
                glBindBuffer(GL_ARRAY_BUFFER, outros_vbo);
                glBufferSubData(GL_ARRAY_BUFFER, static_cast<GLintptr>(outros_control * 5 * sizeof (float)), static_cast<GLintptr>(diff * 5 * sizeof (float)), ps.data());
                outros_control = estado.outros.size();
            }
            
            if (estado.should_recalculate_convex_hull || estado.comecar_passo_a_passo) {
                fecho_calculado = fecho_convexo(estado.cliques);
                std::vector<float> ps {};
                ps.reserve(fecho_calculado.size() * 2 * sizeof (float));
                for (std::size_t i = 0; i < fecho_calculado.size(); ++i) {
                    ps.push_back(fecho_calculado[i][0]);
                    ps.push_back(fecho_calculado[i][1]);
                }
                glBindBuffer(GL_ARRAY_BUFFER, fecho_vbo);
                glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(fecho_calculado.size() * 2 * sizeof (float)), ps.data());
                estado.should_recalculate_convex_hull = false;
            }

            if (estado.should_recalculate_area) {
                double area = area_poligono(fecho_calculado);
                std::cout << "Area atual do fecho convexo: " << area << std::endl;
                estado.should_recalculate_area = false;
            }

            if (estado.should_recalculate_point_in_polygon) {
                for (auto& [ponto, cor] : estado.outros) {
                    if (cor != Cor::DENTRO) {
                        Cor nova_cor = point_in_polygon(ponto, fecho_calculado);
                        cor = nova_cor;
                    }
                }
                // depois disso, por enquanto, recolocamos todos os dados no buffer
                // isso pode ser demorado (?)

                std::vector<float> ps {};
                ps.reserve(outros_control * 5 * sizeof (float));
                for (std::size_t i = 0; i < outros_control; ++i) {
                    auto [ponto, cor] = estado.outros[i];
                    ps.push_back(ponto[0]);
                    ps.push_back(ponto[1]);
                    // agora não vai ter nenhum amarelo
                    if (cor == Cor::DESCONHECIDO) {
                        std::cerr << "não era para ter amarelo" << std::endl;
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.682f); // 174
                        ps.push_back(0.078f); // 20
                    } else if (cor == Cor::DENTRO) {
                        ps.push_back(0.325f); // 83
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.078f); // 20
                    } else {
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.149f); // 38
                        ps.push_back(0.078f); // 20
                    }
                }
                glBindBuffer(GL_ARRAY_BUFFER, outros_vbo);
                glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(outros_control * 5 * sizeof (float)), ps.data());
                
                estado.should_recalculate_point_in_polygon = false;
            }
            
            if (fecho_calculado.size() > 0) {
                program.use();
                program.setFloat("alpha", 0.9f);
                glBindBuffer(GL_ARRAY_BUFFER, fecho_vbo);
                glBindVertexArray(fecho_vao);
                glDrawArrays(GL_LINE_LOOP, 0, fecho_calculado.size());
                
                program.setFloat("alpha", 0.4f);
                glDrawArrays(GL_TRIANGLE_FAN, 0, fecho_calculado.size());
            }
            
            point_program.use();
            point_program.setFloat("pointRadius", estado.pointSize);
            
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBindVertexArray(vao);
            glDrawArrays(GL_POINTS, 0, last_size);

            glBindBuffer(GL_ARRAY_BUFFER, outros_vbo);
            glBindVertexArray(outros_vao);
            glDrawArrays(GL_POINTS, 0, outros_control);

            if (estado.comecar_passo_a_passo) {
                estado.passo_a_passo_em_andamento = true;
                estado.proximo_passo = true;
                estado.comecar_passo_a_passo = false;
                // resultado_ate_agora = {};
                passo_a_passo_manager.reset();
            }
            if (estado.proximo_passo) {
                estado.proximo_passo = false;
                passo_a_passo_manager.proximo_passo(fecho_calculado);
                // auto res = algoritmo_v1_passo_a_passo(fecho_calculado);
                // if (res.etapa_do_passo_executado == Etapa::ETAPA_2) {
                //     resultado_ate_agora = res.resultado_ate_agora;
                // }
                // if (res.acabou) {
                //     estado.passo_a_passo_em_andamento = false;
                //     estado.passo_a_passo_acabou_de_acabar = true;
                //     resultado_arrumado_para_renderizacao = false;
                // }
                passo_a_passo_manager.arruma_renderizacao();
                // std::array<Ponto, 6> pontos {};
                // std::size_t num = 4;
                // pontos[0] = res.colorir_esse;
                // pontos[1] = res.esse_tambem;
                // pontos[2] = res.desenhar_essa[0];
                // pontos[3] = res.desenhar_essa[1];
                // if (res.desenhar_a_outra) {
                //     num += 2;
                //     pontos[4] = res.tambem_desenhar_essa[0];
                //     pontos[5] = res.tambem_desenhar_essa[1];
                // }

                // std::array<std::array<float, 3>, 6> cores {};
                // if (res.etapa_do_passo_executado == Etapa::ETAPA_1) {
                //     // usado como base: #26a6c9
                //     cores[0] = {27, 181, 224};
                //     cores[1] = {101, 197, 224};
                //     cores[2] = {38, 166, 201};
                //     cores[3] = {38, 166, 201};
                //     if (res.desenhar_a_outra) {
                //         cores[4] = {111, 182, 201};
                //         cores[5] = {111, 182, 201};
                //     }
                // } else {
                //     // usado como base: #c9262b
                //     cores[0] = {230, 32, 39};
                //     cores[1] = {230, 78, 83};
                //     cores[2] = {201, 38, 43};
                //     cores[3] = {201, 38, 43};
                //     if (res.desenhar_a_outra) {
                //         cores[4] = {201, 71, 75};
                //         cores[5] = {201, 71, 75};
                //     }
                // }
                // std::vector<float> ps {};
                // ps.reserve(num * 5 * sizeof (float));
                // for (std::size_t i = 0; i < num; ++i) {
                //     ps.push_back(pontos[i][0]);
                //     ps.push_back(pontos[i][1]);
                //     ps.push_back(cores[i][0] / 255.0f);
                //     ps.push_back(cores[i][1] / 255.0f);
                //     ps.push_back(cores[i][2] / 255.0f);
                // }
                // glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
                // glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(num * 5 * sizeof (float)), ps.data());

                // passo_a_passo_manager.renderiza_passo();
                // glBindVertexArray(passo_vao);
                // point_program.use();
                // glDrawArrays(GL_POINTS, 0, 2);
                // color_line_program.use();
                // passo_quantos = num - 2;
                // glDrawArrays(GL_LINES, 2, passo_quantos);
            }
            if (estado.passo_a_passo_em_andamento || estado.passo_a_passo_acabou_de_acabar) {
                passo_a_passo_manager.renderiza_passo();
                // glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
                // glBindVertexArray(passo_vao);
                // point_program.use();
                // glDrawArrays(GL_POINTS, 0, 2);
                // color_line_program.use();
                // glDrawArrays(GL_LINES, 2, passo_quantos);
            }
            if (estado.mostrar_resultado_passo_a_passo) {
                estado.mostrando_resultado_passo_a_passo = true;
                passo_a_passo_manager.renderiza_resultado();
                // if (!resultado_arrumado_para_renderizacao) {
                //     std::array<Ponto, 4> pontos {};
                //     std::size_t num = 4;
                //     pontos[0] = resultado_ate_agora.p;
                //     pontos[1] = resultado_ate_agora.intersecao_encontrada;
                //     pontos[2] = resultado_ate_agora.r[0];
                //     pontos[3] = resultado_ate_agora.r[1];
                //     std::vector<float> ps {};
                //     ps.reserve(num * 5 * sizeof (float));
                //     for (std::size_t i = 0; i < num; ++i) {
                //         ps.push_back(pontos[i][0]);
                //         ps.push_back(pontos[i][1]);
                //         ps.push_back(0.149f); // 38
                //         ps.push_back(0.788f); // 201
                //         ps.push_back(0.682f); // 174
                //     }
                //     glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
                //     glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(num * 5 * sizeof (float)), ps.data());

                //     resultado_arrumado_para_renderizacao = true;
                // }

                // glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
                // glBindVertexArray(passo_vao);
                // point_program.use();
                // glDrawArrays(GL_POINTS, 0, 4);
                // color_line_program.use();
                // glDrawArrays(GL_LINES, 0, 4);
            } else {
                estado.mostrando_resultado_passo_a_passo = false;
            }
            
        } else if (estado.tela == Tela::OPERACOES_BOOLEANAS) {
            glClearColor(0.2f, 0.3f, 0.2f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            
        } else if (estado.tela == Tela::ATIVIDADE) {

            glClearColor(0.3f, 0.2f, 0.3f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            if (estado.resetar_pontos) {
                estado.entrada.clear();
                estado.cores_entrada.clear();
                atividade_size = 0;
                estado.resetar_pontos = false;
            }

            if (estado.entrada.size() > atividade_size) {
                std::size_t diff = estado.entrada.size() - atividade_size;
                std::vector<float> ps {};
                ps.reserve(diff * 5 * sizeof (float));
                for (std::size_t i = atividade_size; i < estado.entrada.size(); ++i) {
                    auto ponto = estado.entrada[i];
                    auto cor = estado.cores_entrada[i];
                    ps.push_back(ponto[0]);
                    ps.push_back(ponto[1]);
                    // sempre começa com amarelo
                    ps.push_back(0.788f);
                    ps.push_back(0.682f);
                    ps.push_back(0.078f);
                }
                glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
                glBufferSubData(GL_ARRAY_BUFFER, static_cast<GLintptr>(atividade_size * 5 * sizeof (float)), static_cast<GLintptr>(diff * 5 * sizeof (float)), ps.data());
                atividade_size = estado.entrada.size();
            }

            if (estado.recalcular_orientacao) {
                auto& v = estado.entrada;
                int curvas_a_esquerda = 0;
                for (std::size_t i = 1; i < v.size() - 1; ++i) {
                    auto& p1 = v[i-1];
                    auto& p2 = v[i];
                    auto& p3 = v[i+1];
                    if (left(p1, p2, p3)) {
                        ++curvas_a_esquerda;
                    } else if (area_orientada(p1, p2, p3) != 0.) {
                        --curvas_a_esquerda;
                    }
                }
                auto& p1 = v[0];
                auto& p2 = v[1];
                auto& pn_2 = v[v.size() - 2];
                auto& pn_1 = v[v.size() - 1];
                if (left(pn_2, pn_1, p1)) {
                    ++curvas_a_esquerda;
                } else if (area_orientada(pn_2, pn_1, p1) != 0.) {
                    --curvas_a_esquerda;
                }
                if (left(pn_1, p1, p2)) {
                    ++curvas_a_esquerda;
                } else if (area_orientada(pn_1, p1, p2) != 0.) {
                    --curvas_a_esquerda;
                }
                if (curvas_a_esquerda > 0) {
                    std::cout << "orientação anti-horária" << std::endl;
                } else {
                    std::cout << "orientação horária" << std::endl;
                }
                estado.recalcular_orientacao = false;
            }

            if (estado.recalcular_convexidade_dos_vertices) {
                auto& v = estado.entrada;
                for (std::size_t i = 0; i < v.size(); ++i) {
                    std::size_t prev = (i == 0) ? (v.size()-1) : (i-1);
                    std::size_t prox = (i+1 >= v.size()) ? (0) : (i+1);
                    auto& p1 = v[prev];
                    auto& p2 = v[i];
                    auto& p3 = v[prox];
                    Cor nova_cor {};
                    if (area_orientada(p1, p2, p3) >= 0.) {
                        nova_cor = Cor::DENTRO;
                    } else {
                        nova_cor = Cor::FORA;
                    }
                    estado.cores_entrada[i] = nova_cor;
                }

                std::vector<float> ps {};
                ps.reserve(atividade_size * 5 * sizeof (float));
                for (std::size_t i = 0; i < atividade_size; ++i) {
                    auto ponto = estado.entrada[i];
                    auto cor = estado.cores_entrada[i];
                    ps.push_back(ponto[0]);
                    ps.push_back(ponto[1]);
                    // agora não vai ter nenhum amarelo
                    if (cor == Cor::DESCONHECIDO) {
                        std::cerr << "não era para ter amarelo" << std::endl;
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.682f); // 174
                        ps.push_back(0.078f); // 20
                    } else if (cor == Cor::DENTRO) {
                        ps.push_back(0.325f); // 83
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.078f); // 20
                    } else {
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.149f); // 38
                        ps.push_back(0.078f); // 20
                    }
                }
                glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
                glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(atividade_size * 5 * sizeof (float)), ps.data());
                
                estado.recalcular_convexidade_dos_vertices = false;
            }
            
            if (estado.recalcular_orelhas) {
                auto& v = estado.entrada;
                for (std::size_t i = 0; i < v.size(); ++i) {
                    Cor nova_cor {};
                    if (orelha(v, i)) {
                        nova_cor = Cor::DENTRO;
                    } else {
                        nova_cor = Cor::FORA;
                    }
                    estado.cores_entrada[i] = nova_cor;
                }

                std::vector<float> ps {};
                ps.reserve(atividade_size * 5 * sizeof (float));
                for (std::size_t i = 0; i < atividade_size; ++i) {
                    auto ponto = estado.entrada[i];
                    auto cor = estado.cores_entrada[i];
                    ps.push_back(ponto[0]);
                    ps.push_back(ponto[1]);
                    // agora não vai ter nenhum amarelo
                    if (cor == Cor::DESCONHECIDO) {
                        std::cerr << "não era para ter amarelo" << std::endl;
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.682f); // 174
                        ps.push_back(0.078f); // 20
                    } else if (cor == Cor::DENTRO) {
                        ps.push_back(0.325f); // 83
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.078f); // 20
                    } else {
                        ps.push_back(0.788f); // 201
                        ps.push_back(0.149f); // 38
                        ps.push_back(0.078f); // 20
                    }
                }
                glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
                glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(atividade_size * 5 * sizeof (float)), ps.data());
                
                estado.recalcular_orelhas = false;
            }
            
            color_line_program.use();
            color_line_program.setFloat("alpha", 0.2f);
            glDrawArrays(GL_TRIANGLE_FAN, 0, atividade_size);
            color_line_program.setFloat("alpha", 1.0f);
            glDrawArrays(GL_LINE_LOOP, 0, atividade_size);

            point_program.use();
            point_program.setFloat("pointRadius", estado.pointSize);
            glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
            glBindVertexArray(atividade_vao);
            glDrawArrays(GL_POINTS, 0, atividade_size);

        }
        //////////////////////////////////////////

        glfwSwapBuffers(window);
        session.pollEvents();
    }
    
    return 0;
}
/*
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
*/