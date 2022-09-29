#include <iostream>
#include <array>
#include <vector>
#include <deque>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <functional>
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

// checa se há interseção entre o segmento p1p2 e o segmento p3p4
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

bool point_in_polygon(Ponto ponto, std::vector<Ponto> poligono);

// recebe um polígono sem o último ponto repetido (fecha automático)
bool point_in_polygon(Ponto ponto, std::vector<Ponto> poligono) {
    Ponto auxiliar {ponto[0] + 1.0, ponto[1]};
    std::size_t n = poligono.size();
    if (poligono[0] != poligono.back()) {
        poligono.push_back(poligono[0]);
    } else {
        --n;
    }
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
        return true;
    } else {
        return false;
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

bool orientado_antihorario(const std::vector<Ponto>& poligono);

bool orientado_antihorario(const std::vector<Ponto>& poligono) {
    auto& v = poligono;
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
    return curvas_a_esquerda > 0;
}

bool abaixo(Ponto p, Ponto q);

// retorna true se 'p' está abaixo de 'q'
bool abaixo(Ponto p, Ponto q) {
    return (p[1] < q[1]) || (p[1] == q[1] && p[0] > q[0]);
}

enum class Categoria {
    START,
    SPLIT,
    MERGE,
    END,
    REGULAR
};

Categoria categoriza_ponto(std::vector<Ponto> poligono, std::size_t indice_ponto);

// categoriza o ponto nas categorias acima, usando as funções acima
Categoria categoriza_ponto(std::vector<Ponto> poligono, std::size_t indice_ponto) {
    // isso fica mais feio mas evita usar o resto da divisão (que não é otimizado nem com -O3)
    std::size_t next = (indice_ponto + 1 >= poligono.size()) ? 0 : indice_ponto + 1;
    std::size_t prev = (indice_ponto == 0) ? poligono.size() - 1 : indice_ponto - 1;

    bool anterior_abaixo = abaixo(poligono[prev], poligono[indice_ponto]);
    bool proximo_abaixo = abaixo(poligono[next], poligono[indice_ponto]);

    double angulo = angulo_interno(poligono[indice_ponto], poligono[prev], poligono[next]);

    if (anterior_abaixo && proximo_abaixo) {
        if (angulo < 3.14159265358979323846) {
            return Categoria::START;
        } else {
            return Categoria::SPLIT;
        }
    } else if (!anterior_abaixo && !proximo_abaixo) {
        if (angulo < 3.14159265358979323846) {
            return Categoria::END;
        } else {
            return Categoria::MERGE;
        }
    }
    return Categoria::REGULAR;
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

// class Compara {
// public:
//     bool operator()(const std::tuple<std::size_t, std::size_t, double, double, Ponto>& a, const std::tuple<std::size_t, std::size_t, double, double, Ponto>& b) {
//         return std::get<2>(a) < std::get<2>(b);
//     }
// };

template<typename T>
struct Par {
    T a;
    T b;
    T operator[](std::size_t i) { if (i == 0) return a; else return b; }
};

std::vector<PoligonoComFuros> op_booleana_poligonos(PoligonoComFuros poly1, PoligonoComFuros poly2, bool calcular_intersecao = true);

std::vector<PoligonoComFuros> op_booleana_poligonos(PoligonoComFuros poly1, PoligonoComFuros poly2, bool calcular_intersecao) {
    // considerando que cada componente já tem o primeiro e último ponto iguais
    // para poder iterar por todas as arestas dentro do loop

    // considerando que cada poligono é composto por um vetor de sequências de pontos,
    // onde a primeira é a única sequência anti-horária, e as seguintes são os buracos,
    // que devem estar inteiramente dentro do primeiro
    // bool compare(const std::tuple<std::size_t, std::size_t, bool, double>& a, const std::tuple<std::size_t, std::size_t, bool, double>& b) {
    //     return std::get<3>(a) < std::get<3>(b);
    // }

    std::vector<std::multimap<std::size_t, std::tuple<std::size_t, std::size_t, double, double, Ponto>>> idas(poly1.size());
    std::vector<std::multimap<std::size_t, std::tuple<std::size_t, std::size_t, double, double, Ponto>>> voltas(poly2.size());
    // std::vector<std::size_t> num_intersecoes;
    std::size_t num_intersecoes_geral = 0;

    // // std::vector<std::vector<std::pair<Ponto, double>>> intersecoes;
    // // std::map<std::pair<std::size_t, std::size_t>, std::vector<std::tuple<Ponto, double, std::size_t, std::size_t, std::size_t>>> intersecoes;

    for (std::size_t p1_idx = 0; p1_idx < poly1.size(); ++p1_idx) {
        auto& comp1 = poly1[p1_idx];
        for (std::size_t i = 0; i < comp1.size() - 1; ++i) {
            
            for (std::size_t p2_idx = 0; p2_idx < poly2.size(); ++p2_idx) {
                auto& comp2 = poly2[p2_idx];

                for (std::size_t j = 0; j < comp2.size() - 1; ++j) {
                    Ponto& p1 = comp1[i];
                    Ponto& p2 = comp1[i+1];
                    Ponto& p3 = comp2[j];
                    Ponto& p4 = comp2[j+1];
                    // std::cout << i << ' ' << comp1.size() << std::endl;
                    // std::cout << j << ' ' << comp2.size() << std::endl;
                    auto [s, t] = intersecao(p1, p2, p3, p4);
                    // std::cout << p1[0] << ' ' << p1[1] << " - " << p2[0] << ' ' << p2[1] << std::endl;
                    // std::cout << p3[0] << ' ' << p3[1] << " - " << p4[0] << ' ' << p4[1] << std::endl;
                    // std::cout << s << ' ' << t << std::endl;
                    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
                        Ponto inter {p3[0] * (1. - t) + p4[0] * t, p3[1] * (1. - t) + p4[1] * t};
                        bool entrando = produto_escalar_com_ortogonal(p1, p2, p3, p4) < 0;
                        if (calcular_intersecao) {
                            if (entrando) {
                                voltas[p2_idx].insert({j, {p1_idx, i, t, s, inter}});
                                std::cout << "voltas " << p2_idx << ' ' << j << ' ' << p1_idx << ' ' << i << ' ' << t << ' ' << s << ' ' << inter[0] << ' ' << inter[1] << std::endl;
                            } else {
                                idas[p1_idx].insert({i, {p2_idx, j, s, t, inter}});
                                std::cout << "idas " << p1_idx << ' ' << i << ' ' << p2_idx << ' ' << j << ' ' << s << ' ' << t << ' ' << inter[0] << ' ' << inter[1] << std::endl;
                            }
                        } else {
                            if (entrando) {
                                idas[p1_idx].insert({i, {p2_idx, j, s, t, inter}});
                                std::cout << "idas " << p1_idx << ' ' << i << ' ' << p2_idx << ' ' << j << ' ' << s << ' ' << t << ' ' << inter[0] << ' ' << inter[1] << std::endl;
                            } else {
                                voltas[p2_idx].insert({j, {p1_idx, i, t, s, inter}});
                                std::cout << "voltas " << p2_idx << ' ' << j << ' ' << p1_idx << ' ' << i << ' ' << t << ' ' << s << ' ' << inter[0] << ' ' << inter[1] << std::endl;
                            }
                        }
                        ++num_intersecoes_geral;
                    }
                }
            }
        }
    }
    // std::cout << num_intersecoes_geral;
    
    Par<PoligonoComFuros&> polys {poly1, poly2};
    Par<std::vector<std::multimap<std::size_t, std::tuple<std::size_t, std::size_t, double, double, Ponto>>>&> inters {idas, voltas};
    /*
    for (auto aa : inters[0][0]) {
        // bool x = decltype(aa)::nothing;
        
        auto& [nao, val] = aa;
        auto& [px_idx, ix, pos_orig, pos_novo, p_inter] = val;
        std::cout << "idas 0 " << nao << ' ' << px_idx << ' ' << ix << ' ' << ' ' << pos_orig << ' ' << pos_novo << ' ' << p_inter[0] << ' ' << p_inter[1] << std::endl;
    }
    auto aa = inters[1][0].equal_range(0);
    for (auto i = aa.first; i != aa.second; ++i) {
        // bool x = decltype(*i)::nothing;
        
        auto& [nao, val] = *i;
        auto& [px_idx, ix, pos_orig, pos_novo, p_inter] = val;
        std::cout << "voltas 0 " << nao << ' '<< px_idx << ' ' << ix << ' ' << ' ' << pos_orig << ' ' << pos_novo << ' ' << p_inter[0] << ' ' << p_inter[1] << std::endl;
    }
    return {};*/
    if (num_intersecoes_geral == 0) {
        std::size_t um = 0;
        std::size_t dois = 1;
        if (point_in_polygon(poly2[0][0], poly1[0])) {
            // nada
            ;
        } else if (point_in_polygon(poly1[0][0], poly2[0])) {
            // fazer o mesmo que acima mas trocado
            um = 1;
            dois = 0;
        } else {
            // nesse caso, nenhum está dentro de nenhum, e a interseção é nula;
            // ainda não sei o que retornar nesse caso;
            return {};
        }
        
        // poly2 (polys[dois]) estando dentro de poly1 (polys[um])
        bool dentro_de_algum_buraco = false;
        for (std::size_t i = 1; i < polys[um].size(); ++i) {
            if (point_in_polygon(polys[dois][0][0], polys[um][i])) {
                dentro_de_algum_buraco = true;
                break;
            }
        }
        if (dentro_de_algum_buraco) {
            // não há interseção
            return {};
        }
        // começa com parte externa e furos de poly2
        PoligonoComFuros retorno = polys[dois];
        std::set<std::size_t, std::greater<std::size_t>> remover;
        
        // adiciona furos de poly1 dentro de poly2 para o retorno
        // quando não estiverem dentro de um furo de poly2, e se algum
        // furo de poly2 estiver dentro de um furo de poly1, marca para
        // remover esse furo de poly2
        for (std::size_t i = 1; i < polys[um].size(); ++i) {
            if (point_in_polygon(polys[um][i][0], polys[dois][0])) {
                bool adicionar = true;
                for (std::size_t j = 1; j < polys[dois].size(); ++j) {
                    if (point_in_polygon(polys[um][i][0], polys[dois][j])) {
                        // não adicionar
                        adicionar = false;
                        break;
                    }
                    if (point_in_polygon(polys[dois][j][0], polys[um][i])) {
                        // o furo 'j' está dentro de um furo de poly1
                        remover.insert(j);
                    }
                }
                if (adicionar) {
                    retorno.push_back(polys[um][i]);
                }
            }
        }
        // remover os que foram marcados
        for (auto indice : remover) {
            // já está na ordem do maior para o menor índice
            // isso funciona mesmo se for preciso excluir todos os últimos
            std::swap(retorno[indice], retorno[retorno.size()-1]);
            retorno.pop_back();
        }
        return {retorno};
    }
    // std::cout << num_intersecoes_geral << std::endl;
    
    // para testar isso, começamos do índice 1, pois só precisamos
    // da informação para os furos
    std::set<std::size_t> sem_intersecoes_poly1;
    std::set<std::size_t> sem_intersecoes_poly2;
    for (std::size_t p1_idx = 1; p1_idx < poly1.size(); ++p1_idx) {
        if (idas[p1_idx].size() == 0) {
            sem_intersecoes_poly1.insert(p1_idx);
        }
    }
    for (std::size_t p2_idx = 1; p2_idx < poly2.size(); ++p2_idx) {
        if (voltas[p2_idx].size() == 0) {
            sem_intersecoes_poly2.insert(p2_idx);
        }
    }
    // depois disso, já não tem problema remover as interseções dos multimaps
    // à medida em que forem sendo percorridas
    
    // começa a percorrer
    std::size_t intersecoes_percorridas {0};
    
    std::vector<PoligonoComFuros> retorno;
    std::vector<std::vector<Ponto>> externos;
    std::vector<std::vector<Ponto>> furos;
    
    // auto pr = [](const std::vector<Ponto>& c) {
    //     std::cout << "{";
    //     for (auto p : c) {
    //         std::cout << " (" << p[0] << ", " << p[1] << ")";
    //     }
    //     std::cout << " }" << std::endl;
    // };
    
    // todos os caminhos começam em uma interseção
    for (std::size_t poly_sel = 0; poly_sel <= 1; ++poly_sel) {
        for (std::size_t idx_comeco = 0; idx_comeco < polys[poly_sel].size(); ++idx_comeco) {
            while (inters[poly_sel][idx_comeco].size() > 0) {
                std::size_t sel = (poly_sel + 1) % 2;
                auto& [nao, val] = *(inters[poly_sel][idx_comeco].begin());
                auto& [px_idx, ix, pos_orig, pos_novo, p_inter] = val;
                // essa é uma interseção de saída, ao contrário do planejado na aula;
                // mas funciona igual, e do jeito que eu já estou fazendo fica mais fácil
                
                // a interseção inicial de cada caminho não é removida, pra que
                // seja possível ver quando chegou ao início de novo.
                // as outras interseções encontradas no caminho são removidas
                
                Ponto inicio = p_inter;
                std::vector<Ponto> caminho;
                bool passou_por_fora = false;
                bool passou_por_buraco = false;
                // talvez isso seja desnecessário
                // caminho.reserve(poly2.size());
                
                caminho.push_back(inicio);
                // pr(caminho);
                
                // começa a procurar próximas interseções e próximos pontos
                do {
                    if (px_idx == 0) {
                        passou_por_fora = true;
                    } else {
                        passou_por_buraco = true;
                    }
                    bool acabou_segmento = true;
                    // std::cout << sel << ' ' << px_idx << ' ' << ix << std::endl;
                    auto range_outras = inters[sel][px_idx].equal_range(ix);
                    std::vector<std::pair<std::tuple<std::size_t, std::size_t, double, double, Ponto>, decltype (range_outras.first)>> ordenar;
                    for (auto outra = range_outras.first; outra != range_outras.second; ++outra) {
                        ordenar.push_back(std::make_pair((*outra).second, outra));
                    }
                    std::sort(ordenar.begin(), ordenar.end(), [](auto a, auto b) {
                        return std::get<2>(a.first) < std::get<2>(b.first);
                    });
                    std::cout << "bla     :" << ' ' << px_idx << ' ' << ix << ' ' << pos_orig << ' ' << pos_novo << ' ' << caminho.back()[0] << ' ' << caminho.back()[1] << std::endl;
                    // for (auto outra = range_outras.first; outra != range_outras.second; ++outra) {
                    //     auto& [px_idx2, ix2, pos_orig2, pos_novo2, p_inter2] = (*outra).second;
                    for (auto outra : ordenar) {
                        auto& [px_idx2, ix2, pos_orig2, pos_novo2, p_inter2] = outra.first;
                        if (px_idx == 0 && ix == 0 && px_idx2 == 1 && ix2 == 0) {
                            std::cout << "esto aqui" << ' ' << pos_orig << ' ' << pos_novo << std::endl;
                            std::cout << "     aqui" << ' ' << pos_orig2 << ' ' << pos_novo2 << std::endl;
                        }
                        if (pos_orig2 > pos_novo) {
                            // essa interseção é a próxima;
                            // remove interseção encontrada
                            inters[sel][px_idx].erase(outra.second);
                            
                            // muda pro outro polígono
                            ++sel;
                            sel %= 2;
                            px_idx = px_idx2;
                            ix = ix2;
                            pos_orig = pos_orig2;
                            pos_novo = pos_novo2;
                            
                            // coloca ponto no caminho
                            caminho.push_back(p_inter2);
                            
                            // como foi encontrada uma interseção, não chegamos ao outro
                            // ponto do segmento
                            acabou_segmento = false;
                            break;
                        }
                    }
                    if (acabou_segmento) {
                        pos_novo = 0.0;
                        caminho.push_back(polys[sel][px_idx][ix+1]);
                        ++ix;
                        if (ix >= polys[sel][px_idx].size() - 1) {
                            ix = 0;
                        }
                    }
                    // std::cout << "hmm" << std::endl;
                // pr(caminho);
                    // sai desse loop quando chegar de volta ao inicio do caminho
                } while (caminho.back() != inicio);
                
                // std::cout << "mas nao chega aqui?" << std::endl;
                
                // com o caminho completo, falta saber se é um caminho externo ou um furo
                if ((passou_por_fora && calcular_intersecao) || (!passou_por_buraco == !calcular_intersecao)) {
                    // std::cout << "aa" << std::endl;
                    externos.push_back(caminho);
                    // pr(externos.back());
                } else {
                    furos.push_back(caminho);
                }
            }
        }
    }
    
    // pr(externos.back());
    
    // se não foram encontrados caminhos exteriores, um dos dois está dentro do outro
    if (externos.size() == 0) {
        // isso significa que não existem interseções passando por nenhuma das
        // duas partes externas, e por isso sabemos que uma das duas está dentro
        // da outra (já que também sabemos que existe pelo menos uma interseção)
        
        if (point_in_polygon(poly2[0][0], poly1[0])) {
            externos.push_back(poly2[0]);
        } else {
            externos.push_back(poly1[0]);
        }
    }
    
    // para cada caminho exterior, haverá um PoligonoComFuros
    // aqui nós incluímos os furos corretos para completar
    for (std::size_t i = 0; i < externos.size(); ++i) {
        retorno.push_back({externos[i]});
        
        // adiciona todos os furos que ficam dentro dessa extremidade
        for (std::size_t j = 0; j < furos.size(); ++j) {
            if (point_in_polygon(furos[j][0], externos[i])) {
                retorno[i].push_back(furos[j]);
            }
        }
        
        // número de furos novos, usado abaixo
        std::size_t num_furos = retorno[i].size() - 1;
        
        // para um dos polígonos, só é necessário testar se os furos não
        // estão dentro de algum dos novos furos
        for (auto p1_idx : sem_intersecoes_poly1) {
            if (point_in_polygon(poly1[p1_idx][0], externos[i])) {
                bool dentro_de_algum = false;
                for (std::size_t furos_idx = 1; furos_idx < num_furos + 1; ++furos_idx) {
                    if (point_in_polygon(poly1[p1_idx][0], retorno[i][furos_idx])) {
                        dentro_de_algum = true;
                        break;
                    }
                }
                if (!dentro_de_algum) {
                    retorno[i].push_back(poly1[p1_idx]);
                }
            }
        }
        
        num_furos = retorno[i].size() - 1;
        
        // para o outro, é preciso testar além disso se algum dos furos já adicionados
        // não são internos ao furo sendo percorrido
        std::set<std::size_t, std::greater<std::size_t>> remover;
        for (auto p2_idx : sem_intersecoes_poly2) {
            if (point_in_polygon(poly2[p2_idx][0], externos[i])) {
                bool dentro_de_algum = false;
                for (std::size_t furos_idx = 1; furos_idx < num_furos + 1; ++furos_idx) {
                    if (point_in_polygon(poly2[p2_idx][0], retorno[i][furos_idx])) {
                        dentro_de_algum = true;
                        break;
                    }
                    if (point_in_polygon(retorno[i][furos_idx][0], poly2[p2_idx])) {
                        remover.insert(furos_idx);
                    }
                }
                if (!dentro_de_algum) {
                    retorno[i].push_back(poly2[p2_idx]);
                }
            }
        }
        
        // por fim, só falta remover os furos que descobrimos que estavam
        // dentro de outro furo (os índices já estão do maior ao menor)
        for (auto idx : remover) {
            // std::cout << idx << std::endl;
            std::swap(retorno[i][idx], retorno[i][retorno[i].size() - 1]);
            retorno[i].pop_back();
        }
    }
    
    return retorno;
}

/*
PoligonoComFuros op_booleana_poligonos(PoligonoComFuros poly1, PoligonoComFuros poly2) {
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

using Aresta = std::array<Ponto, 2>;

class DCEL {
public:
    // poligono_simples pode ser uma sequência anti-horária de vértices
    // ou um ciclo (com o primeiro e último vértices sendo o mesmo)
    DCEL(std::vector<Ponto> poligono_simples) {
        auto& p = poligono_simples;
        std::size_t n = p.size();

        if (p.back() == p[0]) {
            --n;
        }
        if (n < 3) {
            // o caso de n == 2 pode ser feito em outra função
            return;
        }

        vertices.reserve(n);
        edges.reserve(n * 2);
        faces.reserve(2);
        faces.push_back({nullptr});
        faces.push_back({nullptr});
        Face* outside_face = (faces.data());
        Face* inside_face = (faces.data()) + 1;
        for (std::size_t i = 0; i < n; ++i) {
            vertices.push_back({p[i], (edges.data()) + (2*i)});
            edges.push_back({(edges.data()) + (2*i) + 1, (vertices.data()) + i, nullptr, nullptr, inside_face});
            edges.push_back({(edges.data()) + (2*i), (vertices.data()) + i + 1, nullptr, nullptr, outside_face});
        }
        edges.back().origin = (vertices.data());

        edges[0].next = &edges[2];
        edges[0].prev = &edges[2*n - 2];

        edges[2*n - 2].next = &edges[0];
        edges[2*n - 2].prev = &edges[2*n - 4];
        for (std::size_t i = 1; i < n - 1; ++i) {
            edges[2*i].next = &edges[2*(i+1)];
            edges[2*i].prev = &edges[2*(i-1)];
        }

        for (std::size_t i = 0; i < n; ++i) {
            edges[2*i + 1].prev = edges[2*i].next->twin;
            edges[2*i + 1].next = edges[2*i].prev->twin;
        }

        outside_face->edge = &edges[1];
        inside_face->edge = &edges[0];
    }

    // recebe índice de uma face, e retorna um vetor com as arestas dela
    // percorridas de acordo com a DCEL
    std::vector<Aresta> arestas_de_uma_face(std::size_t face) {
        if (face >= faces.size()) {
            return {};
        }
        Edge* e = faces[face].edge;
        Vertex* start = e->origin;
        std::vector<Aresta> retorno;
        retorno.push_back(Aresta{e->origin->xy, e->twin->origin->xy});
        e = e->next;
        while (e->face == &faces[face] && e->origin != start) {
            retorno.push_back(Aresta{e->origin->xy, e->twin->origin->xy});
            e = e->next;
        }
        if (e->face != &faces[face]) {
            // aqui já sabemos que houve um erro kk
            std::cerr << "erro em 'arestas_de_uma_face'" << std::endl;
        }
        return retorno;
    }

    // recebe índice de um vértice, e retorna um vetor com as arestas que o orbitam
    // percorridas de acordo com a DCEL
    std::vector<Aresta> orbita_de_um_vertice(std::size_t vertice) {
        if (vertice >= vertices.size()) {
            return {};
        }
        Edge* e = vertices[vertice].edge;
        Edge* start = e;
        std::vector<Aresta> retorno;
        retorno.push_back(Aresta{e->origin->xy, e->twin->origin->xy});
        e = e->prev->twin;
        while (e->origin == &vertices[vertice] && e != start) {
            retorno.push_back(Aresta{e->origin->xy, e->twin->origin->xy});
            e = e->prev->twin;
        }
        if (e->origin != &vertices[vertice]) {
            // o mesmo do de cima
            std::cerr << "erro em 'orbita_de_um_vertice'" << std::endl;
        }
        return retorno;
    }

    // recebe dois índices de vértice, e adiciona à DCEL a aresta entre os vértices
    // se ainda não existe
    void inclui_aresta(std::size_t v1_i, std::size_t v2_i) {
        if (v1_i >= vertices.size() || v2_i >= vertices.size() || v1_i == v2_i) {
            return;
        }
        Vertex* v1 = &vertices[v1_i];
        Vertex* v2 = &vertices[v2_i];

        // aqui faço uma reimplementação da órbita do vértice
        // por que? não sei
        Edge* e = v1->edge;
        Edge* start = e;
        Edge* last = e;
        if (e->twin->origin == v2) {
            // aresta já existe
            return;
        }
        e = e->prev->twin;
        Edge* found = nullptr;
        while (e != start) {
            if (e->twin->origin == v2) {
                // aresta já existe
                return;
            }
            if (!left(e->twin->origin->xy, v1->xy, last->origin->xy)) {
                // isso significa que v1 "entre" 'last' e 'e' é "reflexo" kkk
                if (left(v1->xy, last->twin->origin->xy, v2->xy) || !left(v1->xy, e->twin->origin->xy, v2->xy)) {
                    // encontramos a aresta que deve compartilhar a face com o v2
                    found = last;
                    break;
                }
            } else {
                // no caso do vértice ser convexo é normal
                if (left(v1->xy, last->twin->origin->xy, v2->xy) && !left(v1->xy, e->twin->origin->xy, v2->xy)) {
                    // encontramos a aresta que deve compartilhar a face com o v2
                    found = last;
                    break;
                }
            }
            last = e;
            e = e->prev->twin;
        }
        if (!found) {
            // houve algum erro, o ponto tem que estar entre duas arestas
            std::cerr << "erro em 'inclui_aresta'" << std::endl;
            return;
        }

        // agora temos que encontrar a aresta nessa face que vai até v2
        start = found;
        e = start->next;
        found = nullptr;
        while (e != start) {
            if (e->origin == v2) {
                // a aresta anterior é a que vai para v2
                found = e->prev;
                break;
            }
            e = e->next;
        }
        if (!found) {
            // nesse caso não necessariamente houve um erro;
            // o caso em que v1 e v2 não compartilham uma face cai aqui
            return;
        }

        // reserva espaço para uma nova face e duas novas arestas
        reserva_espacos(1, 2, 0);
        faces.push_back({nullptr});

        std::size_t idx = edges.size();
        edges.push_back({(edges.data()) + idx + 1, v1, (edges.data()) + idx + 1, (edges.data()) + idx + 1, found->face});
        edges.push_back({(edges.data()) + idx, v2, (edges.data()) + idx, (edges.data()) + idx, found->face});

        conecta_arestas(found, &edges[idx]);
        conecta_arestas(start->prev, &edges[idx + 1]);

        edges[idx + 1].face = &faces.back();
        faces.back().edge = &edges[idx + 1];
        e = edges[idx + 1].next;
        while (e != &edges[idx + 1]) {
            e->face = &faces.back();
        }

        edges[idx].face->edge = &edges[idx];
        edges[idx + 1].face->edge = &edges[idx + 1];

    }

    // recebe índice de uma aresta e um valor de 0 a 1 indicando onde colocar o
    // novo vértice
    void inclui_vertice_em_aresta(std::size_t aresta, double onde) {
        if (onde <= 0 || onde >= 1 || aresta >= edges.size()) {
            return;
        }
        Edge* e = &edges[aresta];
        Vertex* v1 = e->origin;
        Vertex* v2 = e->twin->origin;
        Ponto p = {(v2->xy[0]-v1->xy[0])*onde+v1->xy[0], (v2->xy[1]-v1->xy[1])*onde+v1->xy[1]};

        // reserva espaço para um novo vértice e duas novas arestas
        reserva_espacos(0, 2, 1);

        vertices.push_back({p, e});
        Vertex* v3 = &vertices.back();
        Edge* e8 = e->prev;
        Edge* e9 = e->twin->next;

        e->origin = v3;

        std::size_t idx = edges.size();
        edges.push_back({(edges.data()) + idx + 1, v1, e, e8, e->face});
        edges.push_back({(edges.data()) + idx, v3, e9, e->twin, e->twin->face});

        Edge* e2 = &edges[idx];

        v1->edge = e2;

        e8->next = e2;
        e9->prev = e2->twin;
        e->prev = e2;
        e->twin->next = e2->twin;
    }

    // tenta encontrar a face em que 'p' está
    std::size_t qual_face(Ponto p) {
        // começa em um vértice qualquer
        Vertex* v1 = &vertices[0];

        // aqui faço novamente uma reimplementação da órbita do vértice
        // por que? nesse caso acho que é porque
        Face* found_face = nullptr;
        while (!found_face) {
            
            Edge* e = v1->edge;
            Edge* start = e;
            Edge* last = nullptr;
            e = e->prev->twin;
            Edge* found = nullptr;
            while (last != start) {
                last = e->twin->next;
                if (!left(e->twin->origin->xy, v1->xy, last->twin->origin->xy)) {
                    // isso significa que v1 "entre" 'last' e 'e' é "reflexo" kkk
                    if (left(v1->xy, last->twin->origin->xy, p) || !left(v1->xy, e->twin->origin->xy, p)) {
                        // encontramos uma aresta candidata a compartilhar a face com p
                        found = last;
                        break;
                    }
                } else {
                    // no caso do vértice ser convexo é normal
                    if (left(v1->xy, last->twin->origin->xy, p) && !left(v1->xy, e->twin->origin->xy, p)) {
                        // encontramos uma aresta candidata a compartilhar a face com p
                        found = last;
                        break;
                    }
                }
                last = e;
                e = e->prev->twin;
            }
            if (!found) {
                // houve algum erro, o ponto tem que estar entre duas arestas
                std::cerr << "erro1 em 'qual_face'" << std::endl;
                return 0;
            }

            // agora vamos calcular um segmento de 'v1' até 'p'
            // se 'p' estiver na face de 'v1', vai existir um vértice 'v' (não necessariamente 'v1')
            // tal que o segmento entre 'v' e 'p' não tem interseção com nenhuma aresta da face

            // assim, percorremos a face de 'v1' enquanto não houver interseção entre a
            // a aresta atual e o segmento.
            // toda vez que encontrarmos uma interseção com uma aresta 'a', o vértice de destino
            // de 'a' se tornará o nosso novo 'v1', começando de novo com a parte de encontrar a nova
            // face candidata e calculando o segmento.
            // só paramos quando uma face inteira tiver sido percorrida sem encontrar interseções.
            Edge* a = found->next;

            // Aresta segmento_v1_p = {v1->xy, p};
            // while (intersecao_com_left(a->origin->xy, a->twin->origin->xy, v1->xy, p) == Intersecao::NAO) {
            while (true) {
                std::cout << "testando aresta " << a - edges.data() << " e " << (intersecao_com_left(a->origin->xy, a->twin->origin->xy, v1->xy, p) == Intersecao::NAO) << std::endl;
                if (intersecao_com_left(a->origin->xy, a->twin->origin->xy, v1->xy, p) != Intersecao::NAO) {
                    break;
                }
                a = a->next;
                if (a == found || a->next == found) {
                    found_face = found->face;
                }
            }
            v1 = a->twin->origin;
            // se não tiver encontrado, no próximo loop recomeçaremos orbitando ao
            // redor do vértice v1 acima
        }
        return static_cast<std::size_t>(found_face - faces.data());
    }

private:

    struct Face;
    struct Edge;
    struct Vertex;

    struct Face {
        Edge* edge;
    };
    struct Edge {
        Edge* twin;
        Vertex* origin;
        Edge* next;
        Edge* prev;
        Face* face;
    };
    struct Vertex {
        Ponto xy;
        Edge* edge;
    };

    std::vector<Face> faces;
    std::vector<Edge> edges;
    std::vector<Vertex> vertices;

    void reserva_espacos(std::size_t n_faces, std::size_t n_edges, std::size_t n_vertices) {
        std::size_t faces_cap = faces.capacity();
        std::size_t edges_cap = edges.capacity();
        std::size_t vertices_cap = vertices.capacity();

        if (n_faces) {
            Face* base = faces.data();
            faces.reserve(faces.size() + n_faces);
            if (faces.capacity() != faces_cap) {
                // recalcula todos os ponteiros para face
                for (auto& edge : edges) {
                    edge.face = (faces.data()) + (edge.face - base);
                }
            }
        }

        if (n_vertices) {
            Vertex* base = vertices.data();
            vertices.reserve(vertices.size() + n_vertices);
            if (vertices.capacity() != vertices_cap) {
                // recalcula todos os ponteiros para vértice
                for (auto& edge : edges) {
                    edge.origin = (vertices.data()) + (edge.origin - base);
                }
            }
        }

        if (n_edges) {
            Edge* base = edges.data();
            edges.reserve(edges.size() + n_edges);
            if (edges.capacity() != edges_cap) {
                // recalcula todos os ponteiros para aresta
                for (auto& face : faces) {
                    face.edge = (edges.data()) + (face.edge - base);
                }
                for (auto& vertex : vertices) {
                    vertex.edge = (edges.data()) + (vertex.edge - base);
                }
                for (auto& edge : edges) {
                    edge.twin = (edges.data()) + (edge.twin - base);
                    edge.next = (edges.data()) + (edge.next - base);
                    edge.prev = (edges.data()) + (edge.prev - base);
                }
            }
        }
    }

    void conecta_arestas(Edge* e1, Edge* e2) {
        // conecta e1, e1->next, e2 e e2->twin corretamente, como na aula
        e2->twin->next = e1->next;
        e2->prev = e1;
        e1->next->prev = e2->twin;
        e1->next = e2;
    }
};



enum class Tela {
    ORIGINAL,
    OPERACOES_BOOLEANAS,
    TRIANGULACAO,
    ATIVIDADE,
    DCEL_TESTE,
};

enum class Dcel_Data {
    RESETANDO,
    RECEBENDO,
    CRIANDO_DCEL,
    DCEL_PRONTA,
    PISCANDO,
};

enum class Dcel_Op {
    ENCONTRAR_FACE,
    PISCAR_FACE,
};

struct Dcel_Args {
    Ponto ponto;
};

struct Dcel_Ops {
    Dcel_Args args;
    Dcel_Op op;
};

struct Dcel_Teste_State {
    Dcel_Data estado;
    std::vector<Ponto> poly;
    bool ponto_adicionado;
    bool poligono_fechado;
    std::deque<Dcel_Ops> operacoes;
};

struct State {
    Dcel_Teste_State estado_dcel_teste;

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
    std::vector<Categoria> cores_entrada;
    bool resetar_pontos;
    bool recebendo_pontos;
    bool recalcular_orientacao;
    bool recalcular_convexidade_dos_vertices;
    bool recalcular_orelhas;
    bool recalcular_visivel;
    bool visivel_pronto;
    Ponto observador;

    // parte das operações booleanas
    std::array<PoligonoComFuros, 2> polys;
    std::vector<PoligonoComFuros> intersecoes;
    std::array<bool, 2> mostrando_polys;
    bool mostrando_intersecoes;
    std::array<bool, 2> limpar_polys;
    std::array<bool, 2> limpar_ultimo_polys;
    bool limpar_intersecoes;
    bool recalcular_intersecoes;

    std::array<bool, 2> recebendo_polys;
    std::array<std::size_t, 2> polys_prontos;
};

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
    auto ponto_xy = [window]() -> Ponto {
        double xpos {};
        double ypos {};
        glfwGetCursorPos(window, &xpos, &ypos);
        int width {};
        int height {};
        glfwGetWindowSize(window, &width, &height);
        double x {xpos / static_cast<double> (width) * 2. - 1.};
        double y {1. - ypos / static_cast<double> (height) * 2.};
        return {x, y};
    };
    State& estado = *(static_cast<State*> (glfwGetWindowUserPointer(window)));
    auto coloca_ponto_dcel = [ponto_xy, &estado]() {
        auto& estad = estado.estado_dcel_teste;
        auto& p = estad.poly;
        Ponto ponto = ponto_xy();
        if (p.size() >= 3) {
            for (std::size_t i = 0; i < p.size() - 2; ++i) {
                if (intersecao_com_left(p[i], p[i+1], p[p.size()-1], ponto) != Intersecao::NAO) {
                    return;
                }
            }
        }
        p.push_back(ponto);
        estad.ponto_adicionado = true;
    };
    auto fecha_poligono_dcel = [&estado]() {
        auto& estad = estado.estado_dcel_teste;
        auto& p = estad.poly;
        if (p.size() <= 2) {
            return;
        }
        for (std::size_t i = 1; i < p.size() - 2; ++i) {
            if (intersecao_com_left(p[i], p[i+1], p[p.size()-1], p[0]) != Intersecao::NAO) {
                return;
            }
        }
        bool orientado_certo = orientado_antihorario(p);
        if (!orientado_certo) {
            estad.estado = Dcel_Data::RESETANDO;
            return;
        }
        estad.poligono_fechado = true;
        estad.estado = Dcel_Data::CRIANDO_DCEL;
    };
    auto coloca_ponto_poly = [ponto_xy, &estado](std::size_t qual) {
        Ponto ponto = ponto_xy();
        std::size_t px_idx = estado.polys_prontos[qual];
        if (px_idx != 0) {
            if (!point_in_polygon(ponto, estado.polys[qual][0])) {
                return;
            }
            for (std::size_t i = 1; i < px_idx; ++i) {
                if (point_in_polygon(ponto, estado.polys[qual][i])) {
                    return;
                }
            }
        }
        if (estado.polys[qual][px_idx].size() >= 3) {
            for (std::size_t i = 0; i < estado.polys[qual][px_idx].size() - 2; ++i) {
                auto& p = estado.polys[qual][px_idx];
                if (intersecao_com_left(p[i], p[i+1], p[p.size()-1], ponto) != Intersecao::NAO) {
                    return;
                }
            }
        }
        auto& p = estado.polys[qual][px_idx];
        if (p.size() >= 1) {
            for (std::size_t j = 0; j < estado.polys_prontos[qual]; ++j) {
                for (std::size_t i = 0; i < estado.polys[qual][j].size() - 2; ++i) {
                    auto& q = estado.polys[qual][j];
                    if (intersecao_com_left(q[i], q[i+1], p[p.size()-1], ponto) != Intersecao::NAO) {
                        return;
                    }
                }
            }
        }
        estado.polys[qual][px_idx].push_back(ponto);
        estado.recebendo_polys[qual] = true;
        estado.mostrando_polys[qual] = true;
    };
    auto fecha_poligono_poly = [&estado](std::size_t qual) {
        std::size_t px_idx = estado.polys_prontos[qual];
        auto& vec = estado.polys[qual][px_idx];
        if (vec.size() <= 2) {
            return;
        }
        auto& p = estado.polys[qual][px_idx];
        for (std::size_t j = 0; j <= px_idx; ++j) {
            auto& q = estado.polys[qual][j];
            std::size_t start = px_idx == j;
            std::size_t end = (px_idx == j) + 1;
            for (std::size_t i = start; i < q.size() - end; ++i) {
                if (intersecao_com_left(q[i], q[i+1], p[p.size()-1], p[0]) != Intersecao::NAO) {
                    return;
                }
            }
        }
        bool orientado_errado = orientado_antihorario(estado.polys[qual][px_idx]);
        if (px_idx == 0) {
            orientado_errado = !orientado_errado;
        }
        if (orientado_errado) {
            estado.limpar_ultimo_polys[qual] = true;
            return;
        }
        for (std::size_t j = 1; j < px_idx; ++j) {
            auto& q = estado.polys[qual][j];
            if (point_in_polygon(q[0], vec)) {
                estado.limpar_ultimo_polys[qual] = true;
                return;
            }
        }
        vec.push_back(vec[0]);
        ++estado.polys_prontos[qual];
        // estado.polys[qual].push_back({});
        estado.recebendo_polys[qual] = false;
    };
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
        if (action != GLFW_RELEASE) return;
        if (button == GLFW_MOUSE_BUTTON_MIDDLE && mods == (GLFW_MOD_SHIFT | GLFW_MOD_CONTROL)) {
            estado.limpar_intersecoes = true;
        }
        if (!estado.recebendo_polys[0]) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && !mods) {
                estado.mostrando_polys[0] = !estado.mostrando_polys[0];
            } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && mods == GLFW_MOD_SHIFT) {
                estado.limpar_intersecoes = true;
                estado.limpar_polys[0] = true;
            }
        } else if (button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_SHIFT) {
            coloca_ponto_poly(0);
        } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && mods == GLFW_MOD_SHIFT) {
            // fecha o polígono
            fecha_poligono_poly(0);
        } else if (button == GLFW_MOUSE_BUTTON_RIGHT && mods == GLFW_MOD_SHIFT) {
            estado.limpar_ultimo_polys[0] = true;
            estado.recebendo_polys[0] = false;
        }
        if (!estado.recebendo_polys[1]) {
            if (button == GLFW_MOUSE_BUTTON_MIDDLE && !mods) {
                estado.mostrando_polys[1] = !estado.mostrando_polys[1];
            } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && mods == GLFW_MOD_CONTROL) {
                estado.limpar_intersecoes = true;
                estado.limpar_polys[1] = true;
            }
        } else if (button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_CONTROL) {
            coloca_ponto_poly(1);
        } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && mods == GLFW_MOD_CONTROL) {
            // fecha o polígono
            fecha_poligono_poly(1);
        } else if (button == GLFW_MOUSE_BUTTON_RIGHT && mods == GLFW_MOD_CONTROL) {
            estado.limpar_ultimo_polys[1] = true;
            estado.recebendo_polys[1] = false;
        }
        if (!estado.recebendo_polys[0] && !estado.recebendo_polys[1]) {
            if (button == GLFW_MOUSE_BUTTON_RIGHT && !mods) {
                estado.mostrando_intersecoes = !estado.mostrando_intersecoes;
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_SHIFT) {
                estado.polys[0].push_back({});
                coloca_ponto_poly(0);
                if (estado.polys[0].back().size() == 0) {
                    estado.polys[0].pop_back();
                }
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_CONTROL) {
                estado.polys[1].push_back({});
                coloca_ponto_poly(1);
                if (estado.polys[1].back().size() == 0) {
                    estado.polys[1].pop_back();
                }
            } else if (button == GLFW_MOUSE_BUTTON_RIGHT && mods == (GLFW_MOD_SHIFT | GLFW_MOD_CONTROL)) {
                estado.recalcular_intersecoes = true;
            }
        }
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
                estado.cores_entrada.push_back(Categoria::REGULAR);
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
            } else if (mods == GLFW_MOD_ALT) {
                estado.observador = ponto_xy();
                estado.visivel_pronto = false;
                estado.recalcular_visivel = true;
            } else if (mods == (GLFW_MOD_ALT | GLFW_MOD_SHIFT)) {
                estado.visivel_pronto = !estado.visivel_pronto;
            }
        }
    } else if (estado.tela == Tela::DCEL_TESTE) {
        if (action != GLFW_RELEASE) return;
        auto& estad = estado.estado_dcel_teste;
        auto& p = estad.poly;
        if (p.size() == 0) {
            estad.estado = Dcel_Data::RECEBENDO;
        }
        if (estad.estado == Dcel_Data::RECEBENDO) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && !mods) {
                coloca_ponto_dcel();
            } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && !mods) {
                fecha_poligono_dcel();
            } else if (button == GLFW_MOUSE_BUTTON_RIGHT && !mods) {
                estad.estado = Dcel_Data::RESETANDO;
            }
        } else if (estad.estado == Dcel_Data::DCEL_PRONTA) {
            if (button == GLFW_MOUSE_BUTTON_RIGHT && !mods) {
                estad.estado = Dcel_Data::RESETANDO;
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && !mods) {
                Ponto clicado = ponto_xy();
                estad.operacoes.push_back({{clicado}, Dcel_Op::ENCONTRAR_FACE});
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_SHIFT) {
                Ponto clicado = ponto_xy();
                estad.operacoes.push_back({{clicado}, Dcel_Op::PISCAR_FACE});
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
            case GLFW_KEY_5:
                estado.tela = Tela::DCEL_TESTE;
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
    // std::cout << l_atual << ' ' << r_atual << ' ' << i << ' ' << m << std::endl;
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

std::vector<Ponto> regiao_visivel(std::vector<Ponto> poligono, Ponto p);

std::vector<Ponto> regiao_visivel(std::vector<Ponto> poligono, Ponto p) {
    // para facilitar nos loops, o primeiro vértice do polígono
    // é copiado para o final do vetor
    poligono.push_back(poligono[0]);

    std::vector<bool> visiveis(poligono.size());
    std::multimap<std::size_t, std::pair<Ponto, double>> inters;
    // std::vector<std::tuple<Ponto, double, std::size_t>> inters;
    for (std::size_t i = 1; i < poligono.size(); ++i) {
        bool visivel = true;
        for (std::size_t j = 0; j < poligono.size() - 1; ++j) {
            if (j == i || j == i - 1 || (j == 0 && i == poligono.size() - 1)) {
                continue;
            }
            auto [s, t] = intersecao(p, poligono[i], poligono[j], poligono[j + 1]);

            // if (intersecao_com_left(p, poligono[i], poligono[j], poligono[j + 1]) != Intersecao::NAO) {
            if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
                visivel = false;
                break;
            }
        }
        if (visivel) {
            visiveis[i] = true;
            auto& poly = poligono;
            // std::size_t prev = (i == 0) ? (poly.size()-1) : (i-1);
            std::size_t prev = i - 1;
            std::size_t prox = (i+1 >= poly.size()) ? (1) : (i+1);
            if (!left(poly[prev], poly[i], poly[prox]) && !(left(poly[prev], poly[i], p) && left(poly[i], poly[prox], p))) {
                // vai ter interseção com outra aresta
                std::size_t menor_aresta = poly.size();
                double menor_distancia = std::numeric_limits<double>::infinity();
                double t_do_menor = std::numeric_limits<double>::infinity();
                for (std::size_t j = 0; j < poligono.size() - 1; ++j) {
                    if (j == i || j == i - 1 || (j == 0 && i == poligono.size() - 1)) {
                        continue;
                    }
                    auto [s, t] = intersecao(p, poly[i], poly[j], poly[j + 1]);
                    if (t > 0 && t < 1 && s > 0 && s < menor_distancia) {
                        menor_distancia = s;
                        menor_aresta = j;
                        t_do_menor = t;
                    }
                }
                if (menor_aresta == poly.size()) {
                    // deveria ter pelo menos uma interseção
                    std::cerr << "erro em 'regiao_visivel" << std::endl;
                    return {};
                }
                double dx = menor_distancia * (poly[i][0] - p[0]);
                double dy = menor_distancia * (poly[i][1] - p[1]);
                Ponto inter = {p[0] + dx, p[1] + dy};
                inters.insert({menor_aresta, {inter, t_do_menor}});
            }
        }
    }
    std::vector<Ponto> retorno;
    for (std::size_t i = 0; i < poligono.size(); ++i) {
        if (visiveis[i]) {
            retorno.push_back(poligono[i]);
        }
        std::vector<std::pair<Ponto, double>> sub_inters;
        auto [begin, end] = inters.equal_range(i);
        while (begin != end) {
            sub_inters.push_back(begin->second);
            ++begin;
        }
        std::sort(sub_inters.begin(), sub_inters.end(), [](auto a, auto b) {
            return a.second < b.second;
        });
        for (auto& [ponto, t] : sub_inters) {
            retorno.push_back(ponto);
        }
    }
    return retorno;
}

struct Coisas_Para_Piscar {
    std::size_t ticks;
    std::size_t ticks_por_aresta;
    std::size_t atual;
    std::vector<Aresta> arestas;
};

struct CoisasDCEL {
    unsigned vao;
    unsigned vbo;
    unsigned ebo;
    std::size_t last_size;
    std::size_t edge_count;
    std::unique_ptr<DCEL> dcel_ptr;
    Coisas_Para_Piscar coisas_piscar;
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
    // SEMPRE LEMBRAR DE REPETIR O PRIMEIRO E ULTIMO VERTICE
    // int meeeeeeee = 0;
    // int meu = 0;
    /*auto res = op_booleana_poligonos(
        {
            {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}, {-1, -1}},
            {{0.1, 0.9}, {0.35, 0.93}, {0.75, 0.13}, {0.1, 0.9}},
            {{0.679, 0.415}, {0.838, 0.445}, {0.835, 0.252}, {0.679, 0.415}},
            {{0.632, 0.718}, {0.749, 0.943}, {0.924, 0.688}, {0.632, 0.718}}
        },
        {
            {{0, 0}, {2, 0}, {2, 2}, {0, 2}, {0, 0}},
            {{0.5, 0.05}, {0.5, 0.5}, {0.9, 0.5}, {0.9, 0.05}, {0.5, 0.05}},
            {{0.535, 0.688}, {0.596, 0.960}, {0.945, 0.967}, {0.968, 0.625}, {0.535, 0.688}}
        }
    );*/
    /*auto res = op_booleana_poligonos(
        {
            {{0.679, 0.415}, {0.838, 0.445}, {0.835, 0.252}, {0.679, 0.415}}
        },
        {
            {{0.5, 0.05}, {0.5, 0.5}, {0.9, 0.5}, {0.9, 0.05}, {0.5, 0.05}}
        }
    );*/
    // for (auto poligono : res) {
    //     std::cout << "poligono " << meeeeeeee << std::endl;
    //     for (auto parte : poligono) {
    //         std::cout << "    parte " << meu << std::endl;
    //         for (auto ponto : parte) {
    //             std::cout << "        " << ponto[0] << ' ' << ponto[1] << std::endl;
    //         }
    //         ++meu;
    //     }
    //     ++meeeeeeee;
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
    glBufferData(GL_ARRAY_BUFFER, 16*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
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
    glBufferData(GL_ARRAY_BUFFER, 16*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
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
    glBufferData(GL_ARRAY_BUFFER, 16*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
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
    glBufferData(GL_ARRAY_BUFFER, 16*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
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
    std::array<bool, 2> polys_pronto = {false, false};
    std::array<std::size_t, 2> ultimo_polys_idx = {0, 0};
    std::array<std::size_t, 2> ultimo_polys_pos = {0, 0};
    std::array<std::size_t, 2> polys_buffer_pos = {0, 0};
    std::array<std::vector<std::size_t>, 2> polys_indices_inicio;
    std::array<std::vector<std::size_t>, 2> polys_indices_fim;

    // estado.polys[0].push_back({});
    // estado.polys[1].push_back({});
    estado.mostrando_polys[0] = true;
    estado.mostrando_polys[1] = true;

    std::array<unsigned, 2> polys_vbo {};
    glGenBuffers(1, &polys_vbo[0]);
    glBindBuffer(GL_ARRAY_BUFFER, polys_vbo[0]);
    glBufferData(GL_ARRAY_BUFFER, 16*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    std::array<unsigned, 2> polys_vao {};
    glGenVertexArrays(1, &polys_vao[0]);
    glBindVertexArray(polys_vao[0]);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    
    glGenBuffers(1, &polys_vbo[1]);
    glBindBuffer(GL_ARRAY_BUFFER, polys_vbo[1]);
    glBufferData(GL_ARRAY_BUFFER, 16*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    glGenVertexArrays(1, &polys_vao[1]);
    glBindVertexArray(polys_vao[1]);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    bool inter_resul_pronto = false;
    unsigned inter_resul_vbo {};
    glGenBuffers(1, &inter_resul_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, inter_resul_vbo);
    glBufferData(GL_ARRAY_BUFFER, 16*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    unsigned inter_resul_vao {};
    glGenVertexArrays(1, &inter_resul_vao);
    glBindVertexArray(inter_resul_vao);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    CoisasDCEL coisas_dcel {};
    glGenVertexArrays(1, &coisas_dcel.vao);
    glBindVertexArray(coisas_dcel.vao);
    
    glGenBuffers(1, &coisas_dcel.vbo);
    glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
    glBufferData(GL_ARRAY_BUFFER, 16*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
    glEnableVertexAttribArray(1);

    glGenBuffers(1, &coisas_dcel.ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 16*1024*sizeof (unsigned), nullptr, GL_DYNAMIC_DRAW);
    glBindVertexArray(0);

    
    std::vector<std::size_t> inter_indices_inicio;
    std::vector<std::size_t> inter_indices_fim;

    std::size_t mano = 0;
    std::size_t atividade_size = 0;
    // bool visivel_pronto = false;
    std::vector<Ponto> area_visivel {};
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
                        Cor nova_cor = point_in_polygon(ponto, fecho_calculado) ? Cor::DENTRO : Cor::FORA;
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
            glClearColor(0.34f, 0.34f, 0.34f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            

            for (std::size_t poly_sel = 0; poly_sel < 2; ++poly_sel) {

                if (estado.limpar_ultimo_polys[poly_sel]) {
                    if (estado.polys_prontos[poly_sel] == 0) {
                        estado.limpar_polys[poly_sel] = true;
                    } else {
                        if (ultimo_polys_idx[poly_sel] == estado.polys_prontos[poly_sel]) {
                            ultimo_polys_pos[poly_sel] = 0;
                        }
                        estado.polys[poly_sel].pop_back();
                        // estado.polys[poly_sel].push_back({});
                        polys_buffer_pos[poly_sel] = polys_indices_fim[poly_sel][polys_indices_fim[poly_sel].size()-2];
                        polys_indices_inicio[poly_sel].back() = 0;
                        polys_indices_fim[poly_sel].back() = 0;
                    }
                    estado.limpar_ultimo_polys[poly_sel] = false;
                }
                if (estado.limpar_polys[poly_sel]) {
                    estado.polys[poly_sel] = {{}};
                    estado.mostrando_polys[poly_sel] = false;
                    estado.polys_prontos[poly_sel] = 0;
                    ultimo_polys_idx[poly_sel] = 0;
                    ultimo_polys_pos[poly_sel] = 0;
                    polys_indices_inicio[poly_sel].clear();
                    polys_indices_fim[poly_sel].clear();
                    polys_pronto[poly_sel] = false;
                    polys_buffer_pos[poly_sel] = 0;

                    estado.limpar_polys[poly_sel] = false;
                }
                std::size_t max_idx = estado.polys_prontos[poly_sel];
                std::size_t tamanho = 0;
                if (estado.polys[poly_sel].size() == max_idx) {
                    --max_idx;
                } else {
                    tamanho = estado.polys[poly_sel][max_idx].size();
                }
                if (estado.polys_prontos[poly_sel] > ultimo_polys_idx[poly_sel] || tamanho > ultimo_polys_pos[poly_sel]) {
                    // std::size_t diff = estado.polys[poly_sel][estado.polys_idx[poly_sel]].size() - ultimo_polys_pos[poly_sel];
                    std::vector<float> ps {};
                    // ps.reserve(diff * 5 * sizeof (float));
                    while (polys_indices_inicio[poly_sel].size() > estado.polys[poly_sel].size()) {
                        polys_indices_inicio[poly_sel].pop_back();
                        polys_indices_fim[poly_sel].pop_back();
                        polys_buffer_pos[poly_sel] = polys_indices_fim[poly_sel].back();
                    }
                    std::size_t old_buffer_pos = polys_buffer_pos[poly_sel];
                    while (polys_indices_inicio[poly_sel].size() < estado.polys[poly_sel].size()) {
                        polys_indices_inicio[poly_sel].push_back(0);
                        polys_indices_fim[poly_sel].push_back(0);
                    }
                    
                    if (estado.polys_prontos[poly_sel] > ultimo_polys_idx[poly_sel]) {
                        for (std::size_t i = ultimo_polys_idx[poly_sel]; i < estado.polys_prontos[poly_sel]; ++i) {
                            for (std::size_t j = ultimo_polys_pos[poly_sel]; j < estado.polys[poly_sel][i].size(); ++j) {
                                if (j == 0) {
                                    polys_indices_inicio[poly_sel][i] = polys_buffer_pos[poly_sel];
                                }
                                ++polys_buffer_pos[poly_sel];
                                auto ponto = estado.polys[poly_sel][i][j];
                                ps.push_back(ponto[0]);
                                ps.push_back(ponto[1]);
                                
                                if (poly_sel == 0) {
                                    if (i == 0) {
                                        ps.push_back(0.788f);
                                        ps.push_back(0.682f);
                                        ps.push_back(0.078f);
                                    } else {
                                        ps.push_back(0.325f); // 83
                                        ps.push_back(0.788f); // 201
                                        ps.push_back(0.078f); // 20
                                    }
                                } else {
                                    if (i == 0) {
                                        ps.push_back(0.078f); // 20
                                        ps.push_back(0.325f); // 83
                                        ps.push_back(0.788f); // 201
                                    } else {
                                        ps.push_back(0.682f);
                                        ps.push_back(0.078f);
                                        ps.push_back(0.788f);
                                    }
                                }
                            }
                            polys_indices_fim[poly_sel][i] = polys_buffer_pos[poly_sel];
                            ultimo_polys_pos[poly_sel] = 0;
                        }
                    }

                    for (std::size_t i = ultimo_polys_pos[poly_sel]; i < tamanho; ++i) {
                        // std::cout << i << ' ' << tamanho << std::endl;
                        if (i == 0) {
                            polys_indices_inicio[poly_sel][estado.polys_prontos[poly_sel]] = polys_buffer_pos[poly_sel];
                        }
                        ++polys_buffer_pos[poly_sel];
                        auto ponto = estado.polys[poly_sel][estado.polys_prontos[poly_sel]][i];
                        ps.push_back(ponto[0]);
                        ps.push_back(ponto[1]);
                        
                        if (poly_sel == 0) {
                            if (estado.polys_prontos[poly_sel] == 0) {
                                ps.push_back(0.788f);
                                ps.push_back(0.682f);
                                ps.push_back(0.078f);
                            } else {
                                ps.push_back(0.325f); // 83
                                ps.push_back(0.788f); // 201
                                ps.push_back(0.078f); // 20
                            }
                        } else {
                            if (estado.polys_prontos[poly_sel] == 0) {
                                ps.push_back(0.078f); // 20
                                ps.push_back(0.325f); // 83
                                ps.push_back(0.788f); // 201
                            } else {
                                ps.push_back(0.682f);
                                ps.push_back(0.078f);
                                ps.push_back(0.788f);
                            }
                        }
                    }
                    
                    polys_indices_fim[poly_sel].back() = polys_buffer_pos[poly_sel];

                    glBindBuffer(GL_ARRAY_BUFFER, polys_vbo[poly_sel]);
                    glBufferSubData(GL_ARRAY_BUFFER, static_cast<GLintptr>(old_buffer_pos * 5 * sizeof (float)), static_cast<GLintptr>(ps.size() * 5 * sizeof (float)), ps.data());
                    
                    ultimo_polys_pos[poly_sel] = tamanho;
                    ultimo_polys_idx[poly_sel] = estado.polys_prontos[poly_sel];
                    polys_pronto[poly_sel] = true;
                }
                if (polys_pronto[poly_sel] && estado.mostrando_polys[poly_sel]) {
                    
                    glBindBuffer(GL_ARRAY_BUFFER, polys_vbo[poly_sel]);
                    glBindVertexArray(polys_vao[poly_sel]);
                    std::size_t ate = polys_indices_fim[poly_sel].size();
                    // if (estado.polys[poly_sel][estado.polys_prontos[poly_sel]].size() == 0) {
                    //     ate--;
                    // }

                    color_line_program.use();
                    color_line_program.setFloat("alpha", 0.2f);
                    for (std::size_t i = 0; i < ate; ++i) {
                        glDrawArrays(GL_TRIANGLE_FAN, polys_indices_inicio[poly_sel][i], polys_indices_fim[poly_sel][i] - polys_indices_inicio[poly_sel][i]);
                    }
                    color_line_program.setFloat("alpha", 1.0f);
                    for (std::size_t i = 0; i < ate; ++i) {
                        glDrawArrays(GL_LINE_STRIP, polys_indices_inicio[poly_sel][i], polys_indices_fim[poly_sel][i] - polys_indices_inicio[poly_sel][i]);
                        // if (mano == 80 || mano == 0) {
                        //     std::cout << polys_indices_inicio[poly_sel][i] << ' ' << polys_indices_fim[poly_sel][i] << std::endl;
                        //     mano = 0;
                        // }
                    }
                    ++mano;

                    point_program.use();
                    point_program.setFloat("pointRadius", estado.pointSize);
                    glDrawArrays(GL_POINTS, 0, polys_buffer_pos[poly_sel]);
                }
            }

            if (estado.recalcular_intersecoes) {
                inter_indices_inicio.clear();
                inter_indices_fim.clear();
                // int meeeeeeee = 0;
                // int meu = 0;
                // for (auto poligono : std::array{estado.polys[0], estado.polys[1]}) {
                //     std::cout << "poligono " << meeeeeeee << std::endl;
                //     for (auto parte : poligono) {
                //         std::cout << "    parte " << meu << std::endl;
                //         for (auto ponto : parte) {
                //             std::cout << "        " << ponto[0] << ' ' << ponto[1] << std::endl;
                //         }
                //         ++meu;
                //     }
                //     ++meeeeeeee;
                // }
                estado.intersecoes = op_booleana_poligonos(estado.polys[0], estado.polys[1], false);
                int meeeeeeee = 0;
                int meu = 0;
                for (auto poligono : estado.intersecoes) {
                    std::cout << "poligono " << meeeeeeee << std::endl;
                    for (auto parte : poligono) {
                        std::cout << "    parte " << meu << std::endl;
                        for (auto ponto : parte) {
                            std::cout << "        " << ponto[0] << ' ' << ponto[1] << std::endl;
                        }
                        ++meu;
                    }
                    ++meeeeeeee;
                }
                // estado.intersecoes.pop_back();
                std::vector<float> ps {};
                for (std::size_t k = 0; k < estado.intersecoes.size(); ++k) {
                    for (std::size_t j = 0; j < estado.intersecoes[k].size(); ++j) {
                        inter_indices_inicio.push_back(ps.size()/5);
                        for (std::size_t i = 0; i < estado.intersecoes[k][j].size(); ++i) {
                            auto ponto = estado.intersecoes[k][j][i];
                            ps.push_back(ponto[0]);
                            ps.push_back(ponto[1]);
                            
                            if (j == 0) {
                                ps.push_back(0.788f); // 201
                                ps.push_back(0.325f); // 83
                                ps.push_back(0.078f); // 20
                            } else {
                                ps.push_back(0.788f); // 201
                                ps.push_back(0.078f); // 20
                                ps.push_back(0.325f); // 83
                            }
                        }
                        inter_indices_fim.push_back(ps.size()/5);

                    }
                }
                if (estado.intersecoes.size() == 0) {
                    inter_resul_pronto = false;
                } else {
                    glBindBuffer(GL_ARRAY_BUFFER, inter_resul_vbo);
                    glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(ps.size() * 5 * sizeof (float)), ps.data());
                    inter_resul_pronto = true;
                    estado.mostrando_intersecoes = true;
                }
                estado.recalcular_intersecoes = false;
            }
            if (inter_resul_pronto && estado.mostrando_intersecoes) {
                glBindBuffer(GL_ARRAY_BUFFER, inter_resul_vbo);
                glBindVertexArray(inter_resul_vao);

                color_line_program.use();
                color_line_program.setFloat("alpha", 0.2f);
                for (std::size_t i = 0; i < inter_indices_fim.size(); ++i) {
                    glDrawArrays(GL_TRIANGLE_FAN, inter_indices_inicio[i], inter_indices_fim[i] - inter_indices_inicio[i]);
                }
                color_line_program.setFloat("alpha", 1.0f);
                for (std::size_t i = 0; i < inter_indices_fim.size(); ++i) {
                    glDrawArrays(GL_LINE_STRIP, inter_indices_inicio[i], inter_indices_fim[i] - inter_indices_inicio[i]);
                    // if (mano == 80 || mano == 0) {
                    //     std::cout << inter_indices_inicio[i] << ' ' << inter_indices_fim[i] << std::endl;
                    //     mano = 0;
                    // }
                }
                // ++mano;

                point_program.use();
                point_program.setFloat("pointRadius", estado.pointSize);
                glDrawArrays(GL_POINTS, 0, inter_indices_fim.back());
            }

        } else if (estado.tela == Tela::ATIVIDADE) {

            glClearColor(0.3f, 0.2f, 0.3f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            if (estado.resetar_pontos) {
                estado.entrada.clear();
                estado.cores_entrada.clear();
                atividade_size = 0;
                area_visivel = {};
                estado.visivel_pronto = false;
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
                if (orientado_antihorario(estado.entrada)) {
                    std::cout << "orientação anti-horária" << std::endl;
                } else {
                    std::cout << "orientação horária" << std::endl;
                }
                // auto& v = estado.entrada;
                // int curvas_a_esquerda = 0;
                // for (std::size_t i = 1; i < v.size() - 1; ++i) {
                //     auto& p1 = v[i-1];
                //     auto& p2 = v[i];
                //     auto& p3 = v[i+1];
                //     if (left(p1, p2, p3)) {
                //         ++curvas_a_esquerda;
                //     } else if (area_orientada(p1, p2, p3) != 0.) {
                //         --curvas_a_esquerda;
                //     }
                // }
                // auto& p1 = v[0];
                // auto& p2 = v[1];
                // auto& pn_2 = v[v.size() - 2];
                // auto& pn_1 = v[v.size() - 1];
                // if (left(pn_2, pn_1, p1)) {
                //     ++curvas_a_esquerda;
                // } else if (area_orientada(pn_2, pn_1, p1) != 0.) {
                //     --curvas_a_esquerda;
                // }
                // if (left(pn_1, p1, p2)) {
                //     ++curvas_a_esquerda;
                // } else if (area_orientada(pn_1, p1, p2) != 0.) {
                //     --curvas_a_esquerda;
                // }
                // if (curvas_a_esquerda > 0) {
                //     std::cout << "orientação anti-horária" << std::endl;
                // } else {
                //     std::cout << "orientação horária" << std::endl;
                // }
                estado.recalcular_orientacao = false;
            }

            if (estado.recalcular_convexidade_dos_vertices) {
                estado.recalcular_convexidade_dos_vertices = false;

                // auto& v = estado.entrada;
                // for (std::size_t i = 0; i < v.size(); ++i) {
                //     std::size_t prev = (i == 0) ? (v.size()-1) : (i-1);
                //     std::size_t prox = (i+1 >= v.size()) ? (0) : (i+1);
                //     auto& p1 = v[prev];
                //     auto& p2 = v[i];
                //     auto& p3 = v[prox];
                //     Cor nova_cor {};
                //     if (area_orientada(p1, p2, p3) >= 0.) {
                //         nova_cor = Cor::DENTRO;
                //     } else {
                //         nova_cor = Cor::FORA;
                //     }
                //     estado.cores_entrada[i] = nova_cor;
                // }

                // std::vector<float> ps {};
                // ps.reserve(atividade_size * 5 * sizeof (float));
                // for (std::size_t i = 0; i < atividade_size; ++i) {
                //     auto ponto = estado.entrada[i];
                //     auto cor = estado.cores_entrada[i];
                //     ps.push_back(ponto[0]);
                //     ps.push_back(ponto[1]);
                //     // agora não vai ter nenhum amarelo
                //     if (cor == Cor::DESCONHECIDO) {
                //         std::cerr << "não era para ter amarelo" << std::endl;
                //         ps.push_back(0.788f); // 201
                //         ps.push_back(0.682f); // 174
                //         ps.push_back(0.078f); // 20
                //     } else if (cor == Cor::DENTRO) {
                //         ps.push_back(0.325f); // 83
                //         ps.push_back(0.788f); // 201
                //         ps.push_back(0.078f); // 20
                //     } else {
                //         ps.push_back(0.788f); // 201
                //         ps.push_back(0.149f); // 38
                //         ps.push_back(0.078f); // 20
                //     }
                // }
                // glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
                // glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(atividade_size * 5 * sizeof (float)), ps.data());
                
                // estado.recalcular_convexidade_dos_vertices = false;
            }
            
            if (estado.recalcular_orelhas) {
                auto& v = estado.entrada;
                for (std::size_t i = 0; i < v.size(); ++i) {
                    Categoria nova_cat = categoriza_ponto(v, i);
                    estado.cores_entrada[i] = nova_cat;
                }

                std::vector<float> ps {};
                ps.reserve(atividade_size * 5 * sizeof (float));
                for (std::size_t i = 0; i < atividade_size; ++i) {
                    auto ponto = estado.entrada[i];
                    auto cor = estado.cores_entrada[i];
                    ps.push_back(ponto[0]);
                    ps.push_back(ponto[1]);
                    switch (cor) {
                        case Categoria::REGULAR:
                            ps.push_back(static_cast<float>(201) / 255.f); // 201
                            ps.push_back(static_cast<float>(174) / 255.f); // 174
                            ps.push_back(static_cast<float>(20) / 255.f); // 20
                            break;
                        case Categoria::START:
                            ps.push_back(static_cast<float>(174) / 255.f); // 174
                            ps.push_back(static_cast<float>(201) / 255.f); // 201
                            ps.push_back(static_cast<float>(20) / 255.f); // 20
                            break;
                        case Categoria::END:
                            ps.push_back(static_cast<float>(201) / 255.f); // 201
                            ps.push_back(static_cast<float>(56) / 255.f); // 56
                            ps.push_back(static_cast<float>(20) / 255.f); // 20
                            break;
                        case Categoria::SPLIT:
                            ps.push_back(static_cast<float>(20) / 255.f); // 20
                            ps.push_back(static_cast<float>(122) / 255.f); // 122
                            ps.push_back(static_cast<float>(201) / 255.f); // 201
                            break;
                        case Categoria::MERGE:
                            ps.push_back(static_cast<float>(20) / 255.f); // 20
                            ps.push_back(static_cast<float>(56) / 255.f); // 56
                            ps.push_back(static_cast<float>(201) / 255.f); // 201
                            break;
                        default:
                            ps.push_back(static_cast<float>(201) / 255.f); // 201
                            ps.push_back(static_cast<float>(174) / 255.f); // 174
                            ps.push_back(static_cast<float>(20) / 255.f); // 20
                    }
                }
                glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
                glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(atividade_size * 5 * sizeof (float)), ps.data());
                
                estado.recalcular_orelhas = false;
            }

            if (estado.recalcular_visivel) {
                area_visivel = regiao_visivel(estado.entrada, estado.observador);

                std::vector<float> ps {};
                ps.reserve(area_visivel.size() * 5 * sizeof (float));
                for (std::size_t i = 0; i < area_visivel.size(); ++i) {
                    auto ponto = area_visivel[i];
                    ps.push_back(ponto[0]);
                    ps.push_back(ponto[1]);
                    // sempre tem verde
                    ps.push_back(static_cast<float>(119) / 255.f); // 119
                    ps.push_back(static_cast<float>(201) / 255.f); // 201
                    ps.push_back(static_cast<float>(20) / 255.f); // 20
                }
                glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
                glBufferSubData(GL_ARRAY_BUFFER, static_cast<GLintptr>(512 * 5 * sizeof (float)), static_cast<GLintptr>(area_visivel.size() * 5 * sizeof (float)), ps.data());

                std::array<float, 5> obs {
                    static_cast<float>(estado.observador[0]),
                    static_cast<float>(estado.observador[1]),
                    static_cast<float>(20) / 255.f,  // 20
                    static_cast<float>(201) / 255.f, // 201
                    static_cast<float>(168) / 255.f  // 168
                };

                glBufferSubData(GL_ARRAY_BUFFER, static_cast<GLintptr>(511 * 5 * sizeof (float)), static_cast<GLintptr>(5 * sizeof (float)), obs.data());

                estado.visivel_pronto = true;
                estado.recalcular_visivel = false;
            }
            
            glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
            glBindVertexArray(atividade_vao);

            color_line_program.use();
            color_line_program.setFloat("alpha", 0.2f);
            glDrawArrays(GL_TRIANGLE_FAN, 0, atividade_size);
            color_line_program.setFloat("alpha", 1.0f);
            glDrawArrays(GL_LINE_LOOP, 0, atividade_size);

            point_program.use();
            point_program.setFloat("pointRadius", estado.pointSize);
            glDrawArrays(GL_POINTS, 0, atividade_size);

            if (estado.visivel_pronto) {
                glBindBuffer(GL_ARRAY_BUFFER, atividade_vbo);
                glBindVertexArray(atividade_vao);

                color_line_program.use();
                // color_line_program.setFloat("alpha", 0.2f);
                // glDrawArrays(GL_TRIANGLE_FAN, 512, area_visivel.size());
                color_line_program.setFloat("alpha", 1.0f);
                glDrawArrays(GL_LINE_LOOP, 512, area_visivel.size());

                point_program.use();
                glDrawArrays(GL_POINTS, 511, area_visivel.size());
            }

        } else if (estado.tela == Tela::DCEL_TESTE) {

            glClearColor(0.4f, 0.3f, 0.4f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            auto& estad = estado.estado_dcel_teste;
            if (estad.estado == Dcel_Data::RESETANDO) {
                estad.poly.clear();
                estad.operacoes.clear();
                estad.ponto_adicionado = false;
                estad.poligono_fechado = false;

                coisas_dcel.last_size = 0;
                coisas_dcel.edge_count = 0;
                coisas_dcel.dcel_ptr.reset();
                coisas_dcel.coisas_piscar = {};
                estad.estado = Dcel_Data::RECEBENDO;
            }

            if (estad.ponto_adicionado) {
                // antigamente era como na linha a seguir, mas são a mesma coisa
                // if (estad.poly.size() > coisas_dcel.last_size) {
                std::size_t diff = estad.poly.size() - coisas_dcel.last_size;
                
                std::vector<float> ps {};
                ps.reserve(diff * 5 * sizeof (float));

                std::vector<unsigned> is {};
                is.reserve(diff * 2 * sizeof (unsigned));
                for (std::size_t i = coisas_dcel.last_size; i < estad.poly.size(); ++i) {
                    auto ponto = estad.poly[i];
                    ps.push_back(ponto[0]);
                    ps.push_back(ponto[1]);
                    // sempre começa com amarelo
                    ps.push_back(0.788f);
                    ps.push_back(0.682f);
                    ps.push_back(0.078f);
                    
                    if (i >= 1) {
                        is.push_back(i-1);
                        is.push_back(i);
                    }
                }
                glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                glBufferSubData(GL_ARRAY_BUFFER, static_cast<GLintptr>(coisas_dcel.last_size * 5 * sizeof (float)), static_cast<GLintptr>(diff * 5 * sizeof (float)), ps.data());
                
                // atualiza arestas a serem desenhadas:
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLintptr>(coisas_dcel.edge_count * 2 * sizeof (unsigned)), static_cast<GLintptr>(is.size() * sizeof (unsigned)), is.data());
                
                // atualiza contagem de arestas
                std::size_t edge_diff = diff;
                if (coisas_dcel.last_size == 0) {
                    --edge_diff;
                }
                coisas_dcel.edge_count += edge_diff;

                coisas_dcel.last_size = estad.poly.size();
                estad.ponto_adicionado = false;
            }

            if (estad.poligono_fechado) {
                // adiciona a última aresta no EBO:
                std::array<unsigned, 2> ultima_aresta {static_cast<unsigned>(estad.poly.size() - 1), 0};
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLintptr>(coisas_dcel.edge_count * 2 * sizeof (unsigned)), static_cast<GLintptr>(2 * sizeof (unsigned)), ultima_aresta.data());
                
                // atualiza contagem de arestas
                ++coisas_dcel.edge_count;
                
                estad.poligono_fechado = false;
            }

            if (estad.estado == Dcel_Data::CRIANDO_DCEL) {
                coisas_dcel.dcel_ptr = std::make_unique<DCEL>(estad.poly);
                estad.estado = Dcel_Data::DCEL_PRONTA;
            }

            if (estad.estado == Dcel_Data::DCEL_PRONTA) {
                while (estad.operacoes.size() > 0) {
                    auto op = estad.operacoes.front();
                    estad.operacoes.pop_front();
                    if (op.op == Dcel_Op::ENCONTRAR_FACE) {
                        // auto b = coisas_dcel.dcel_ptr->arestas_de_uma_face(0);
                        // for (auto c : b) {
                        //     std::cout << c[0][0] << ',' << c[0][1] << " -> " << c[1][0] << ',' << c[1][1] << std::endl;
                        // }

                        std::cout << "ala ó" << std::endl;
                        auto a = coisas_dcel.dcel_ptr->qual_face(op.args.ponto);
                        std::cout << "Indice da face: " << a << std::endl;
                    } else if (op.op == Dcel_Op::PISCAR_FACE) {
                        auto a = coisas_dcel.dcel_ptr->qual_face(op.args.ponto);
                        auto b = coisas_dcel.dcel_ptr->arestas_de_uma_face(a);
                        coisas_dcel.coisas_piscar.arestas = b;
                        coisas_dcel.coisas_piscar.ticks = 0;
                        coisas_dcel.coisas_piscar.ticks_por_aresta = 10;
                        coisas_dcel.coisas_piscar.atual = 0;
                        estad.estado = Dcel_Data::PISCANDO;
                        break;
                    }
                }
            }

            if (estad.estado == Dcel_Data::PISCANDO) {
                auto& c = coisas_dcel.coisas_piscar;
                if (c.ticks == 0) {
                    estad.estado = Dcel_Data::DCEL_PRONTA;
                    coisas_dcel.coisas_piscar = {};
                }
            }

            glBindVertexArray(coisas_dcel.vao);
            glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

            color_line_program.use();
            // color_line_program.setFloat("alpha", 0.2f);
            // glDrawArrays(GL_TRIANGLE_FAN, 0, coisas_dcel.last_size);
            // color_line_program.setFloat("alpha", 1.0f);
            // glDrawArrays(GL_LINE_LOOP, 0, coisas_dcel.last_size);
            color_line_program.setFloat("alpha", 1.0f);
            glDrawElements(GL_LINES, coisas_dcel.edge_count*2, GL_UNSIGNED_INT, nullptr);

            point_program.use();
            point_program.setFloat("pointRadius", estado.pointSize);
            glDrawArrays(GL_POINTS, 0, coisas_dcel.last_size);
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