#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
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
};

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
    State& estado = *(static_cast<State*> (glfwGetWindowUserPointer(window)));
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
    // corajosamente vou supor que r_atual nunca será menor que l_atual
    if (!encontrado && l_atual <= r_atual) {
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

        return {Etapa::ETAPA_1, {fecho[i], fecho[m]}, false, {}, fecho[i], fecho[m], false, {}};
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

int main() {

    {
        // std::vector<Ponto> pontos = {{15, 15}, {8, -9}, {-2, 3}, {-13, -5}, {-17, 20}};
        std::vector<Ponto> pontos = {{15, 15}, {8, -9}, {-2, 3}, {-13, -5}, {-17, 20}, {-20, -10}, {-10, 10}, {25, -10}, {-10, -15}, {10, -20}};
        // std::vector<Ponto> pontos = {{2, 3}, {2, -3}, {-2, 3}, {-2, -3}};
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
    }

    //return 0;

    GLFW::Session session {};

    // tamanho padrão de tela: 800x600
    GLFW::Window::options opts;
    opts.version_major = 3;
    opts.version_minor = 3;
    opts.decorated = 1;
    opts.samples = 128;
    GLFW::Window win {session, "Geometria", opts};

    GLFWwindow* window = win.justGimmeTheWindow();

    glfwMakeContextCurrent(window);
    
    State estado {};
    glfwSetWindowUserPointer(window, &estado);
    
    // glfwSetCursorEnterCallback(window, cursor_enter_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    // glfwSetCursorPosCallback(window, cursor_position_callback);
    // glfwSetScrollCallback(window, scroll_callback);
    glfwSetScrollCallback(window, scroll_callback);

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
    unsigned passo_vbo {};
    glGenBuffers(1, &passo_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
    glBufferData(GL_ARRAY_BUFFER, 2*1024*sizeof (float), nullptr, GL_DYNAMIC_DRAW);
    
    unsigned passo_vao {};
    glGenVertexArrays(1, &passo_vao);
    glBindVertexArray(passo_vao);
    
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

    RetornoAlg resultado_ate_agora {};
    bool resultado_arrumado_para_renderizacao = false;
    std::size_t passo_quantos {};
    while (!glfwWindowShouldClose(window)) {
        // win.processInput();
        // processar entradas??

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
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
            resultado_ate_agora = {};
        }
        if (estado.proximo_passo) {
            estado.proximo_passo = false;
            auto res = algoritmo_v1_passo_a_passo(fecho_calculado);
            if (res.etapa_do_passo_executado == Etapa::ETAPA_2) {
                resultado_ate_agora = res.resultado_ate_agora;
            }
            if (res.acabou) {
                estado.passo_a_passo_em_andamento = false;
                estado.passo_a_passo_acabou_de_acabar = true;
                resultado_arrumado_para_renderizacao = false;
            }
            std::array<Ponto, 6> pontos {};
            std::size_t num = 4;
            pontos[0] = res.colorir_esse;
            pontos[1] = res.esse_tambem;
            pontos[2] = res.desenhar_essa[0];
            pontos[3] = res.desenhar_essa[1];
            if (res.desenhar_a_outra) {
                num += 2;
                pontos[4] = res.tambem_desenhar_essa[0];
                pontos[5] = res.tambem_desenhar_essa[1];
            }

            std::array<std::array<float, 3>, 6> cores {};
            if (res.etapa_do_passo_executado == Etapa::ETAPA_1) {
                cores[0] = {};
                cores[1] = {};
                cores[2] = {};
                cores[3] = {};
                if (res.desenhar_a_outra) {
                    cores[4] = {};
                    cores[5] = {};
                }
            } else {
                cores[0] = {};
                cores[1] = {};
                cores[2] = {};
                cores[3] = {};
                if (res.desenhar_a_outra) {
                    cores[4] = {};
                    cores[5] = {};
                }
            }
            std::vector<float> ps {};
            ps.reserve(num * 5 * sizeof (float));
            for (std::size_t i = 0; i < num; ++i) {
                ps.push_back(pontos[i][0]);
                ps.push_back(pontos[i][1]);
                if (res.etapa_do_passo_executado == Etapa::ETAPA_1) {
                    ps.push_back(0.149f); // 38
                    ps.push_back(0.149f); // 38
                    ps.push_back(0.788f); // 201
                } else {
                    ps.push_back(0.149f); // 38
                    ps.push_back(0.682f); // 174
                    ps.push_back(0.788f); // 201
                }
            }
            glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
            glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(num * 5 * sizeof (float)), ps.data());

            glBindVertexArray(passo_vao);
            point_program.use();
            glDrawArrays(GL_POINTS, 0, 2);
            color_line_program.use();
            passo_quantos = num - 2;
            glDrawArrays(GL_LINES, 2, passo_quantos);
        }
        if (estado.passo_a_passo_em_andamento || estado.passo_a_passo_acabou_de_acabar) {
            glBindBuffer(GL_ARRAY_BUFFER, passo_vbo);
            glBindVertexArray(passo_vao);
            point_program.use();
            glDrawArrays(GL_POINTS, 0, 2);
            color_line_program.use();
            glDrawArrays(GL_LINES, 2, passo_quantos);
        }
        if (estado.mostrar_resultado_passo_a_passo) {
            estado.mostrando_resultado_passo_a_passo = true;
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
            point_program.use();
            glDrawArrays(GL_POINTS, 0, 4);
            color_line_program.use();
            glDrawArrays(GL_LINES, 0, 4);
        } else {
            estado.mostrando_resultado_passo_a_passo = false;
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