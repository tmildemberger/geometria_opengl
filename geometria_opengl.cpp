#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <unordered_set>
#include <deque>
#include <algorithm>
#include <memory>
#include <random>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <functional>

#include <limits>
#include <glad/glad.h>

#include "shader.hpp"
#include "Window.hpp"

#define STB_IMAGE_IMPLEMENTATION
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-overflow"
#pragma GCC diagnostic ignored "-Wswitch-default"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wmissing-declarations"
#include "stb_image.h"
#pragma GCC diagnostic pop

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

struct Ponto {
    double x;
    double y;
    bool operator==(const Ponto& rhs) const {
        return this->x == rhs.x && this->y == rhs.y;
    }
    bool operator<(const Ponto& rhs) const {
        if (this->x < rhs.x) return true;
        else if (this->x == rhs.x) return this->y > rhs.y;
        else return false;
    }
};

double area_orientada(Ponto p1, Ponto p2, Ponto p3);
bool left(Ponto p1, Ponto p2, Ponto p3);
std::vector<Ponto> fecho_convexo(std::vector<Ponto> pontos);
double dist(Ponto p1, Ponto p2);
double angulo_interno(Ponto p1, Ponto p2, Ponto p3);

double area_orientada(Ponto p1, Ponto p2, Ponto p3) {

    if (p1 == p2 || p1 == p3 || p2 == p3) return 0.0;
    return (p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y);
}

bool left(Ponto p1, Ponto p2, Ponto p3) {
    return area_orientada(p1, p2, p3) > 0.;
}

double dist(Ponto p1, Ponto p2) {
    return std::sqrt((p2.x - p1.x)*(p2.x - p1.x) + (p2.y - p1.y)*(p2.y - p1.y));
}

using Reta = std::array<Ponto, 2>;

double sombra_reta_ponto(Ponto p, Reta r);

double sombra_reta_ponto(Ponto p, Reta r) {
    double x3_x1 = p.x - r[0].x;
    double x2_x1 = r[1].x - r[0].x;
    double y3_y1 = p.y - r[0].y;
    double y2_y1 = r[1].y - r[0].y;

    double c = (x2_x1*x3_x1 + y2_y1*y3_y1) / (x2_x1*x2_x1 + y2_y1*y2_y1);
    return c;
}

class Cor {
public:
    Cor(std::string hex) {
        std::size_t start = 0;
        if (hex[0] == '#') ++start;

        r_ = std::stoul(hex.substr(start, 2), nullptr, 16);
        g_ = std::stoul(hex.substr(start + 2, 2), nullptr, 16);
        b_ = std::stoul(hex.substr(start + 4, 2), nullptr, 16);
    }
    Cor() {
        r_ = 0;
        g_ = 0;
        b_ = 0;
    }
    Cor(unsigned char r, unsigned char g, unsigned char b) {
        r_ = r;
        g_ = g;
        b_ = b;
    }
    float r() const { return r_/255.f; }
    float g() const { return g_/255.f; }
    float b() const { return b_/255.f; }
private:
    unsigned long r_;
    unsigned long g_;
    unsigned long b_;
};

enum class Intersecao {
    PROPRIA,
    IMPROPRIA,
    NAO,
};

Intersecao intersecao_com_left(Ponto p1, Ponto p2, Ponto p3, Ponto p4);

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

double distancia_ponto_reta_com_area(Ponto p1, Ponto p2, Ponto p);

double distancia_ponto_reta_com_area(Ponto p1, Ponto p2, Ponto p) {
    double area = std::abs(area_orientada(p1, p2, p));
    double base = dist(p1, p2);
    return area / base;
}

double distancia_ponto_segmento(Ponto p1, Ponto p2, Ponto p);

double distancia_ponto_segmento(Ponto p1, Ponto p2, Ponto p) {
    double altura = distancia_ponto_reta_com_area(p1, p2, p);
    double dist_p1 = dist(p1, p);
    double dist_p2 = dist(p2, p);
    double base = dist(p1, p2);
    double dist_max = std::sqrt(base * base + altura * altura);
    if (dist_p1 > dist_max) {
        return dist_p2;
    } if (dist_p2 > dist_max) {
        return dist_p1;
    }
    return altura;
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

double in_circle(Ponto a, Ponto b, Ponto c, Ponto d);

double in_circle(Ponto a, Ponto b, Ponto c, Ponto d) {

    double m_11 = a.x;
    double m_12 = a.y;
    double m_13 = a.x * a.x + a.y * a.y;
    double m_14 = 1;
    double m_21 = b.x;
    double m_22 = b.y;
    double m_23 = b.x * b.x + b.y * b.y;
    double m_24 = 1;
    double m_31 = c.x;
    double m_32 = c.y;
    double m_33 = c.x * c.x + c.y * c.y;
    double m_34 = 1;
    double m_41 = d.x;
    double m_42 = d.y;
    double m_43 = d.x * d.x + d.y * d.y;
    double m_44 = 1;

    double res_1 = m_11 * (m_22 * m_33 * m_44 + m_23 * m_34 * m_42 + m_32 * m_43 * m_24 - m_24 * m_33 * m_42 - m_23 * m_32 * m_44 - m_34 * m_43 * m_22);
    double res_2 = m_12 * (m_21 * m_33 * m_44 + m_23 * m_34 * m_41 + m_31 * m_43 * m_24 - m_24 * m_33 * m_41 - m_23 * m_31 * m_44 - m_34 * m_43 * m_21);
    double res_3 = m_13 * (m_21 * m_32 * m_44 + m_22 * m_34 * m_41 + m_31 * m_42 * m_24 - m_24 * m_32 * m_41 - m_22 * m_31 * m_44 - m_34 * m_42 * m_21);
    double res_4 = m_14 * (m_21 * m_32 * m_43 + m_22 * m_33 * m_41 + m_31 * m_42 * m_23 - m_23 * m_32 * m_41 - m_22 * m_31 * m_43 - m_33 * m_42 * m_21);

    return res_1 - res_2 + res_3 - res_4;
}

using Aresta = std::array<Ponto, 2>;

class DCEL {
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
        std::size_t componente;
    };
    struct Vertex {
        Ponto xy;
        Edge* edge;
    };

    struct EnganaCompilador;

public:
    DCEL(std::vector<Ponto> poligono_simples) : geracao_atual{0} {
        auto& p = poligono_simples;
        std::size_t n = p.size();

        if (p.back() == p[0]) {
            --n;
        }
        if (n < 3) {

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
            edges.push_back({(edges.data()) + (2*i) + 1, (vertices.data()) + i, nullptr, nullptr, inside_face, 0});
            edges.push_back({(edges.data()) + (2*i), (vertices.data()) + i + 1, nullptr, nullptr, outside_face, 0});
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

        vertice_valido = &vertices[0];

        geracao_atual = 1;
    }

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

            std::cerr << "erro em 'arestas_de_uma_face'" << std::endl;
        }
        return retorno;
    }

    std::vector<std::size_t> indices_das_arestas_de_uma_face(std::size_t face) {
        if (face >= faces.size()) {
            return {};
        }
        Edge* e = faces[face].edge;

        Edge* start = e;
        std::vector<std::size_t> retorno;
        retorno.push_back(static_cast<std::size_t>(e - edges.data()));
        e = e->next;

        while (e->face == &faces[face] && e != start) {
            retorno.push_back(static_cast<std::size_t>(e - edges.data()));
            e = e->next;
        }
        if (e->face != &faces[face]) {

            std::cerr << "erro em 'arestas_de_uma_face'" << std::endl;
        }
        return retorno;
    }

    std::vector<std::size_t> indices_dos_vertices_de_uma_face(std::size_t face) {
        if (face >= faces.size()) {
            return {};
        }
        Edge* e = faces[face].edge;

        Edge* start = e;
        std::vector<std::size_t> retorno;
        retorno.push_back(static_cast<std::size_t>(e->origin - vertices.data()));
        e = e->next;

        while (e->face == &faces[face] && e != start) {
            retorno.push_back(static_cast<std::size_t>(e->origin - vertices.data()));
            e = e->next;
        }
        if (e->face != &faces[face]) {

            std::cerr << "erro em 'indices_dos_vertices_de_uma_face'" << std::endl;
        }
        return retorno;
    }

    std::vector<Aresta> orbita_de_um_vertice(std::size_t vertice) {
        if (vertice >= vertices.size()) {
            return {};
        }
        Edge* e = vertices[vertice].edge;
        Edge* start = e;
        Edge* last = nullptr;
        e = e->prev->twin;
        std::vector<Aresta> retorno;
        while (e->origin == &vertices[vertice] && last != start) {
            last = e->twin->next;
            retorno.push_back(Aresta{e->origin->xy, e->twin->origin->xy});
            last = e;
            e = e->prev->twin;
        }
        if (e->origin != &vertices[vertice]) {

            std::cerr << "erro em 'orbita_de_um_vertice'" << std::endl;
        }
        return retorno;
    }

    std::vector<std::size_t> indices_orbita_de_um_vertice(std::size_t vertice) {
        if (vertice >= vertices.size()) {
            return {};
        }
        std::vector<std::size_t> retorno;
        Edge* e = vertices[vertice].edge;
        retorno.push_back(static_cast<std::size_t>(e - edges.data()));
        Edge* start = e;

        e = e->prev->twin;
        while (e->origin == &vertices[vertice] && e != start) {

            retorno.push_back(static_cast<std::size_t>(e - edges.data()));

            e = e->prev->twin;
        }
        if (e->origin != &vertices[vertice]) {

            std::cerr << "erro em 'orbita_de_um_vertice'" << std::endl;
        }
        return retorno;
    }

    void inclui_aresta(std::size_t v1_i, std::size_t v2_i) {
        if (v1_i >= vertices.size() || v2_i >= vertices.size() || v1_i == v2_i) {
            return;
        }

        reserva_espacos(1, 2, 0);

        Vertex* v1 = &vertices[v1_i];
        Vertex* v2 = &vertices[v2_i];

        Edge* e = v1->edge;
        Edge* start = e;
        Edge* last = nullptr;
        if (e->twin->origin == v2) {

            return;
        }
        e = e->prev->twin;
        Edge* found = nullptr;
        while (last != start) {
            last = e->twin->next;
            if (e->twin->origin == v2) {

                return;
            }
            if (!left(e->twin->origin->xy, v1->xy, last->twin->origin->xy)) {

                if (left(v1->xy, last->twin->origin->xy, v2->xy) || !left(v1->xy, e->twin->origin->xy, v2->xy)) {

                    found = last;
                    break;
                }
            } else {

                if (left(v1->xy, last->twin->origin->xy, v2->xy) && !left(v1->xy, e->twin->origin->xy, v2->xy)) {

                    found = last;
                    break;
                }
            }
            last = e;
            e = e->prev->twin;
        }
        if (!found) {

            std::cerr << "erro em 'inclui_aresta'" << std::endl;
            return;
        }
        start = found;
        e = start->prev;
        found = nullptr;
        while (e != start) {
            if (e->origin == v2) {

                found = e->prev;
            }
            if (intersecao_com_left(e->origin->xy, e->twin->origin->xy, v1->xy, v2->xy) != Intersecao::NAO) {
                if (e->origin != v1 && e->twin->origin != v1 && e->origin != v2 && e->twin->origin != v2) {
                    return;
                }
            }
            e = e->prev;
        }
        if (!found) {
            return;
        }

        std::cout << found->origin->xy.x << ',' << found->origin->xy.y << " -> " << found->twin->origin->xy.x << ',' << found->twin->origin->xy.y << std::endl;
        ++geracao_atual;
        faces.push_back({nullptr});

        std::size_t idx = edges.size();
        edges.push_back({(edges.data()) + idx + 1, v2, (edges.data()) + idx + 1, (edges.data()) + idx + 1, found->face, 0});
        edges.push_back({(edges.data()) + idx, v1, (edges.data()) + idx, (edges.data()) + idx, found->face, 0});
        conecta_arestas(found, &edges[idx]);
        conecta_arestas(start->prev, &edges[idx + 1]);

        if (found->face != faces.data()) {
            edges[idx + 1].face = &faces.back();
            faces.back().edge = &edges[idx + 1];
            e = edges[idx + 1].next;
            while (e != &edges[idx + 1]) {
                e->face = &faces.back();
                e = e->next;
            }
        } else {
            std::cout << "k" << std::endl;
            long long curvas_a_esquerda = 0;
            e = edges[idx + 1].next;
            while (e != &edges[idx + 1]) {
                double area = area_orientada(e->prev->origin->xy, e->origin->xy, e->twin->origin->xy);
                if (area > 0.) {
                    ++curvas_a_esquerda;
                } else if (area < 0.) {
                    --curvas_a_esquerda;
                }
                e = e->next;
            }
            start = &edges[idx + 1];
            if (curvas_a_esquerda < 0) {
                start = &edges[idx];
            }
            start->face = &faces.back();
            faces.back().edge = start;
            e = start->next;
            while (e != start) {
                e->face = &faces.back();
                e = e->next;
            }
        }

        edges[idx].face->edge = &edges[idx];
        edges[idx + 1].face->edge = &edges[idx + 1];
    }

    void inclui_vertice_em_aresta(std::size_t aresta, double onde) {
        if (onde <= 0 || onde >= 1 || aresta >= edges.size()) {
            return;
        }

        reserva_espacos(0, 2, 1);

        Edge* e = &edges[aresta];
        Vertex* v1 = e->origin;
        Vertex* v2 = e->twin->origin;
        Ponto p = {(v2->xy.x-v1->xy.x)*onde+v1->xy.x, (v2->xy.y-v1->xy.y)*onde+v1->xy.y};

        ++geracao_atual;

        vertices.push_back({p, e});
        Vertex* v3 = &vertices.back();
        Edge* e8 = e->prev;
        Edge* e9 = e->twin->next;

        e->origin = v3;

        std::size_t idx = edges.size();
        edges.push_back({(edges.data()) + idx + 1, v1, e, (edges.data()) + idx + 1, e->face, 0});
        edges.push_back({(edges.data()) + idx, v3, (edges.data()) + idx, e->twin, e->twin->face, 0});

        Edge* e2 = &edges[idx];

        v1->edge = e2;
        if (e->prev != e->twin) {
            edges[idx].prev = e8;
            edges[idx+1].next = e9;

            e8->next = e2;
            e9->prev = e2->twin;
        }

        e->prev = e2;
        e->twin->next = e2->twin;

    }

    void deleta_aresta(std::size_t aresta, bool atualiza_geracao = true) {
        if (aresta >= edges.size() || edges_invalidas.count(aresta)) {
            return;
        }

        Edge* e = &edges[aresta];
        Vertex* v1 = e->origin;
        Vertex* v2 = e->twin->origin;
        Edge* e_next = e->next;
        Edge* e_prev = e->prev;
        interno_deleta_aresta(aresta, false);
        while (e_next->prev == e_next->twin && !vazia()) {
            e_next = e_next->next;
            interno_deleta_aresta(static_cast<std::size_t>(e_next->prev - edges.data()), false);
        }
        while (e_prev->next == e_prev->twin && !vazia()) {
            e_prev = e_prev->prev;
            interno_deleta_aresta(static_cast<std::size_t>(e_prev->next - edges.data()), false);
        }

        if (atualiza_geracao) {
            ++geracao_atual;
        }
    }

    void deleta_vertice(std::size_t vertice, bool atualiza_geracao = true) {
        if (vertice >= vertices.size() || vertices_invalidas.count(vertice)) {
            return;
        }

        auto arestas = indices_orbita_de_um_vertice(vertice);
        for (auto aresta : arestas) {
            deleta_aresta(aresta, false);
        }

        if (atualiza_geracao) {
            ++geracao_atual;
        }
    }

    std::size_t qual_face(Ponto p) {

        Vertex* v1 = vertice_valido;
        if (!vertice_valido) {

            return 0;
        }

        Face* found_face = nullptr;
        
        std::size_t tentativas = 80000;
        while (!found_face) {
            if (tentativas == 0) {
                return 123456789;
            } else {
                --tentativas;
            }
            Edge* e = v1->edge;
            Edge* start = e;
            Edge* last = nullptr;
            e = e->prev->twin;
            Edge* found = nullptr;
            while (last != start) {
                last = e->twin->next;
                if (!left(e->twin->origin->xy, v1->xy, last->twin->origin->xy)) {

                    if (left(v1->xy, last->twin->origin->xy, p) || !left(v1->xy, e->twin->origin->xy, p)) {

                        found = last;
                        break;
                    }
                } else {

                    if (left(v1->xy, last->twin->origin->xy, p) && !left(v1->xy, e->twin->origin->xy, p)) {

                        found = last;
                        break;
                    }
                }
                last = e;
                e = e->prev->twin;
            }
            if (!found) {

                std::cerr << "erro1 em 'qual_face'" << std::endl;
                return 0;
            }

            Edge* a = found->next;

            while (true) {
                if (intersecao_com_left(a->origin->xy, a->twin->origin->xy, v1->xy, p) != Intersecao::NAO) {
                    if (a->origin != v1 && a->twin->origin != v1) {
                        break;
                    }
                }
                a = a->next;
                if (a == found || a->next == found) {
                    found_face = found->face;
                    break;
                }
            }
            v1 = a->twin->origin;

        }
        return static_cast<std::size_t>(found_face - faces.data());
    }

    std::pair<bool, std::size_t> em_alguma_aresta(Ponto p) {
        for (std::size_t i = 0; i < edges.size(); i += 2) {
            if (edges_invalidas.count(i)) {
                continue;
            }
            auto& p2 = edges[i].origin->xy;
            auto& p3 = edges[i+1].origin->xy;
            if (area_orientada(p, p2, p3) == 0.0 && dist(p, p2) + dist(p, p3) <= dist(p2, p3)) {
                return {true, i};
            }
        }
        return {false, 0};
    }

    std::pair<std::reference_wrapper<const std::vector<Vertex>>, const std::unordered_set<std::size_t>> vec_vertices() {
        return std::make_pair(std::cref(vertices), vertices_invalidas);
    }

    std::pair<std::reference_wrapper<const std::vector<Edge>>, const std::unordered_set<std::size_t>> vec_edges() {
        return std::make_pair(std::cref(edges), edges_invalidas);
    }

    std::pair<std::reference_wrapper<const std::vector<Face>>, const std::unordered_set<std::size_t>> vec_faces() {
        return std::make_pair(std::cref(faces), faces_invalidas);
    }

    std::size_t gen() {
        return geracao_atual;
    }

    std::size_t indice_vertice_valido() {
        return static_cast<std::size_t>(vertice_valido - vertices.data());
    }

    bool vazia() {
        return vertices.size() == 0 || vertices.size() == vertices_invalidas.size();
    }

    DCEL(EnganaCompilador engana_compilador, std::vector<Ponto> pontos) : geracao_atual{0} {

        if (engana_compilador) {
            geracao_atual = 0;
        }
        std::size_t n = pontos.size();

        vertices.reserve(20*n);
        edges.reserve(120*n);
        faces.reserve(120*n);
        faces.push_back({nullptr});
        Face* outside_face = (faces.data());
        for (std::size_t i = 0; i < n; ++i) {
            vertices.push_back({pontos[i], nullptr});
        }

        vertice_valido = &vertices[0];
    }
private:

    friend class CoisasDelaunay;
    friend class DelaunayPassoAPasso;

    friend class CoisasTrabalho;

    struct EnganaCompilador {
        explicit EnganaCompilador() = default;
        explicit operator bool() const { return true; }
    };
    bool novo_inclui_aresta(std::size_t v1_i, std::size_t v2_i) {
        if (v1_i >= vertices.size() || v2_i >= vertices.size() || v1_i == v2_i) {
            return false;
        }
        std::size_t menor = std::min(v1_i, v2_i);
        reserva_espacos(1, 2, 0);

        Vertex* v1 = &vertices[v1_i];
        Vertex* v2 = &vertices[v2_i];

        Edge* e = v1->edge;
        Edge* start = e;
        Edge* last = nullptr;
        Edge* found_v1 = nullptr;
        if (!e) {
        } else {
            if (e->twin->origin == v2) {

                return false;
            }
            e = e->prev->twin;
            while (last != start) {
                last = e->twin->next;
                if (e->twin->origin == v2) {

                    return false;
                }
                if (!left(e->twin->origin->xy, v1->xy, last->twin->origin->xy)) {

                    if (left(v1->xy, last->twin->origin->xy, v2->xy) || !left(v1->xy, e->twin->origin->xy, v2->xy)) {

                        found_v1 = last->prev;
                        break;
                    }
                } else {

                    if (left(v1->xy, last->twin->origin->xy, v2->xy) && !left(v1->xy, e->twin->origin->xy, v2->xy)) {

                        found_v1 = last->prev;
                        break;
                    }
                }
                last = e;
                e = e->prev->twin;
            }
            if (!found_v1) {

                std::cerr << "erro v1 em 'novo_inclui_aresta'" << std::endl;
                return false;
            }
        }

        e = v2->edge;
        start = e;
        last = nullptr;
        Edge* found_v2 = nullptr;
        if (!e) {
        } else {
            if (e->twin->origin == v1) {

                return false;
            }
            e = e->prev->twin;
            while (last != start) {
                last = e->twin->next;
                if (e->twin->origin == v1) {

                    return false;
                }
                if (!left(e->twin->origin->xy, v2->xy, last->twin->origin->xy)) {

                    if (left(v2->xy, last->twin->origin->xy, v1->xy) || !left(v2->xy, e->twin->origin->xy, v1->xy)) {

                        found_v2 = last->prev;
                        break;
                    }
                } else {

                    if (left(v2->xy, last->twin->origin->xy, v1->xy) && !left(v2->xy, e->twin->origin->xy, v1->xy)) {

                        found_v2 = last->prev;
                        break;
                    }
                }
                last = e;
                e = e->prev->twin;
            }
            if (!found_v2) {

                std::cerr << "erro v2 em 'novo_inclui_aresta'" << std::endl;
                return false;
            }
        }

        std::size_t idx = edges.size();
        edges.push_back({(edges.data()) + idx + 1, v2, (edges.data()) + idx + 1, (edges.data()) + idx + 1, faces.data(), menor});
        edges.push_back({(edges.data()) + idx, v1, (edges.data()) + idx, (edges.data()) + idx, faces.data(), menor});

        vertice_valido = v1;

        std::size_t v2_c = 0;
        std::size_t v1_c = 0;
        if (found_v2) {

            conecta_arestas(found_v2, &edges[idx]);
            edges[idx].face = found_v2->face;
            edges[idx + 1].face = found_v2->face;
            v2_c = found_v2->componente;
            if (v2_c < menor) {
                for (auto& ed : edges) {
                    if (ed.componente == menor) {
                        ed.componente = v2_c;
                    }
                }
                menor = v2_c;
            } else if (v2_c > menor) {
                for (auto& ed : edges) {
                    if (ed.componente == v2_c) {
                        ed.componente = menor;
                    }
                }
            }
        }
        v2->edge = &edges[idx];
        if (found_v1) {

            conecta_arestas(found_v1, &edges[idx + 1]);
            edges[idx + 1].face = found_v1->face;
            edges[idx].face = found_v1->face;
            v1_c = found_v1->componente;
            if (v1_c < menor) {
                for (auto& ed : edges) {
                    if (ed.componente == menor) {
                        ed.componente = v1_c;
                    }
                }
                menor = v1_c;
            } else if (v1_c > menor) {
                for (auto& ed : edges) {
                    if (ed.componente == v1_c) {
                        ed.componente = menor;
                    }
                }
            }
        }
        v1->edge = &edges[idx + 1];

        if (found_v1 && found_v2) {
            if (found_v1->face != faces.data()) {
                if (v1_c != v2_c) {
                    std::cout << "ta errado isso ai" << std::endl;
                }
                faces.push_back({nullptr});
                edges[idx + 1].face = &faces.back();
                faces.back().edge = &edges[idx + 1];
                e = edges[idx + 1].next;
                while (e != &edges[idx + 1]) {
                    e->face = &faces.back();
                    e = e->next;
                }
            } else {
                long long curvas_a_esquerda_e1 = 0;
                long long curvas_a_esquerda_e2 = 0;
                e = &edges[idx + 1];
                do {
                    double area = area_orientada(e->prev->origin->xy, e->origin->xy, e->twin->origin->xy);
                    if (area > 0.) {
                        ++curvas_a_esquerda_e2;
                    } else if (area < 0.) {
                        --curvas_a_esquerda_e2;
                    }
                    e = e->next;
                }while (e != &edges[idx + 1]);
                e = &edges[idx];
                do {
                    double area = area_orientada(e->prev->origin->xy, e->origin->xy, e->twin->origin->xy);
                    if (area > 0.) {
                        ++curvas_a_esquerda_e1;
                    } else if (area < 0.) {
                        --curvas_a_esquerda_e1;
                    }
                    e = e->next;
                } while (e != &edges[idx]);
                if (curvas_a_esquerda_e1 > 0 || curvas_a_esquerda_e2 > 0) {
                    if (v1_c != v2_c) {
                        std::cout << "ta muito errado isso ai" << std::endl;
                    }

                    faces.push_back({nullptr});

                    start = &edges[idx + 1];
                    if (curvas_a_esquerda_e2 < 0) {
                        start = &edges[idx];
                    }
                    start->face = &faces.back();
                    faces.back().edge = start;
                    e = start->next;
                    while (e != start) {
                        e->face = &faces.back();
                        e = e->next;
                    }
                } else {

                    if (v1_c == v2_c) {

                        std::cout << "ta extremamente errado isso ai" << std::endl;
                    }
                }
            }
        }
        edges[idx].face->edge = &edges[idx];
        edges[idx + 1].face->edge = &edges[idx + 1];
        return true;
    }

    std::vector<Face> faces;
    std::vector<Edge> edges;
    std::vector<Vertex> vertices;
    std::unordered_set<std::size_t> faces_invalidas;
    std::unordered_set<std::size_t> edges_invalidas;
    std::unordered_set<std::size_t> vertices_invalidas;

    Vertex* vertice_valido;

    std::size_t geracao_atual;

    void interno_deleta_aresta(std::size_t aresta, bool atualiza_geracao = true) {
        if (aresta >= edges.size() || edges_invalidas.count(aresta)) {
            return;
        }
        bool deleta_v1 = false;
        bool deleta_v2 = false;
        Edge* e = &edges[aresta];
        Vertex* v1 = e->origin;
        Vertex* v2 = e->twin->origin;
        Face* f = e->face;
        Face* f_outra = e->twin->face;
        Edge* f_aux = e;
        if (f_outra != f) {

            if (f == faces.data()) {
                f = e->twin->face;
                f_outra = e->face;
                f_aux = e->twin;
            }
            Edge* a = f_aux->next;
            while (a != f_aux) {
                a->face = f_outra;
                a = a->next;
            }
            faces_invalidas.insert(static_cast<std::size_t>(f - faces.data()));
        }

        if (e->next == e->twin) {
            deleta_v2 = true;
            vertices_invalidas.insert(static_cast<std::size_t>(v2 - vertices.data()));
        }

        if (e->prev == e->twin) {
            deleta_v1 = true;
            vertices_invalidas.insert(static_cast<std::size_t>(v1 - vertices.data()));
        }

        if (!deleta_v2) {
            e->next->prev = e->twin->prev;
            e->twin->prev->next = e->next;

            f_outra->edge = e->next;
            vertice_valido = v2;
            v2->edge = e->next;
        }

        if (!deleta_v1) {
            e->prev->next = e->twin->next;
            e->twin->next->prev = e->prev;

            f_outra->edge = e->prev;
            vertice_valido = v1;
            v1->edge = e->twin->next;
        }

        edges_invalidas.insert(static_cast<std::size_t>(e - edges.data()));
        edges_invalidas.insert(static_cast<std::size_t>(e->twin - edges.data()));

        if (deleta_v1 && deleta_v2) {
            if (f_outra->edge == e || f_outra->edge == e->twin) {
                Edge* found = nullptr;
                for (std::size_t i = 0; i < edges.size(); ++i) {
                    if (edges_invalidas.count(i)) {
                        continue;
                    }
                    if (edges[i].face == f_outra) {
                        found = &edges[i];
                        break;
                    }
                }
                if (!found) {

                    faces_invalidas.insert(static_cast<std::size_t>(f_outra - faces.data()));
                } else {
                    f_outra->edge = found;
                }
            }
            if (vertice_valido == v1 || vertice_valido == v2) {
                for (std::size_t i = 0; i < vertices.size(); ++i) {
                    if (vertices_invalidas.count(i)) {
                        continue;
                    }
                    vertice_valido = &vertices[i];
                    break;
                }
            }
        }

        if (atualiza_geracao) {
            ++geracao_atual;
        }
    }

    void interno_deleta_vertice(std::size_t vertice, bool atualiza_geracao = true) {
        if (vertice >= vertices.size() || vertices_invalidas.count(vertice)) {
            return;
        }

        auto arestas = indices_orbita_de_um_vertice(vertice);
        for (auto aresta : arestas) {
            interno_deleta_aresta(aresta, false);
        }

        if (atualiza_geracao) {
            ++geracao_atual;
        }
    }
    void reserva_espacos(std::size_t n_faces, std::size_t n_edges, std::size_t n_vertices) {
        std::size_t faces_cap = faces.capacity();
        std::size_t edges_cap = edges.capacity();
        std::size_t vertices_cap = vertices.capacity();

        if (n_faces) {
            Face* base = faces.data();
            faces.reserve(faces.size() + n_faces);
            if (faces.capacity() != faces_cap) {
                std::cout << "faaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << std::endl;

                for (auto& edge : edges) {
                    if (edge.face) { edge.face = (faces.data()) + (edge.face - base); }
                }
            }
        }

        if (n_vertices) {
            Vertex* base = vertices.data();
            vertices.reserve(vertices.size() + n_vertices);
            if (vertices.capacity() != vertices_cap) {
                std::cout << "vaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << std::endl;

                for (auto& edge : edges) {
                    if (edge.origin) { edge.origin = (vertices.data()) + (edge.origin - base); }
                }
                vertice_valido = (vertices.data()) + (vertice_valido - base);
            }
        }

        if (n_edges) {
            Edge* base = edges.data();
            edges.reserve(edges.size() + n_edges);

            if (edges.capacity() != edges_cap) {
                std::cout << "eaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << std::endl;

                for (auto& face : faces) {
                    if (face.edge) { face.edge = (edges.data()) + (face.edge - base); }
                }
                for (auto& vertex : vertices) {
                    if (vertex.edge) { vertex.edge = (edges.data()) + (vertex.edge - base); }
                }
                for (auto& edge : edges) {
                    if (edge.twin) { edge.twin = (edges.data()) + (edge.twin - base); }
                    if (edge.next) { edge.next = (edges.data()) + (edge.next - base); }
                    if (edge.prev) { edge.prev = (edges.data()) + (edge.prev - base); }
                }
            }
        }
    }

    void conecta_arestas(Edge* e1, Edge* e2) {

        e2->twin->next = e1->next;
        e2->prev = e1;
        e1->next->prev = e2->twin;
        e1->next = e2;
    }
};

enum class Tela {
    DCEL_TESTE,
    DELAUNAY,

};

enum class Dcel_Data {
    RESETANDO,
    RECEBENDO,
    CRIANDO_DCEL,
    DCEL_PRONTA,
    PISCANDO,
    ADICIONANDO_ARESTA,
    ADICIONANDO_VERTICE,
    ESPERANDO_ORBITA,
    DELETANDO_ARESTA,
    DELETANDO_VERTICE,

};

enum class Dcel_Op {
    ENCONTRAR_FACE,
    PISCAR_FACE,
    PONTO_SELECIONADO,
    CLIQUE_VERTICE,
    PISCAR_ORBITA,
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

enum class General_Op {
    CLIQUE,
    TECLA,
};

struct General_Msg {
    Ponto p;
    int button_key;
    int mods;
    General_Op op;
};

struct Delaunay_State {
    std::deque<General_Msg> eventos;
};

struct Trabalho_State {
    std::deque<General_Msg> eventos;
};

struct State {
    Dcel_Teste_State estado_dcel_teste;
    Delaunay_State estado_delaunay;
    Trabalho_State estado_trabalho;
    float pointSize;
    Tela tela;

};

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods);

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

    if (estado.tela == Tela::DCEL_TESTE) {
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
            } else if (button == GLFW_MOUSE_BUTTON_MIDDLE && !mods) {
                Ponto clicado = ponto_xy();
                estad.operacoes.push_back({{clicado}, Dcel_Op::CLIQUE_VERTICE});
            }
        } else if (estad.estado == Dcel_Data::ADICIONANDO_ARESTA) {
            if (button == GLFW_MOUSE_BUTTON_MIDDLE && !mods) {
                estad.estado = Dcel_Data::DCEL_PRONTA;
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && !mods) {
                Ponto clicado = ponto_xy();
                estad.operacoes.push_back({{clicado}, Dcel_Op::PONTO_SELECIONADO});
            }
        } else if (estad.estado == Dcel_Data::ADICIONANDO_VERTICE) {
            if (button == GLFW_MOUSE_BUTTON_MIDDLE && !mods) {
                estad.estado = Dcel_Data::DCEL_PRONTA;
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && !mods) {
                Ponto clicado = ponto_xy();
                estad.operacoes.push_back({{clicado}, Dcel_Op::CLIQUE_VERTICE});
            }
        } else if (estad.estado == Dcel_Data::ESPERANDO_ORBITA) {
            if (button == GLFW_MOUSE_BUTTON_MIDDLE && !mods) {
                estad.estado = Dcel_Data::DCEL_PRONTA;
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && !mods) {
                Ponto clicado = ponto_xy();
                estad.operacoes.push_back({{clicado}, Dcel_Op::PISCAR_ORBITA});
            }
        } else if (estad.estado == Dcel_Data::DELETANDO_ARESTA) {
            if (button == GLFW_MOUSE_BUTTON_MIDDLE && !mods) {
                estad.estado = Dcel_Data::DCEL_PRONTA;
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && !mods) {
                Ponto clicado = ponto_xy();
                estad.operacoes.push_back({{clicado}, Dcel_Op::CLIQUE_VERTICE});
            }
        } else if (estad.estado == Dcel_Data::DELETANDO_VERTICE) {
            if (button == GLFW_MOUSE_BUTTON_MIDDLE && !mods) {
                estad.estado = Dcel_Data::DCEL_PRONTA;
            } else if (button == GLFW_MOUSE_BUTTON_LEFT && !mods) {
                Ponto clicado = ponto_xy();
                estad.operacoes.push_back({{clicado}, Dcel_Op::CLIQUE_VERTICE});
            }
        }
    } else if (estado.tela == Tela::DELAUNAY) {
        if (action != GLFW_RELEASE) return;
        auto& estad = estado.estado_delaunay;
        Ponto clicado = ponto_xy();
        estad.eventos.push_back({clicado, button, mods, General_Op::CLIQUE});
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
    if (action == GLFW_RELEASE) {
        switch (key) {
            case GLFW_KEY_5:
                if (!mods) estado.tela = Tela::DCEL_TESTE;
                break;
            case GLFW_KEY_6:
                if (!mods) estado.tela = Tela::DELAUNAY;
                break;

            case GLFW_KEY_R:

                if (estado.tela == Tela::DELAUNAY) {
                    estado.estado_delaunay.eventos.push_back({{}, key, mods, General_Op::TECLA});
                }

                break;
            case GLFW_KEY_T:
                if (estado.tela == Tela::DCEL_TESTE && !mods) {
                    if (estado.estado_dcel_teste.estado == Dcel_Data::DCEL_PRONTA) {
                        estado.estado_dcel_teste.estado = Dcel_Data::DELETANDO_VERTICE;
                    }
                } else if (estado.tela == Tela::DELAUNAY) {
                    estado.estado_delaunay.eventos.push_back({{}, key, mods, General_Op::TECLA});
                }

                break;
            case GLFW_KEY_A:
                if (estado.tela == Tela::DCEL_TESTE && !mods) {
                    if (estado.estado_dcel_teste.estado == Dcel_Data::DCEL_PRONTA) {
                        estado.estado_dcel_teste.estado = Dcel_Data::ADICIONANDO_ARESTA;
                    }
                } else if (estado.tela == Tela::DELAUNAY) {
                    estado.estado_delaunay.eventos.push_back({{}, key, mods, General_Op::TECLA});
                }
                break;
            case GLFW_KEY_C:
                if (estado.tela == Tela::DCEL_TESTE && !mods) {
                    if (estado.estado_dcel_teste.estado == Dcel_Data::DCEL_PRONTA) {
                        estado.estado_dcel_teste.estado = Dcel_Data::DELETANDO_ARESTA;
                    }
                } else if (estado.tela == Tela::DELAUNAY) {
                    estado.estado_delaunay.eventos.push_back({{}, key, mods, General_Op::TECLA});
                }
                break;
            case GLFW_KEY_S:
                if (estado.tela != Tela::DELAUNAY) {
                    return;
                }
                estado.estado_delaunay.eventos.push_back({{}, key, mods, General_Op::TECLA});
                break;
            case GLFW_KEY_E:
                if (estado.tela != Tela::DELAUNAY) {
                    return;
                }
                estado.estado_delaunay.eventos.push_back({{}, key, mods, General_Op::TECLA});
                break;
            case GLFW_KEY_O:
                if (estado.tela != Tela::DCEL_TESTE || mods) {
                    return;
                }
                if (estado.estado_dcel_teste.estado == Dcel_Data::DCEL_PRONTA) {
                    estado.estado_dcel_teste.estado = Dcel_Data::ESPERANDO_ORBITA;
                }
                break;
            case GLFW_KEY_V:
                if (estado.tela != Tela::DCEL_TESTE || mods) {
                    return;
                }
                if (estado.estado_dcel_teste.estado == Dcel_Data::DCEL_PRONTA) {
                    estado.estado_dcel_teste.estado = Dcel_Data::ADICIONANDO_VERTICE;
                }
                break;
            default:
                break;
        }
    } else if (action == GLFW_REPEAT && key == GLFW_KEY_R) {

        if (estado.tela == Tela::DELAUNAY) {
            estado.estado_delaunay.eventos.push_back({{}, key, mods, General_Op::TECLA});
        }

    }
}

struct Coisas_Para_Piscar {
    std::size_t ticks;
    std::size_t ticks_por_aresta;
    std::size_t atual;
    std::vector<std::size_t> arestas;
};

struct Coisas_Para_Adicionar_Aresta {
    std::size_t p1_idx;
    std::size_t p2_idx;
    std::size_t indices_obtidos;
};

struct Coisas_Para_Adicionar_Vertice {
    std::size_t aresta_idx;
    bool aresta_selecionada;
};

struct CoisasDCEL {
    unsigned vao;
    unsigned vbo;
    unsigned ebo;
    unsigned extra_vao;
    unsigned extra_ebo;
    std::size_t last_size;
    std::size_t edge_count;
    std::size_t last_gen;
    std::unique_ptr<DCEL> dcel_ptr;
    Coisas_Para_Piscar coisas_piscar;
    Coisas_Para_Adicionar_Aresta coisas_aresta;
    Coisas_Para_Adicionar_Vertice coisas_vertice;
};

enum class EstadoDelaunay {
    INICIANDO,
    TRIANGULANDO,
    OK,
};

enum class EntradaDelaunay {
    NORMAL,
    TROCANDO_ARESTA,
};

bool aconteceu_aquilo = false;

struct CoisasDelaunay {
    CoisasDelaunay(std::string imagem) : x{0}, y{0}, n{0} {
        glGenBuffers(1, &extra_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, extra_vbo);
        glBufferData(GL_ARRAY_BUFFER, 128*sizeof (float), nullptr, GL_DYNAMIC_DRAW);

        glGenVertexArrays(1, &extra_vao);
        glBindVertexArray(extra_vao);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
        glEnableVertexAttribArray(1);
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, max_floats*sizeof (float), nullptr, GL_DYNAMIC_DRAW);

        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
        glEnableVertexAttribArray(1);

        glGenBuffers(1, &ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, max_floats*sizeof (unsigned), nullptr, GL_DYNAMIC_DRAW);

        glGenVertexArrays(1, &faces_vao);
        glBindVertexArray(faces_vao);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
        glEnableVertexAttribArray(1);

        glGenBuffers(1, &faces_ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, faces_ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, max_floats*sizeof (unsigned), nullptr, GL_DYNAMIC_DRAW);

        glBindVertexArray(0);

        mostrando_linhas = false;
        estado = EstadoDelaunay::OK;
        estado_entrada = EntradaDelaunay::NORMAL;
        last_size = 0;
        edge_count = 0;
        triangle_count = 0;
        last_gen = 0;

        image_data = stbi_load(imagem.c_str(), &x, &y, &n, 0);
        if (!image_data) {
            std::cerr << "nao carregou imagem \"" << imagem << "\"" << std::endl;
            std::exit(1);
        }
        if (n != 3 || x != y) {
            std::cout << "que estranho " << x << ' ' << y << ' ' << n << std::endl;
        }
        inicia();
    }
    ~CoisasDelaunay() {
        stbi_image_free(image_data);
    }

    void reset(std::string imagem) {

        stbi_image_free(image_data);
        mostrando_linhas = false;
        estado = EstadoDelaunay::OK;
        estado_entrada = EntradaDelaunay::NORMAL;
        last_size = 0;
        edge_count = 0;
        triangle_count = 0;
        last_gen = 0;
        dcel.reset();

        image_data = stbi_load(imagem.c_str(), &x, &y, &n, 0);
        if (!image_data) {
            std::cerr << "nao carregou imagem \"" << imagem << "\"" << std::endl;
            std::exit(1);
        }
        inicia();
    }

    Cor encontra_cor(Ponto p) {
        double p_x = std::floor(((p.x + 1.0) / 2.0) * x);
        int i_x = std::min(x, static_cast<int>(p_x));

        double p_y = std::floor(((p.y + 1.0) / 2.0) * y);
        int i_y = std::min(y, static_cast<int>(p_y));
        // std::cout << p.x << ' ' << p.y << " -- " << p_x << ' ' << p_y << std::endl;
        // std::cout << "foi buscada a cor do pixel " << i_x << ' ' << i_y << std::endl;

        unsigned char r = image_data[i_y * y * n + i_x * n + 0];
        unsigned char g = image_data[i_y * y * n + i_x * n + 1];
        unsigned char b = image_data[i_y * y * n + i_x * n + 2];

        return Cor(r, g, b);
    }

    bool adiciona_ponto(Ponto p) {
        double p_x = std::floor(((p.x + 1.0) / 2.0) * x);
        int i_x = std::min(x, static_cast<int>(p_x));

        double p_y = std::floor(((p.y + 1.0) / 2.0) * y);
        int i_y = std::min(y, static_cast<int>(p_y));
        std::cout << "foi adicionado o pixel " << i_x << ' ' << i_y << std::endl;

        double new_x = (((static_cast<double>(i_x)) / x) * 2.0) - 1.0 + (1.0 / x);
        double new_y = (((static_cast<double>(i_y)) / y) * 2.0) - 1.0 + (1.0 / y);

        double lim_inf = (((0) / x) * 2.0) - 1.0 + (1.0 / x);
        double lim_sup = ((static_cast<double>(x-1) / x) * 2.0) - 1.0 + (1.0 / x);
        if (new_x < lim_inf || new_x > lim_sup || new_y < lim_inf || new_y > lim_sup) {
            std::cout << "nao entendi " << new_x << ' ' << new_y << std::endl;
            return false;
        }
        Ponto new_p = {new_x, new_y};
        if (pontos.count(new_p)) {
            return false;
        }
        pontos.insert(new_p);

        triangulacao(new_p);
        std::cout << "supostamente triangulado" << std::endl;

        return true;
    }
private:
    void triangulacao(Ponto novo) {
        // adiciona novo ponto e arruma a triangulação
        // (algoritmo incremental)
        auto [resposta, qual_aresta] = dcel->em_alguma_aresta(novo);
        std::size_t face = 1234567890;
        if (resposta) {
            // std::cout << "s pra saber" << std::endl;
            // std::exit(1);
            auto& es = dcel->edges;
            auto& fs = dcel->faces;

            std::size_t proxima = static_cast<std::size_t>(es[qual_aresta].next - es.data());
            dcel->deleta_aresta(qual_aresta, false);

            face = static_cast<std::size_t>(es[proxima].face - fs.data());
        } else {
            face = dcel->qual_face(novo);
            if (face == 123456789) {
                // isso significa que houve um erro
                std::cout << "deu o problema com o ponto " << novo.x << ' ' << novo.y << std::endl;
                aconteceu_aquilo = true;
                return;
            }
        }
        std::cout << "face do novo ponto: " << face << std::endl;
        auto arestas = dcel->indices_das_arestas_de_uma_face(face);

        dcel->reserva_espacos(0, 0, 1);
        dcel->vertices.push_back(DCEL::Vertex{novo, nullptr});
        auto vertices = dcel->indices_dos_vertices_de_uma_face(face);
        for (auto v : vertices) {
            dcel->novo_inclui_aresta(dcel->vertices.size() - 1, v);
            std::cout << "incluida aresta entre " << dcel->vertices.size() - 1 << " e " << v << std::endl;
        }

        std::deque<std::size_t> lista_de_potencialmente_invalidas;
        for (auto aresta : arestas) {
            lista_de_potencialmente_invalidas.push_back(aresta);
        }

        while (lista_de_potencialmente_invalidas.size() > 0) {
            std::size_t aresta = lista_de_potencialmente_invalidas.front();
            lista_de_potencialmente_invalidas.pop_front();

            auto& vs = dcel->vertices;
            auto& es = dcel->edges;

            auto p3 = static_cast<std::size_t>(es[aresta].origin - vs.data());
            auto p4 = static_cast<std::size_t>(es[aresta].twin->origin - vs.data());

            auto p1 = static_cast<std::size_t>(es[aresta].prev->origin - vs.data());
            auto p2 = static_cast<std::size_t>(es[aresta].twin->prev->origin - vs.data());
            
            std::cout << "aresta " << aresta;
            if (intersecao_com_left(vs[p1].xy, vs[p2].xy, vs[p3].xy, vs[p4].xy) == Intersecao::PROPRIA) {
                if (in_circle(vs[p3].xy, vs[p4].xy, vs[p1].xy, vs[p2].xy) > 0.0) {
                    std::cout << " precisou ser flipada" << std::endl;
                    auto e1 = static_cast<std::size_t>(es[aresta].twin->next - es.data());
                    auto e2 = static_cast<std::size_t>(es[aresta].twin->prev - es.data());
                    lista_de_potencialmente_invalidas.push_back(e1);
                    lista_de_potencialmente_invalidas.push_back(e2);

                    dcel->deleta_aresta(aresta, false);
                    dcel->inclui_aresta(p1, p2);
                } else {
                    std::cout << " ok" << std::endl;
                }
            } else {
                std::cout << " nao pode ser flipada" << std::endl;
            }
        }

    }

    void inicia() {
        double lim_inf = (((0) / x) * 2.0) - 1.0 + (1.0 / x);
        double lim_sup = ((static_cast<double>(x-1) / x) * 2.0) - 1.0 + (1.0 / x);

        pontos.insert({
            {lim_inf, lim_inf},
            {lim_sup, lim_inf},
            {lim_sup, lim_sup},
            {lim_inf, lim_sup}
        });
        dcel = std::make_unique<DCEL>(std::vector<Ponto>{
            {lim_inf, lim_inf},
            {lim_sup, lim_inf},
            {lim_sup, lim_sup},
            {lim_inf, lim_sup},
            {lim_inf, lim_inf}
        });
        dcel->inclui_aresta(0, 2);

    }

public:
    int x;
    int y;
    int n;
    unsigned char* image_data;

    unsigned vbo;
    unsigned vao;
    unsigned ebo;
    unsigned faces_vao;
    unsigned faces_ebo;
    unsigned extra_vbo;
    unsigned extra_vao;

    bool mostrando_linhas;
    std::size_t last_size;
    std::size_t edge_count;
    std::size_t triangle_count;
    std::size_t last_gen;
    EstadoDelaunay estado;
    EntradaDelaunay estado_entrada;
    std::set<Ponto> pontos;
    std::unique_ptr<DCEL> dcel;

    static const std::size_t max_floats = 512*1024;
};

const Cor base_delaunay {"#2b2831"};
const Cor cor_dly {"#4d9184"};
const Cor base_trabalho {"#30272b"};
const Cor cor_trabalho {"#87914d"};

int main() {
    GLFW::Session session {};
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
    estado.tela = Tela::DELAUNAY;

    glfwSetWindowUserPointer(window, &estado);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetKeyCallback(window, key_callback);

    glEnable(GL_MULTISAMPLE);

    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(65535u);

    Shader fixed_color_program {"shaders/color_line_vertex.glsl", "shaders/fixed_color_fragment.glsl"};
    Shader fixed_point_program {"shaders/color_line_vertex.glsl", "shaders/fixed_color_fragment.glsl"};
    Shader point_program {"shaders/point_vertex.glsl", "shaders/point_fragment.glsl"};
    Shader color_line_program {"shaders/color_line_vertex.glsl", "shaders/color_line_fragment.glsl"};

    Shader circle_program {"shaders/color_line_vertex.glsl", "shaders/circle_fragment.glsl", "shaders/halfcircles_geometry.glsl"};
    Shader quad_program {"shaders/color_line_vertex.glsl", "shaders/quad_fragment.glsl", "shaders/quad_geometry.glsl"};

    color_line_program.setFloat("alpha", 1.0f);
    point_program.setFloat("alpha", 1.0f);
    circle_program.setFloat("alpha", 0.8f);
    quad_program.setFloat("alpha", 1.0f);
    estado.pointSize = 4.0f;
    glLineWidth(estado.pointSize / 2.0f);

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
    glGenVertexArrays(1, &coisas_dcel.extra_vao);
    glBindVertexArray(coisas_dcel.extra_vao);

    glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof (float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5 * sizeof (float), reinterpret_cast<void*>(2 * sizeof (float)));
    glEnableVertexAttribArray(1);

    glGenBuffers(1, &coisas_dcel.extra_ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.extra_ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 16*1024*sizeof (unsigned), nullptr, GL_DYNAMIC_DRAW);

    glBindVertexArray(0);
    std::size_t mano = 0;
    std::size_t atividade_size = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(-1000.0, 1000.0);
    std::uniform_real_distribution<> dis_y(-1000.0, 1000.0);

    std::uniform_real_distribution<> ok_x(-1.0, 1.0);
    std::uniform_real_distribution<> ok_y(-1.0, 1.0);

    CoisasDelaunay delaunay {"teste0.bmp"};

    while (!glfwWindowShouldClose(window)) {

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        if (estado.tela == Tela::DCEL_TESTE) {

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
                coisas_dcel.last_gen = 0;
                coisas_dcel.dcel_ptr.reset();
                coisas_dcel.coisas_piscar = {};
                coisas_dcel.coisas_aresta = {};
                coisas_dcel.coisas_vertice = {};
                estad.estado = Dcel_Data::RECEBENDO;
            }

            if (estad.ponto_adicionado) {
                std::size_t diff = estad.poly.size() - coisas_dcel.last_size;

                std::vector<float> ps {};
                ps.reserve(diff * 5 * sizeof (float));

                std::vector<unsigned> is {};
                is.reserve(diff * 2 * sizeof (unsigned));
                for (std::size_t i = coisas_dcel.last_size; i < estad.poly.size(); ++i) {
                    auto ponto = estad.poly[i];
                    ps.push_back(ponto.x);
                    ps.push_back(ponto.y);

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
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLintptr>(coisas_dcel.edge_count * 2 * sizeof (unsigned)), static_cast<GLintptr>(is.size() * sizeof (unsigned)), is.data());
                std::size_t edge_diff = diff;
                if (coisas_dcel.last_size == 0) {
                    --edge_diff;
                }
                coisas_dcel.edge_count += edge_diff;

                coisas_dcel.last_size = estad.poly.size();
                estad.ponto_adicionado = false;
            }

            if (estad.poligono_fechado) {

                std::array<unsigned, 2> ultima_aresta {static_cast<unsigned>(estad.poly.size() - 1), 0};
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLintptr>(coisas_dcel.edge_count * 2 * sizeof (unsigned)), static_cast<GLintptr>(2 * sizeof (unsigned)), ultima_aresta.data());
                ++coisas_dcel.edge_count;

                estad.poligono_fechado = false;
            }

            if (estad.estado == Dcel_Data::CRIANDO_DCEL) {
                coisas_dcel.dcel_ptr = std::make_unique<DCEL>(estad.poly);
                estad.estado = Dcel_Data::DCEL_PRONTA;
            }

            if (coisas_dcel.dcel_ptr && coisas_dcel.last_gen < coisas_dcel.dcel_ptr->gen()) {

                auto [verts_r, v_invs] = coisas_dcel.dcel_ptr->vec_vertices();
                auto [edges_r, e_invs] = coisas_dcel.dcel_ptr->vec_edges();

                auto& verts = verts_r.get();
                auto& edges = edges_r.get();
                std::cout << verts.size() << std::endl;
                std::vector<float> ps {};
                ps.reserve(verts.size() * 5 * sizeof (float));
                for (std::size_t i = 0; i < verts.size(); ++i) {
                    auto ponto = verts[i];
                    if (v_invs.count(i)) {
                        ps.push_back(2.0f);
                        ps.push_back(2.0f);
                    } else {
                        ps.push_back(ponto.xy.x);
                        ps.push_back(ponto.xy.y);
                    }

                    ps.push_back(0.788f);
                    ps.push_back(0.682f);
                    ps.push_back(0.078f);
                }

                std::vector<unsigned> is {};
                is.reserve((edges.size() / 2) * sizeof (unsigned));
                for (std::size_t i = 0; i < (edges.size() / 2); ++i) {
                    if (e_invs.count(2*i)) {
                        is.push_back(65535);
                        is.push_back(65535);
                    } else {
                        unsigned p1 = static_cast<unsigned>(edges[2*i].origin - &verts[0]);
                        unsigned p2 = static_cast<unsigned>(edges[2*i + 1].origin - &verts[0]);
                        is.push_back(p1);
                        is.push_back(p2);
                    }
                }
                glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(ps.size() * sizeof (float)), ps.data());
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, static_cast<GLintptr>(is.size() * sizeof (unsigned)), is.data());
                coisas_dcel.edge_count = edges.size() / 2;
                coisas_dcel.last_size = verts.size();

                coisas_dcel.last_gen = coisas_dcel.dcel_ptr->gen();
            }

            if (estad.estado == Dcel_Data::DCEL_PRONTA) {
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
                Ponto mouse = ponto_xy();
                double menor_d_vertice = std::numeric_limits<double>::infinity();
                std::size_t menor_i_vertice = 0;
                auto [vs_r, v_iv] = coisas_dcel.dcel_ptr->vec_vertices();
                auto& vs = vs_r.get();
                for (std::size_t i = 0; i < vs.size(); ++i) {
                    if (v_iv.count(i)) {
                        continue;
                    }
                    double d = dist(vs[i].xy, mouse);
                    if (d < menor_d_vertice) {
                        menor_d_vertice = d;
                        menor_i_vertice = i;
                    }
                }

                double menor_d_aresta = std::numeric_limits<double>::infinity();
                std::size_t menor_i_aresta = 0;
                auto [es_r, e_iv] = coisas_dcel.dcel_ptr->vec_edges();
                auto& es = es_r.get();
                for (std::size_t i = 0; i < es.size(); i += 2) {
                    if (e_iv.count(i)) {
                        continue;
                    }
                    double d = distancia_ponto_segmento(es[i].origin->xy, es[i+1].origin->xy, mouse);
                    if (d < menor_d_aresta) {
                        menor_d_aresta = d;
                        menor_i_aresta = i;
                    }
                }
                auto& p1 = es[menor_i_aresta].origin->xy;
                auto& p2 = es[menor_i_aresta+1].origin->xy;

                double menor_d = menor_d_vertice;
                std::size_t menor_i = menor_i_vertice;
                bool vertice = true;
                if (menor_d_aresta < menor_d_vertice) {
                    menor_d = menor_d_aresta;
                    menor_i = menor_i_aresta;
                    vertice = false;
                }
                if (menor_d <= 0.05) {
                    if (vertice) {
                        glBindVertexArray(coisas_dcel.vao);
                        glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

                        point_program.use();
                        point_program.setFloat("pointRadius", estado.pointSize + 30.0f);
                        point_program.setFloat("alpha", 0.5f);
                        glDrawArrays(GL_POINTS, menor_i_vertice, 1);
                        point_program.setFloat("alpha", 1.0f);
                    } else {
                        glBindVertexArray(coisas_dcel.vao);
                        glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

                        point_program.use();
                        point_program.setFloat("pointRadius", estado.pointSize + 30.0f);
                        point_program.setFloat("alpha", 0.5f);
                        glDrawElements(GL_POINTS, 2, GL_UNSIGNED_INT, reinterpret_cast<void*>(menor_i_aresta * sizeof (unsigned)));
                        point_program.setFloat("alpha", 1.0f);
                    }
                }

                while (estad.operacoes.size() > 0) {
                    auto op = estad.operacoes.front();
                    estad.operacoes.pop_front();
                    if (op.op == Dcel_Op::ENCONTRAR_FACE) {

                        std::cout << "clicado no ponto " << mouse.x << ' ' << mouse.y << std::endl;
                        auto a = coisas_dcel.dcel_ptr->qual_face(op.args.ponto);
                        std::cout << "Indice da face: " << a << std::endl;
                    } else if (op.op == Dcel_Op::PISCAR_FACE) {
                        auto a = coisas_dcel.dcel_ptr->qual_face(op.args.ponto);
                        auto b = coisas_dcel.dcel_ptr->indices_das_arestas_de_uma_face(a);
                        coisas_dcel.coisas_piscar.arestas = b;
                        coisas_dcel.coisas_piscar.ticks = 0;
                        coisas_dcel.coisas_piscar.ticks_por_aresta = 50;
                        coisas_dcel.coisas_piscar.atual = 0;
                        estad.estado = Dcel_Data::PISCANDO;
                        break;
                    } else if (op.op == Dcel_Op::CLIQUE_VERTICE) {
                        if (menor_d > 0.05) {
                            continue;
                        }
                        if (vertice) {
                            std::cout << "Indice do vertice: " << menor_i << std::endl;
                            std::cout << "Indice do vertice valido: " << coisas_dcel.dcel_ptr->indice_vertice_valido() << std::endl;
                        } else {
                            std::cout << "Indice da aresta: " << menor_i << std::endl;
                            std::cout << "Indice da aresta gemea: " << menor_i+1 << std::endl;
                            std::cout << "Indice da proxima aresta: " << es[menor_i].next - es.data() << std::endl;
                            std::cout << "Indice da aresta anterior: " << es[menor_i].prev - es.data() << std::endl;
                            std::cout << "Indice da proxima aresta da gemea: " << es[menor_i+1].next - es.data() << std::endl;
                            std::cout << "Indice da aresta anterior da gemea: " << es[menor_i+1].prev - es.data() << std::endl;
                        }
                    } else {
                        estad.operacoes.push_front(op);
                        break;
                    }
                }
            }

            if (estad.estado == Dcel_Data::ADICIONANDO_ARESTA) {
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
                Ponto mouse = ponto_xy();
                double menor_d = std::numeric_limits<double>::infinity();
                std::size_t menor_i = 0;
                auto [vs_r, iv] = coisas_dcel.dcel_ptr->vec_vertices();
                auto& vs = vs_r.get();
                for (std::size_t i = 0; i < vs.size(); ++i) {
                    if (iv.count(i)) {
                        continue;
                    }
                    double d = dist(vs[i].xy, mouse);
                    if (d < menor_d) {
                        menor_d = d;
                        menor_i = i;
                    }
                }
                if (menor_d <= 0.05) {

                    glBindVertexArray(coisas_dcel.vao);
                    glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

                    point_program.use();
                    point_program.setFloat("pointRadius", estado.pointSize + 30.0f);
                    point_program.setFloat("alpha", 0.5f);
                    glDrawArrays(GL_POINTS, menor_i, 1);
                    point_program.setFloat("alpha", 1.0f);
                }

                while (estad.operacoes.size() > 0) {
                    auto op = estad.operacoes.front();
                    estad.operacoes.pop_front();
                    if (op.op == Dcel_Op::PONTO_SELECIONADO) {
                        if (menor_d > 0.05) {
                            coisas_dcel.coisas_aresta.indices_obtidos = 0;
                            continue;
                        }
                        std::cout << "bla" << std::endl;
                        if (coisas_dcel.coisas_aresta.indices_obtidos == 0) {
                            coisas_dcel.coisas_aresta.p1_idx = menor_i;
                            ++coisas_dcel.coisas_aresta.indices_obtidos;
                        } else if (coisas_dcel.coisas_aresta.p1_idx == menor_i) {

                            continue;
                        } else {
                            coisas_dcel.coisas_aresta.p2_idx = menor_i;
                            auto& v1_idx = coisas_dcel.coisas_aresta.p1_idx;
                            auto& v2_idx = coisas_dcel.coisas_aresta.p2_idx;

                            coisas_dcel.dcel_ptr->inclui_aresta(v1_idx, v2_idx);
                            estad.estado = Dcel_Data::DCEL_PRONTA;
                            coisas_dcel.coisas_aresta = {};
                            if (coisas_dcel.dcel_ptr->gen() > coisas_dcel.last_gen) {

                                continue;
                            }
                        }
                    } else {
                        estad.operacoes.push_front(op);
                        break;
                    }
                }
            }

            if (estad.estado == Dcel_Data::ADICIONANDO_VERTICE) {
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
                Ponto mouse = ponto_xy();
                double menor_d = std::numeric_limits<double>::infinity();
                std::size_t menor_i = 0;
                auto [es_r, iv] = coisas_dcel.dcel_ptr->vec_edges();
                auto& es = es_r.get();
                for (std::size_t i = 0; i < es.size(); i += 2) {
                    if (iv.count(i)) {
                        continue;
                    }
                    double d = distancia_ponto_segmento(es[i].origin->xy, es[i+1].origin->xy, mouse);
                    if (d < menor_d) {
                        menor_d = d;
                        menor_i = i;
                    }
                }
                auto& p1 = es[menor_i].origin->xy;
                auto& p2 = es[menor_i+1].origin->xy;
                if (menor_d <= 0.05) {

                    glBindVertexArray(coisas_dcel.vao);
                    glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

                    point_program.use();
                    point_program.setFloat("pointRadius", estado.pointSize + 30.0f);
                    point_program.setFloat("alpha", 0.5f);
                    glDrawElements(GL_POINTS, 2, GL_UNSIGNED_INT, reinterpret_cast<void*>(menor_i * sizeof (unsigned)));
                    point_program.setFloat("alpha", 1.0f);
                }

                while (estad.operacoes.size() > 0) {
                    auto op = estad.operacoes.front();
                    estad.operacoes.pop_front();
                    if (op.op == Dcel_Op::CLIQUE_VERTICE) {
                        if (menor_d > 0.05) {
                            coisas_dcel.coisas_vertice.aresta_selecionada = false;
                            continue;
                        }
                        std::cout << "bla2" << std::endl;
                        if (!coisas_dcel.coisas_vertice.aresta_selecionada) {
                            coisas_dcel.coisas_vertice.aresta_idx = menor_i;
                            coisas_dcel.coisas_vertice.aresta_selecionada = true;
                        } else if (coisas_dcel.coisas_vertice.aresta_idx != menor_i) {

                            coisas_dcel.coisas_vertice.aresta_idx = menor_i;
                            continue;
                        } else {
                            double s = sombra_reta_ponto(mouse, {p1, p2});

                            coisas_dcel.dcel_ptr->inclui_vertice_em_aresta(menor_i, s);
                            estad.estado = Dcel_Data::DCEL_PRONTA;
                            coisas_dcel.coisas_vertice = {};
                            if (coisas_dcel.dcel_ptr->gen() > coisas_dcel.last_gen) {

                                continue;
                            }
                        }
                    } else {
                        estad.operacoes.push_front(op);
                        break;
                    }
                }
            }

            if (estad.estado == Dcel_Data::ESPERANDO_ORBITA) {
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
                Ponto mouse = ponto_xy();
                double menor_d = std::numeric_limits<double>::infinity();
                std::size_t menor_i = 0;
                auto [vs_r, iv] = coisas_dcel.dcel_ptr->vec_vertices();
                auto& vs = vs_r.get();
                for (std::size_t i = 0; i < vs.size(); ++i) {
                    if (iv.count(i)) {
                        continue;
                    }
                    double d = dist(vs[i].xy, mouse);
                    if (d < menor_d) {
                        menor_d = d;
                        menor_i = i;
                    }
                }
                if (menor_d <= 0.05) {

                    glBindVertexArray(coisas_dcel.vao);
                    glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

                    point_program.use();
                    point_program.setFloat("pointRadius", estado.pointSize + 30.0f);
                    point_program.setFloat("alpha", 0.5f);
                    glDrawArrays(GL_POINTS, menor_i, 1);
                    point_program.setFloat("alpha", 1.0f);
                }

                while (estad.operacoes.size() > 0) {
                    auto op = estad.operacoes.front();
                    estad.operacoes.pop_front();
                    if (op.op == Dcel_Op::PISCAR_ORBITA) {
                        if (menor_d > 0.05) {
                            continue;
                        }
                        std::cout << "bla3" << std::endl;
                        auto b = coisas_dcel.dcel_ptr->indices_orbita_de_um_vertice(menor_i);
                        coisas_dcel.coisas_piscar.arestas = b;
                        coisas_dcel.coisas_piscar.ticks = 0;
                        coisas_dcel.coisas_piscar.ticks_por_aresta = 50;
                        coisas_dcel.coisas_piscar.atual = 0;
                        estad.estado = Dcel_Data::PISCANDO;
                    } else {
                        estad.operacoes.push_front(op);
                        break;
                    }
                }
            }

            if (estad.estado == Dcel_Data::DELETANDO_ARESTA) {
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
                Ponto mouse = ponto_xy();
                double menor_d = std::numeric_limits<double>::infinity();
                std::size_t menor_i = 0;
                auto [es_r, iv] = coisas_dcel.dcel_ptr->vec_edges();
                auto& es = es_r.get();
                for (std::size_t i = 0; i < es.size(); i += 2) {
                    if (iv.count(i)) {
                        continue;
                    }
                    double d = distancia_ponto_segmento(es[i].origin->xy, es[i+1].origin->xy, mouse);
                    if (d < menor_d) {
                        menor_d = d;
                        menor_i = i;
                    }
                }
                auto& p1 = es[menor_i].origin->xy;
                auto& p2 = es[menor_i+1].origin->xy;
                if (menor_d <= 0.05) {

                    glBindVertexArray(coisas_dcel.vao);
                    glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

                    point_program.use();
                    point_program.setFloat("pointRadius", estado.pointSize + 30.0f);
                    point_program.setFloat("alpha", 0.5f);
                    glDrawElements(GL_POINTS, 2, GL_UNSIGNED_INT, reinterpret_cast<void*>(menor_i * sizeof (unsigned)));
                    point_program.setFloat("alpha", 1.0f);
                }

                while (estad.operacoes.size() > 0) {
                    auto op = estad.operacoes.front();
                    estad.operacoes.pop_front();
                    if (op.op == Dcel_Op::CLIQUE_VERTICE) {
                        if (menor_d > 0.05) {
                            continue;
                        }
                        coisas_dcel.dcel_ptr->deleta_aresta(menor_i);
                        estad.estado = Dcel_Data::DCEL_PRONTA;
                        if (coisas_dcel.dcel_ptr->vazia()) {
                            estad.estado = Dcel_Data::RESETANDO;
                            continue;
                        }
                    } else {
                        estad.operacoes.push_front(op);
                        break;
                    }
                }
            }

            if (estad.estado == Dcel_Data::DELETANDO_VERTICE) {
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
                Ponto mouse = ponto_xy();
                double menor_d = std::numeric_limits<double>::infinity();
                std::size_t menor_i = 0;
                auto [vs_r, iv] = coisas_dcel.dcel_ptr->vec_vertices();
                auto& vs = vs_r.get();
                for (std::size_t i = 0; i < vs.size(); ++i) {
                    if (iv.count(i)) {
                        continue;
                    }
                    double d = dist(vs[i].xy, mouse);
                    if (d < menor_d) {
                        menor_d = d;
                        menor_i = i;
                    }
                }
                if (menor_d <= 0.05) {

                    glBindVertexArray(coisas_dcel.vao);
                    glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

                    point_program.use();
                    point_program.setFloat("pointRadius", estado.pointSize + 30.0f);
                    point_program.setFloat("alpha", 0.5f);
                    glDrawArrays(GL_POINTS, menor_i, 1);
                    point_program.setFloat("alpha", 1.0f);
                }

                while (estad.operacoes.size() > 0) {
                    auto op = estad.operacoes.front();
                    estad.operacoes.pop_front();
                    if (op.op == Dcel_Op::CLIQUE_VERTICE) {
                        if (menor_d > 0.05) {
                            continue;
                        }
                        coisas_dcel.dcel_ptr->deleta_vertice(menor_i);
                        estad.estado = Dcel_Data::DCEL_PRONTA;
                        if (coisas_dcel.dcel_ptr->vazia()) {
                            estad.estado = Dcel_Data::RESETANDO;
                            continue;
                        }
                    } else {
                        estad.operacoes.push_front(op);
                        break;
                    }
                }
            }

            glBindVertexArray(coisas_dcel.vao);
            glBindBuffer(GL_ARRAY_BUFFER, coisas_dcel.vbo);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, coisas_dcel.ebo);

            color_line_program.use();
            color_line_program.setFloat("alpha", 1.0f);
            glDrawElements(GL_LINES, coisas_dcel.edge_count*2, GL_UNSIGNED_INT, nullptr);

            if (estad.estado == Dcel_Data::PISCANDO) {
                auto& c = coisas_dcel.coisas_piscar;
                if (c.ticks == 0) {
                    auto [edges_r, e_invs] = coisas_dcel.dcel_ptr->vec_edges();
                    auto& edges = edges_r.get();

                    auto& e = edges[c.arestas[c.atual]];
                    std::cout << c.arestas[c.atual] << std::endl;
                }
                fixed_color_program.use();
                fixed_color_program.setFloat("alpha", 0.8f);
                glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, reinterpret_cast<void*>((c.arestas[c.atual] & (~1llu)) * sizeof (unsigned)));

                ++c.ticks;
                if (c.ticks >= c.ticks_por_aresta) {
                    c.ticks = 0;
                    ++c.atual;
                    if (c.atual >= c.arestas.size()) {
                        estad.estado = Dcel_Data::DCEL_PRONTA;
                    }
                }
            }

            point_program.use();
            point_program.setFloat("pointRadius", estado.pointSize);
            glDrawArrays(GL_POINTS, 0, coisas_dcel.last_size);
        } else if (estado.tela == Tela::DELAUNAY) {
            glClearColor(base_trabalho.r(), base_trabalho.g(), base_trabalho.b(), 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            std::size_t novos_pontos_aleatorios = 0;
            bool outra_tela = false;

            while (estado.estado_delaunay.eventos.size() > 0) {
                auto op = estado.estado_delaunay.eventos.front();
                estado.estado_delaunay.eventos.pop_front();
                if (delaunay.estado == EstadoDelaunay::OK) {
                    switch (op.op) {
                        case General_Op::TECLA:
                            if (op.button_key == GLFW_KEY_T && !op.mods) {
                                delaunay.mostrando_linhas = !delaunay.mostrando_linhas;
                            } else if (op.button_key == GLFW_KEY_R) {

                                if (!op.mods) {
                                    ++novos_pontos_aleatorios;
                                } else if (op.mods == GLFW_MOD_SHIFT) {
                                    novos_pontos_aleatorios += 10;
                                } else if (op.mods == GLFW_MOD_CONTROL) {
                                    novos_pontos_aleatorios += 25;
                                } else if (op.mods == (GLFW_MOD_SHIFT | GLFW_MOD_CONTROL)) {
                                    novos_pontos_aleatorios += 500;
                                }
                            } else if (op.button_key == GLFW_KEY_E && !op.mods) {
                                estado.tela = Tela::DCEL_TESTE;

                                auto& estad = estado.estado_dcel_teste;
                                estad.poly = {Ponto{0.0, 0.0}};
                                estad.operacoes.clear();
                                estad.ponto_adicionado = false;
                                estad.poligono_fechado = false;

                                coisas_dcel.vao = delaunay.vao;
                                coisas_dcel.vbo = delaunay.vbo;
                                coisas_dcel.ebo = delaunay.ebo;
                                coisas_dcel.last_size = delaunay.last_size;
                                coisas_dcel.edge_count = delaunay.edge_count;
                                coisas_dcel.last_gen = delaunay.last_gen - 1;
                                coisas_dcel.dcel_ptr = std::move(delaunay.dcel);
                                coisas_dcel.coisas_piscar = {};
                                coisas_dcel.coisas_aresta = {};
                                coisas_dcel.coisas_vertice = {};
                                estad.estado = Dcel_Data::DCEL_PRONTA;
                                outra_tela = true;
                            }
                            break;
                        default:
                            break;
                    }
                }
            }
            if (outra_tela) continue;

            if (novos_pontos_aleatorios > 0) {
                while (novos_pontos_aleatorios --> 0 && !aconteceu_aquilo) {
                    Ponto p {ok_x(gen), ok_y(gen)};
                    bool foi = delaunay.adiciona_ponto(p);
                    if (aconteceu_aquilo) {
                        std::cout << "alerta" << std::endl;
                        std::cout << "alerta" << std::endl;
                        std::cout << "alerta" << std::endl;
                        std::cout << "alerta" << std::endl;
                        std::cout << "alerta" << std::endl;
                        std::cout << "investigar problema que aconteceu" << std::endl;
                        break;
                    }
                    while (!foi) {
                        p = {ok_x(gen), ok_y(gen)};
                        foi = delaunay.adiciona_ponto(p);
                        if (aconteceu_aquilo) {
                            std::cout << "alerta" << std::endl;
                            std::cout << "alerta" << std::endl;
                            std::cout << "alerta" << std::endl;
                            std::cout << "alerta" << std::endl;
                            std::cout << "alerta" << std::endl;
                            std::cout << "investigar problema que aconteceu" << std::endl;
                            break;
                        }
                    }
                }
            }

            if (delaunay.last_gen < delaunay.dcel->gen()) {

                auto [verts_r, v_invs] = delaunay.dcel->vec_vertices();
                auto [edges_r, e_invs] = delaunay.dcel->vec_edges();
                auto [faces_r, f_invs] = delaunay.dcel->vec_faces();

                auto& verts = verts_r.get();
                auto& edges = edges_r.get();
                auto& faces = faces_r.get();

                std::vector<float> ps {};
                ps.reserve(verts.size() * 5 * sizeof (float));
                for (std::size_t i = 0; i < verts.size(); ++i) {
                    auto ponto = verts[i];
                    auto cor_ponto = delaunay.encontra_cor(ponto.xy);
                    if (v_invs.count(i)) {
                        ps.push_back(2.0f);
                        ps.push_back(2.0f);
                    } else {
                        ps.push_back(ponto.xy.x);
                        ps.push_back(ponto.xy.y);
                    }
                    ps.push_back(cor_ponto.r());
                    ps.push_back(cor_ponto.g());
                    ps.push_back(cor_ponto.b());
                }

                std::vector<unsigned> is {};
                is.reserve((edges.size() / 2) * sizeof (unsigned));
                for (std::size_t i = 0; i < (edges.size() / 2); ++i) {
                    if (e_invs.count(2*i)) {
                        is.push_back(65535);
                        is.push_back(65535);
                    } else {
                        unsigned p1 = static_cast<unsigned>(edges[2*i].origin - &verts[0]);
                        unsigned p2 = static_cast<unsigned>(edges[2*i + 1].origin - &verts[0]);
                        is.push_back(p1);
                        is.push_back(p2);
                    }
                }

                glBindVertexArray(delaunay.vao);
                glBindBuffer(GL_ARRAY_BUFFER, delaunay.vbo);
                glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLintptr>(ps.size() * sizeof (float)), ps.data());
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, delaunay.ebo);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, static_cast<GLintptr>(is.size() * sizeof (unsigned)), is.data());

                std::vector<unsigned> fs {};
                fs.reserve(faces.size() * 3 * sizeof (unsigned));
                for (std::size_t i = 1; i < faces.size(); ++i) {
                    if (f_invs.count(i)) {
                        fs.push_back(65535);
                        fs.push_back(65535);
                        fs.push_back(65535);
                    } else {
                        auto vs = delaunay.dcel->indices_dos_vertices_de_uma_face(i);
                        if (vs.size() != 3) {
                            std::cerr << "aviso: vai dar errado" << std::endl;
                        }
                        for (auto v : vs) {
                            unsigned p = static_cast<unsigned>(v);
                            fs.push_back(p);
                        }
                    }
                }
                glBindVertexArray(delaunay.faces_vao);
                glBindBuffer(GL_ARRAY_BUFFER, delaunay.vbo);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, delaunay.faces_ebo);
                glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, static_cast<GLintptr>(fs.size() * sizeof (unsigned)), fs.data());
                delaunay.edge_count = edges.size() / 2;
                delaunay.last_size = verts.size();
                delaunay.triangle_count = faces.size() - 1;

                delaunay.last_gen = delaunay.dcel->gen();
            }
            if (delaunay.estado == EstadoDelaunay::OK) {
                glBindVertexArray(delaunay.faces_vao);
                glBindBuffer(GL_ARRAY_BUFFER, delaunay.vbo);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, delaunay.faces_ebo);

                color_line_program.use();
                color_line_program.setFloat("alpha", 1.0f);
                glDrawElements(GL_TRIANGLES, delaunay.triangle_count*3, GL_UNSIGNED_INT, nullptr);

                if (delaunay.mostrando_linhas) {
                    glBindVertexArray(delaunay.vao);
                    glBindBuffer(GL_ARRAY_BUFFER, delaunay.vbo);
                    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, delaunay.ebo);

                    glLineWidth(4.0f);
                    fixed_color_program.use();
                    fixed_color_program.setFloat("alpha", 1.0f);
                    glDrawElements(GL_LINES, delaunay.edge_count*2, GL_UNSIGNED_INT, nullptr);
                    glLineWidth(std::max(estado.pointSize / 2.0f, 1.0f));

                    fixed_point_program.use();
                    fixed_point_program.setFloat("pointRadius", 8.0f);

                    glDrawArrays(GL_POINTS, 0, delaunay.last_size);
                }
            }

        }

        glfwSwapBuffers(window);
        session.pollEvents();
    }

    return 0;
}