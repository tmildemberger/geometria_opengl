#ifndef DRAWING_HPP
#define DRAWING_HPP

#include <array>
#include <vector>
#include <cstddef>
#include <memory>
#include <glad/glad.h>
#include "renderer.hpp"

struct Attribute {
    unsigned int index;
    int size;
    unsigned int type;
    int normalized;
    int offset;
};

class BufferLayout {
private:
public:
    unsigned int index { 0 };
    int stride { 0 };
    unsigned int number { 0 };
    std::vector<Attribute> attributes;
    BufferLayout() {};

    template <typename T> void push (int count);
};

// #ifndef DRAWING_IMPL
// #define DRAWING_IMPL

// #endif

struct Sharing {
    unsigned int count { 1 };
};

#include <iostream>

class Drawing {
private:
public:
    unsigned int VBO { 0 };
    unsigned int VAO { 0 };
    std::size_t vertices { 0 };
    std::shared_ptr<Sharing> sharing;

    template <typename T, std::size_t N>
    Drawing(std::array<T, N> attributes, BufferLayout layout) : 
        vertices { N / layout.number } {
        glCreateBuffers(1, &this->VBO);
        glNamedBufferStorage(this->VBO, N * sizeof (T),
                             attributes.data(), GL_DYNAMIC_STORAGE_BIT);
        
        glCreateVertexArrays(1, &this->VAO);
        glVertexArrayVertexBuffer(VAO, 0, VBO, 0, layout.stride);

        for (const auto& attr : layout.attributes) {
            glEnableVertexArrayAttrib(VAO, attr.index);
            glVertexArrayAttribFormat(VAO,
                attr.index, attr.size, attr.type, attr.normalized, static_cast<GLuint>(attr.offset));
            glVertexArrayAttribBinding(VAO, attr.index, 0);
        }
        sharing = std::make_shared<Sharing>();
        // std::cout << "Drawing created with new sharing\n";
        // std::cout << "number of vertices:" << this->vertices << "\n\n";
    }

    Drawing(const Drawing& that) : VBO { that.VBO },
                                   VAO { that.VAO }, 
                                   vertices { that.vertices }, 
                                   sharing { that.sharing } {
        ++sharing->count;
        // std::cout << "Drawing created with old sharing\n";
        // std::cout << "sharing count:" << sharing->count << "\n\n";
    }

    Drawing& operator=(const Drawing& that) = delete;
    const Drawing& operator=(const Drawing& that) const = delete;

    ~Drawing() {
        if (sharing->count == 1) {
            glDeleteBuffers(1, &this->VBO);
            glDeleteVertexArrays(1, &this->VAO);
            // std::cout << "they got deleted only one time\n";
        } else {
            --sharing->count;
        }
    }

    void render(const Renderer& renderer);
};

#endif // DRAWING_HPP