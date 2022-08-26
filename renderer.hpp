#ifndef RENDERER_HPP
#define RENDERER_HPP

#include <glad/glad.h>

class Renderer {
private:
public:
    Renderer(){};
    ~Renderer(){};

    void render(unsigned int VAO, GLsizei count) const;
};

#endif // RENDERER_HPP