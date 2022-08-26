#include "renderer.hpp"

#include <glad/glad.h>

void Renderer::render(unsigned int VAO, GLsizei count) const {
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, count);
}