#include "drawing.hpp"

void Drawing::render(const Renderer& renderer) {
    renderer.render(this->VAO, this->vertices);
}

template <> void BufferLayout::push<float> (int count) {
    this->attributes.push_back( {index, count, GL_FLOAT, GL_FALSE, stride} );
    ++index;
    number += static_cast<unsigned>(count);
    stride += count * static_cast<int>(sizeof (float));
}

template <> void BufferLayout::push<unsigned int> (int count) {
    this->attributes.push_back( 
        {index, count, GL_UNSIGNED_INT, GL_FALSE, stride} 
    );
    ++index;
    number += static_cast<unsigned>(count);
    stride += count * static_cast<int>(sizeof (unsigned int));
}