#include "texture.hpp"
#include "stb_img.hpp"

Texture::Texture(std::string filename, GLenum format) {
    STB_IMG::Image img {filename};/*
    glCreateTextures(GL_TEXTURE_2D, 1, &this->textureID);

    glTextureParameteri(this->textureID, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTextureParameteri(this->textureID, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTextureParameteri(this->textureID, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureParameteri(this->textureID, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTextureStorage2D(this->textureID, 1, GL_RGB8, img.width, img.height);
    glTextureSubImage2D(this->textureID, 0, 0, 0, img.width, img.height,
        format, GL_UNSIGNED_BYTE, img.data);

    glGenerateTextureMipmap(this->textureID);*/
}

void Texture::useAsUnit(unsigned int unit) const {
    /*glBindTextureUnit(unit, this->textureID);*/
}