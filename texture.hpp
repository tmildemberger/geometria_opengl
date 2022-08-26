#ifndef TEXTURE_HPP
#define TEXTURE_HPP

#include <string>
#include <glad/glad.h>

class Texture {
private:
public:
    unsigned int textureID { 0 };
    Texture(std::string filename, GLenum format);
    void useAsUnit(unsigned int unit) const;
};

#endif // TEXTURE_HPP