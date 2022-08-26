#ifndef STB_IMG_HPP
#define STB_IMG_HPP

#include <iostream>
#include <cmath>
#include <string>

namespace STB_IMG {
#define STBI_ONLY_PNG
#include "stb_image.h"
    class Image {
    private:
    public:
        int width { 0 };
        int height { 0 };
        int nChannels { 0 };
        unsigned char* data { nullptr };
        Image(::std::string filename);
        ~Image();
    };
}

#endif // STB_IMG_HPP