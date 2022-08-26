#include "stb_img.hpp"

namespace STB_IMG {
#define STBI_ONLY_PNG
// #define STBI_NO_STDIO
#define STB_IMAGE_IMPLEMENTATION
#pragma GCC diagnostic push
//needed to gcc
#pragma GCC diagnostic ignored "-Wstrict-overflow"
#pragma GCC diagnostic ignored "-Wswitch-default"
//needed to clang
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wunused-function"
#include "stb_image.h"
#pragma GCC diagnostic pop
// #include <iostream>
#include <exception>
#include <stdexcept>

    Image::Image(::std::string filename) {
        stbi_set_flip_vertically_on_load(true);
        this->data = stbi_load(
            filename.c_str(), &this->width, &this->height, &this->nChannels, 0);
        if (!this->data) {
            throw ::std::runtime_error {"Error loading image file " + filename};
        }
    }
    Image::~Image() {
        // ::std::cout << "I was released\n";
        stbi_image_free(this->data);
    }
}