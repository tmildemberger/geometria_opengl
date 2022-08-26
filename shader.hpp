#ifndef SHADER_HPP
#define SHADER_HPP

#include <string>
#include <exception>
#include <glad/glad.h>

#include "thing.hpp"

class Shader
{
private:
    unsigned int ID { 0 }; // the shader program ID

    //used only by the constructor
    unsigned int compile(const std::string& path, GLenum type) const;
    unsigned int link(unsigned int vertID, unsigned int fragID) const;
public:
    Shader(const std::string& pathvert, const std::string& pathfrag);

    void use() const;

    // void setBool(const std::string& name, const bool value) const;
    // void setInt(const std::string& name, const int value) const;
    void setFloat(const std::string& name, const float value) const;

    /*void setTransform(float* transformMatrix) const;
    void setModelTransform(float* transformMatrix) const;
    void setViewTransform(float* transformMatrix) const;
    void setProjectionTransform(float* transformMatrix) const;*/
    /*
    template <typename ... Ts>
    void setUniform(const std::string& name, const Ts ... objs) const {
        glUseProgram(this->ID);
        overload(glUniform1f,
                 glUniform1i,
                 glUniform1ui,
                 glUniform2f,
                 glUniform2i,
                 glUniform2ui,
                 glUniform3f,
                 glUniform3i,
                 glUniform3ui,
                 glUniform4f,
                 glUniform4i,
                 glUniform4ui) 
                    (glGetUniformLocation(this->ID, name.c_str()),
                     objs...);
    }*/
};

// not gonna do a new exception class because there will be no
// actions other than warn the user
// class ShaderBuildError : public std::exception {
// private:
//     const std::string error;
// public:
//     ShaderBuildError(const std::string& err) noexcept : error {err} {

//     }
//     const char* what() const noexcept { return this->error.c_str(); }
// };

#endif // SHADER_HPP