#include "shader.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <exception>
#include <stdexcept>

Shader::Shader(const std::string& pathvert, const std::string& pathfrag) {
    try {
        unsigned int vertexID { compile(pathvert, GL_VERTEX_SHADER) };
        unsigned int fragmentID { compile(pathfrag, GL_FRAGMENT_SHADER) };

        this->ID = link(vertexID, fragmentID);

        glDeleteShader(vertexID);
        glDeleteShader(fragmentID);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

Shader::Shader(const std::string& pathvert, const std::string& pathfrag, const std::string& pathgeo) {
    try {
        unsigned int vertexID { compile(pathvert, GL_VERTEX_SHADER) };
        unsigned int fragmentID { compile(pathfrag, GL_FRAGMENT_SHADER) };
        unsigned int geometryID { compile(pathgeo, GL_GEOMETRY_SHADER) };

        this->ID = link(vertexID, fragmentID, geometryID);

        glDeleteShader(vertexID);
        glDeleteShader(fragmentID);
        glDeleteShader(geometryID);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

unsigned int Shader::compile(const std::string& path, GLenum type) const {
    std::ifstream file;
    file.open(path);
    if (!file) {
        throw std::system_error { errno, std::system_category(),
                                  "Failed to open " + path };
    }

    std::stringstream file_contents;
    file_contents << file.rdbuf();
    file.close();

    std::string codeString { file_contents.str() };
    const char* const shaderCode { codeString.c_str() };

    unsigned int shaderID { glCreateShader(type) };
    glShaderSource(shaderID, 1, &shaderCode, NULL);
    glCompileShader(shaderID);
    int success { 0 };
    glGetShaderiv(shaderID, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(shaderID, 512, NULL, infoLog);
        std::string message { "ERROR::SHADER::" };
        message.append((type == GL_VERTEX_SHADER) ? ("VERTEX") : ((type == GL_FRAGMENT_SHADER) ? "FRAGMENT" : "GEOMETRY"));
        message.append("::COMPILATION_FAILED\n");
        message.append(infoLog);
        throw std::runtime_error { message };
    }
    return shaderID;
}

unsigned int Shader::link(unsigned int vertID, unsigned int fragID) const {
    unsigned int shaderProgramID { glCreateProgram() };
    glAttachShader(shaderProgramID, vertID);
    glAttachShader(shaderProgramID, fragID);
    glLinkProgram(shaderProgramID);
    int success { 0 };
    glGetProgramiv(shaderProgramID, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetProgramInfoLog(shaderProgramID, 512, NULL, infoLog);
        throw std::runtime_error { "ERROR::SHADER::PROGRAM::LINK_FAILED\n" +
                                 std::string {infoLog} };
    }
    return shaderProgramID;
}

unsigned int Shader::link(unsigned int vertID, unsigned int fragID, unsigned int geoID) const {
    unsigned int shaderProgramID { glCreateProgram() };
    glAttachShader(shaderProgramID, vertID);
    glAttachShader(shaderProgramID, fragID);
    glAttachShader(shaderProgramID, geoID);
    glLinkProgram(shaderProgramID);
    int success { 0 };
    glGetProgramiv(shaderProgramID, GL_LINK_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetProgramInfoLog(shaderProgramID, 512, NULL, infoLog);
        throw std::runtime_error { "ERROR::SHADER::PROGRAM::LINK_FAILED\n" +
                                 std::string {infoLog} };
    }
    return shaderProgramID;
}

void Shader::use() const {
    glUseProgram(this->ID);
}

// void Shader::setBool(const std::string& name, const bool value) const {
//     glProgramUniform1i(this->ID, glGetUniformLocation(this->ID, name.c_str()),
//         static_cast<int>(value));
// }

// void Shader::setInt(const std::string& name, const int value) const {
//     glProgramUniform1i(this->ID, glGetUniformLocation(this->ID, name.c_str()),
//         value);
// }

void Shader::setFloat(const std::string& name, const float value) const {
    glUseProgram(this->ID);
    glUniform1f(glGetUniformLocation(this->ID, name.c_str()), value);
}

void Shader::setVec3(const std::string& name, const float value1, const float value2, const float value3) const {
    glUseProgram(this->ID);
    glUniform3f(glGetUniformLocation(this->ID, name.c_str()), value1, value2, value3);
}
/*
void Shader::setTransform(float* transformMatrix) const {
    glProgramUniformMatrix4fv(this->ID, glGetUniformLocation(this->ID,
        "transform"), 1, GL_FALSE, transformMatrix);
}

void Shader::setModelTransform(float* transformMatrix) const {
    glProgramUniformMatrix4fv(this->ID, glGetUniformLocation(this->ID,
        "model"), 1, GL_FALSE, transformMatrix);
}

void Shader::setViewTransform(float* transformMatrix) const {
    glProgramUniformMatrix4fv(this->ID, glGetUniformLocation(this->ID,
        "view"), 1, GL_FALSE, transformMatrix);
}

void Shader::setProjectionTransform(float* transformMatrix) const {
    glProgramUniformMatrix4fv(this->ID, glGetUniformLocation(this->ID,
        "projection"), 1, GL_FALSE, transformMatrix);
}
*/