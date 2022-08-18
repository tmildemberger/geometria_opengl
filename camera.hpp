#ifndef CAMERA_HPP
#define CAMERA_HPP

#pragma GCC diagnostic push
//needed to gcc
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wpedantic"
//needed to clang
#pragma GCC diagnostic ignored "-Wgnu-anonymous-struct"
#pragma GCC diagnostic ignored "-Wnested-anon-types"
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#pragma GCC diagnostic pop

class Camera {
private:

public:
    float fov { 45.0f };

    glm::vec3 cameraPos { glm::vec3 {0.0f, 0.0f, 0.0f} };
    glm::vec3 cameraFront { glm::vec3 {0.0f, 0.0f, -1.0f} };
    glm::vec3 cameraUp { glm::vec3 {0.0f, 1.0f, 0.0f} };

    Camera() {};
    ~Camera() {};
};

#endif // CAMERA_HPP