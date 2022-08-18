#version 450 core
layout (location = 0) in vec3 aPos;
// layout (location = 1) in vec3 aColor;
layout (location = 1) in vec2 aTexCoord;

// out vec3 ourColor;
out vec2 TexCoord;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
// uniform float offset_x;

void main() {
    // gl_Position = vec4(aPos.x + offset_x - 0.3f, aPos.yz, 1.0);
    // ourColor = vec3(aPos.x + 0.2f + offset_x, aPos.y + 0.5f, 0.8f);
    gl_Position = projection * view * model * vec4(aPos, 1.0);
    // ourColor = aColor;
    TexCoord = aTexCoord;
}