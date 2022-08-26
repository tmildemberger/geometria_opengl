#version 330 core
layout (location = 0) out vec4 fragColor;

in vec3 theColor;

void main() {
    fragColor = vec4(theColor.rgb, 1.0);
}