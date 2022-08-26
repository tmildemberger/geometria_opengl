#version 330 core
layout (location = 0) out vec4 fragColor;

uniform float alpha;

void main() {
    fragColor = vec4(0.2f, 0.4f, 0.7f, alpha);
}