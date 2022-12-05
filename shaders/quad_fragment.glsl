#version 330 core
layout (location = 0) out vec4 fragColor;

in vec3 outColor;
uniform float alpha;

void main() {
    fragColor = vec4(outColor.rgb, alpha);
}