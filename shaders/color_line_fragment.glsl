#version 330 core
layout (location = 0) out vec4 fragColor;

in vec3 theColor;
uniform float alpha;

void main() {
    fragColor = vec4(theColor.rgb, alpha);
}