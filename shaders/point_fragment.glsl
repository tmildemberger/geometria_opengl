#version 330 core
layout (location = 0) out vec4 fragColor;

in vec3 theColor;

void main() {
    vec2 circCoord = 2.0 * gl_PointCoord - 1.0;
    float squared = dot(circCoord, circCoord);
    if (squared > 1.0) {
        discard;
    //} else if (squared < 0.5) {
      //  fragColor = vec4(0.5, 0.2, 0.3, 0.0);
    } else {
        fragColor = mix(vec4(theColor.rgb, 0.0), vec4(theColor.rgb, 1.0), smoothstep(1.0, 0.85, squared));
    }
}