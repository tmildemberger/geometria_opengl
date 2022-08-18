#version 450 core
layout (location = 0) out vec4 fragColor;

// in vec3 ourColor;
in vec2 TexCoord;

uniform sampler2D texture0;
uniform sampler2D texture1;

// uniform float mixing;
// uniform float otherColorRed;

void main() {
    // fragColor = vec4(ourColor, 1.0f);
    vec2 faceCoord = 1 * TexCoord;
    // fragColor = mix(texture(texture0, TexCoord) * vec4(ourColor, 1.0),
    //                 texture(texture1, TexCoord), 0.2);
    // fragColor = mix(texture(texture0, TexCoord) * vec4(otherColorRed, 0.2, 0.2, 1.0),
    //                 texture(texture1, faceCoord), mixing);
    fragColor = mix(texture(texture0, TexCoord),
                    texture(texture1, faceCoord), 0.2);
    // fragColor = vec4(0.2f, 0.4f, 0.7f, 1.0f);
}