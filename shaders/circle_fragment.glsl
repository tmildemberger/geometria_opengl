#version 330 core
layout (location = 0) out vec4 fragColor;

in vec3 outColor;
in float invPointSize;
in float pointSize;
in float aaaRadius;
in float borderSize;
in vec2 centerPos;
// float triplo = invPointSize * invPointSize;
uniform float alpha;
uniform vec3 color;

void main() {
    // fragColor = vec4(outColor.rgb, alpha / 3.0f);
    vec3 realColor = (color == 0) ? outColor : color;
    vec2 pos = 2.0f * ((gl_FragCoord.xy - 0.5f) / 600.0f) - 1.0f;
    float squared = dot(pos - centerPos, pos - centerPos) / (aaaRadius*aaaRadius);
    // if (squared < 1) {
    //     fragColor = vec4(outColor.rgb, alpha / 1.7f);
    // }
    
    // vec2 circCoord = 2.0 * gl_PointCoord - 1.0;
    // 
    // float squared = dot(circCoord, circCoord);
    // vec2 pos = (gl_FragCoord.xy - 0.5f) / 600.0f;
    // float squared = dot(pos - centerPos, pos - centerPos);
    
    if (squared > 1.0) {
        // fragColor = vec4(outColor.rgb, alpha / 2.0f);
        discard;
    //} else if (squared < 0.5) {
      //  fragColor = vec4(0.5, 0.2, 0.3, 0.0);
    } else {
        // fragColor = mix(vec4(outColor.rgb, 0.0), vec4(outColor.rgb, 1.0), smoothstep(1.0, 0.15, squared));
        // fragColor = mix(vec4(outColor.rgb, 0.0), vec4(outColor.rgb, 1.0), smoothstep(1.0, 1.0 - 1.0 * invPointSize, squared));
        // fragColor = mix(vec4(outColor.rgb, 0.0), vec4(outColor.rgb, 1.0), smoothstep(1.0, clamp((1.0 - 1.0 * invPointSize)*(2.6237*(1.0 - triplo)), squared));
        // smoothstep estava errado kk
        // finalmente: fragColor = mix(vec4(outColor.rgb, 0.0), vec4(outColor.rgb, 1.0), 1.0 - smoothstep(1.0 - (invPointSize)*(1.0 + 1.6237*(1.0 - smoothstep(7, 30, pointSize))), 1.0, squared));
        // fragColor = mix(vec4(outColor.rgb, 0.0), vec4(outColor.rgb, alpha), 1.0 - smoothstep(1.0 - (invPointSize)*(1.0 + 1.41*(1.0 - smoothstep(7, 30, pointSize))), 1.0, squared));
        
        if (squared > (pointSize-2*borderSize)*(pointSize-2*borderSize) / (pointSize*pointSize) ) {
            fragColor = mix(vec4(realColor.rgb, 0.0), vec4(realColor.rgb, alpha), 1.0 - smoothstep(1.0, 0.0, squared));
        } else {
            fragColor = vec4(realColor.rgb, alpha / 4.0f);
        }
        
        // 1.0 = pointSize^2
        // x = (pointSize - borderSize)^2

        // ps - bs = ps * x
        // 
        // if (squared < ()*())
        // fragColor = mix(vec4(outColor.rgb, 0.0), vec4(outColor.rgb, alpha), )
    }
}