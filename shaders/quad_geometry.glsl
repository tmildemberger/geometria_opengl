#version 330 core
layout (lines_adjacency) in;
// layout (points, max_vertices = 1) out;
layout (triangle_strip, max_vertices = 4) out;

in vec3 theColor[];

out vec3 outColor;

float area(vec4 p1, vec4 p2, vec4 p3) {
    return (p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y);
}

void main() {
    int i = 0;
    int i1 = 0;
    int i2 = 0;
    int i3 = 0;
    int i4 = 0;
    for (; i < 4; ++i) {
        vec4 p1 = gl_in[i].gl_Position;
        vec4 p2 = gl_in[(i+1)%4].gl_Position;
        vec4 p3 = gl_in[(i+2)%4].gl_Position;
        vec4 p4 = gl_in[(i+3)%4].gl_Position;

        if ((area(p1, p3, p2) > 0.0f) != (area(p1, p3, p4) > 0.0f)) {
            i2 = i;
            i3 = (i+2)%4;
            if (area(p1, p3, p2) > 0.0f) {
                i4 = (i+3)%4;
                i1 = (i+1)%4;
            } else {
                i1 = (i+3)%4;
                i4 = (i+1)%4;
            }
            break;
        }
    }

    outColor = theColor[i1];
    gl_Position = gl_in[i1].gl_Position;
    EmitVertex();

    outColor = theColor[i2];
    gl_Position = gl_in[i2].gl_Position;
    EmitVertex();

    outColor = theColor[i3];
    gl_Position = gl_in[i3].gl_Position;
    EmitVertex();

    outColor = theColor[i4];
    gl_Position = gl_in[i4].gl_Position;
    EmitVertex();
    EndPrimitive();
}

