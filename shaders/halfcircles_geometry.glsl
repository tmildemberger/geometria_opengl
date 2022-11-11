#version 330 core
layout (triangles) in;
// layout (points, max_vertices = 1) out;
layout (triangle_strip, max_vertices = 4) out;

in vec3 theColor[];

out vec3 outColor;
out float invPointSize;
out float pointSize;
out float aaaRadius;
out float borderSize;
out vec2 centerPos;

void main() {
    // utilizando fórmula do circuncentro de
    // https://www.omnicalculator.com/math/circumcenter-of-a-triangle
    vec2 p1 = gl_in[0].gl_Position.xy;
    vec2 p2 = gl_in[1].gl_Position.xy;
    vec2 p3 = gl_in[2].gl_Position.xy;
    float t = dot(p1, p1) - dot(p2, p2);
    float u = dot(p1, p1) - dot(p3, p3);
    float J = (p1.x - p2.x)*(p1.y - p3.y) - (p1.x - p3.x)*(p1.y - p2.y);
    float x = (-(p1.y - p2.y)*u + (p1.y - p3.y)*t)/(2*J);
    float y = ((p1.x - p2.x)*u - (p1.x - p3.x)*t)/(2*J);

    // float xba, yba, xca, yca;
    // float balength, calength;
    // float denominator;
    // float xcirca, ycirca;

    // /* Use coordinates relative to point `a' of the triangle. */
    // xba = p2.x - p1.x;
    // yba = p2.y - p1.y;
    // xca = p3.x - p1.x;
    // yca = p3.y - p1.y;
    // /* Squares of lengths of the edges incident to `a'. */
    // balength = xba * xba + yba * yba;
    // calength = xca * xca + yca * yca;

    // /* Take your chances with floating-point roundoff. */
    // denominator = 0.5 / (xba * yca - yba * xca);

    // /* Calculate offset (from `a') of circumcenter. */
    // xcirca = (yca * balength - yba * calength) * denominator;
    // ycirca = (xba * calength - xca * balength) * denominator;

    // vec2 centro = vec2(xcirca, ycirca);

    vec2 centro = vec2(x, y);
    float raio = sqrt(dot(centro - p1, centro - p1)) * 1200.0f + 20.0f;
    float rraio = sqrt(dot(centro - p1, centro - p1)) + 20.0f / 1200.0f;
    // raio = 50.0f;

    // gl_Position = vec4(centro, 0.0f, 1.0f);
    // gl_PointSize = raio;

////
    invPointSize = sqrt(1.0 / raio);
    pointSize = raio;
    outColor = vec3(0.0f, 1.0f, 1.0f) - theColor[0];
    aaaRadius = rraio;
    borderSize = 20.0f;
    centerPos = centro;
    gl_Position = vec4(centro, 0.0f, 1.0f) + vec4(-rraio, rraio, 0.0f, 0.0f);

    EmitVertex();

////
    invPointSize = sqrt(1.0 / raio);
    pointSize = raio;
    outColor = vec3(0.0f, 1.0f, 1.0f) - theColor[0];
    aaaRadius = rraio;
    borderSize = 20.0f;
    centerPos = centro;
    gl_Position = vec4(centro, 0.0f, 1.0f) + vec4(rraio, rraio, 0.0f, 0.0f);

    EmitVertex();

////
    invPointSize = sqrt(1.0 / raio);
    pointSize = raio;
    outColor = vec3(0.0f, 1.0f, 1.0f) - theColor[0];
    aaaRadius = rraio;
    borderSize = 20.0f;
    centerPos = centro;
    gl_Position = vec4(centro, 0.0f, 1.0f) + vec4(-rraio, -rraio, 0.0f, 0.0f);

    EmitVertex();

////
    invPointSize = sqrt(1.0 / raio);
    pointSize = raio;
    outColor = vec3(0.0f, 1.0f, 1.0f) - theColor[0];
    aaaRadius = rraio;
    borderSize = 20.0f;
    centerPos = centro;
    gl_Position = vec4(centro, 0.0f, 1.0f) + vec4(rraio, -rraio, 0.0f, 0.0f);

    EmitVertex();

////
    invPointSize = sqrt(1.0 / raio);
    pointSize = raio;
    outColor = vec3(0.0f, 1.0f, 1.0f) - theColor[0];
    borderSize = 20.0f;
    centerPos = centro;
    gl_Position = vec4(centro, 0.0f, 1.0f) + vec4(rraio, rraio, 0.0f, 0.0f);

    EmitVertex();
    
    EndPrimitive();
}

