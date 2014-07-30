#version 150 core
in vec2 frag_tex_coord;

uniform sampler2D tex;

out vec4 OutColor;

void main() {
    vec4 color = texture(tex, frag_tex_coord);
    float alpha = 0.5*step(0.015, color.a);
    float white = 0.0;
    white = step(-0.03, -color.a);
    OutColor = vec4(white, white, 1.0, alpha);
}
