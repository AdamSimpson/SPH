#version 150 core
in vec2 frag_tex_coord;

uniform sampler2D tex;

out vec4 OutColor;

void main() {
    vec4 color = texture(tex, frag_tex_coord);
    float alpha = step(0.01, color.a);
    OutColor = vec4(0.0, 0.0, 1.0, 0.5*alpha);
}
