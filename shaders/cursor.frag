#version 150 core
in vec2 frag_tex_coord;

uniform sampler2D tex;

out vec4 OutColor;

void main() {
    OutColor = texture(tex, frag_tex_coord);
}
