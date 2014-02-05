#version 150 core
in vec2 frag_tex_coord;
out vec4 outColor;

uniform vec4 color;
uniform sampler2D tex;

void main() {
    outColor = vec4(1, 1, 1, texture(tex, frag_tex_coord).r)*color;
}

