#version 150 core
in vec2 frag_tex_coord;
in vec3 text_color;

out vec4 outColor;

uniform sampler2D tex;

void main() {
    outColor = vec4(1, 1, 1, texture(tex, frag_tex_coord).r)*vec4(text_color, 1);
}

