#version 150 core
in vec2 local_frag_coord;

out vec4 out_color;

uniform vec4 color;

void main() {
    out_color = color;
}

