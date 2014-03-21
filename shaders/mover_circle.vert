#version 150 core
in vec2 position;
in vec2 tex_coord;

uniform float radius;
uniform vec2 center;
uniform vec3 color;

out vec2 frag_tex_coord;

void main() {
    frag_tex_coord = tex_coord;
    gl_Position = vec4(position, 0.0, 1.0);
}

