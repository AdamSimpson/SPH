#version 150 core
in vec2 position;
in vec2 tex_coord;

uniform float radius;
uniform vec2 center;
uniform vec3 color;

uniform mat4 view;
uniform mat4 proj;

out vec2 frag_tex_coord;

void main() {
    frag_tex_coord = tex_coord;
    gl_Position = proj*view*vec4(position, 0.0, 1.0);
}

