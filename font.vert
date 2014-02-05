#version 150 core
in vec4 coord;
out vec2 frag_tex_coord;

void main() {
    frag_tex_coord = coord.zw;
    gl_Position = vec4(coord.xy, 0, 1);
}
