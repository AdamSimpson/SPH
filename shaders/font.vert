#version 150 core
in vec4 coord;
in vec3 color;

out vec2 frag_tex_coord;
out vec3 text_color;

void main() {
    frag_tex_coord = coord.zw;
    text_color = color;
    gl_Position = vec4(coord.xy, 0, 1);
}
