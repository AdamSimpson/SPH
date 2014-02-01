#version 150 core
in vec2 position;
in vec2 tex_coord;
in vec4 color;

uniform ivec2 screen_dims;

out vec2 frag_tex_coord;
out vec4 frag_color;

void main() {
    frag_tex_coord = tex_coord;
    frag_color = color/255.0;
    gl_Position = vec4(position.x/(screen_dims.x*0.5), position.y/(screen_dims.y*0.5), 0.0, 1.0);
}
