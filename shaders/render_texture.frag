#version 150 core
in vec2 frag_tex_coord;

uniform sampler2D tex;

out vec4 OutColor;

void main() {
    vec2 coord = vec2(frag_tex_coord.x,-frag_tex_coord.y);
    vec4 color = texture(tex, coord);
    OutColor = (color.a > 0.2) ? vec4(0.0,0.0,1.0,1.0) : vec4(0,0,0,0);;
}
