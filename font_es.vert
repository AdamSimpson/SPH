attribute vec2 position;
attribute vec2 tex_coord;
attribute vec4 color;

uniform ivec2 screen_dims;

varying vec2 frag_tex_coord;
varying vec4 frag_color;

void main() {
    frag_tex_coord = tex_coord;
    frag_color = color/255.0;
    gl_Position = vec4(position.x/(screen_dims.x*0.5), position.y/(screen_dims.y*0.5), 0.0, 1.0);
}
