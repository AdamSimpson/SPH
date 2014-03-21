attribute vec4 coord;
attribute vec3 color;

varying vec2 frag_tex_coord;
varying vec3 text_color;



void main() {
    frag_tex_coord = coord.zw;
    text_color = color;
    gl_Position = vec4(coord.xy, 0, 1);
}
