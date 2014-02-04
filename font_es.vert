attribute vec2 coord;
varying vec2 frag_tex_coord;

void main() {
    frag_tex_coord = coord.zw;
    gl_Position = vec4(coord.xy, 0, 1);
}
