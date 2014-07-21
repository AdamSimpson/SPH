attribute vec2 position;
attribute vec2 tex_coord;

varying vec2 frag_tex_coord;

void main() {
    gl_Position = vec4(position, 0.0, 1.0);
    frag_tex_coord = tex_coord;
}
