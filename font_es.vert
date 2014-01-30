attribute vec2 position;
attribute vec2 tex_coord;
varying vec2 frag_tex_coord;

void main() {
    frag_tex_coord = tex_coord;
    gl_Position = vec4(position.x/800, position.y/-800, 0.0, 1.0);
}
