attribute vec2 position;
attribute vec2 tex_coord;

uniform float radius;
uniform vec2 center;
uniform vec3 color;

varying vec2 frag_tex_coord;

void main() {
    frag_tex_coord = tex_coord;
    gl_Position = vec4(position, 0.0, 1.0);
}

