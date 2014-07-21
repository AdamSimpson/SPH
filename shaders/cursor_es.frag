varying vec2 frag_tex_coord;

uniform sampler2D tex;

void main() {
    gl_FragColor = texture2D(tex, frag_tex_coord);
}
