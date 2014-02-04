precision mediump float;

varying vec2 frag_tex_coord;

uniform sampler2D tex;
uniform vec4 color;

void main() {
    gl_FragColor = vec4(1, 1, 1, texture2D(tex, texcoord).a) * color;
}

