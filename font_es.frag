precision mediump float;

varying vec2 frag_tex_coord;

uniform sampler2D tex;

const vec4 color = vec4(1,1,1,1);

void main() {
    gl_FragColor = vec4(1, 1, 1, texture2D(tex, frag_tex_coord).a) * color;
}

