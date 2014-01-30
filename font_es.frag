precision mediump float;
varying vec2 frag_tex_coord;
attribute vec4 outColor;
uniform sampler2D tex;
void main() {
    float alpha = texture2D(tex, frag_tex_coord).r;
    outColor = vec4(1.0, 1.0, 1.0, alpha);
}

