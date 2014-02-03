precision mediump float;

varying vec2 frag_tex_coord;
varying vec4 frag_color;

uniform sampler2D tex;

void main() {
    float alpha = texture2D(tex, frag_tex_coord).a;
    gl_FragColor = vec4(frag_color.rgb, alpha*frag_color.a);
}

