varying vec2 frag_tex_coord;

uniform sampler2D tex;

void main() {
    vec4 color = texture(tex, frag_tex_coord);
    float alpha = 0.8*step(0.1, color.a);
    float white;
    white = 0.45*step(-0.15, -color.a);
    gl_FragColor = vec4(white, white, 0.5, alpha);
}
