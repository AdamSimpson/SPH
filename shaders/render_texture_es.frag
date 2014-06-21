varying vec2 frag_tex_coord;

uniform sampler2D tex;

void main() {
    vec4 color = texture2D(tex, frag_tex_coord);
    float alpha = step(0.15, color.a);
    float white;
    white = step(-0.2, -color.a);
    gl_FragColor = vec4(white, white, 10, 0.8*alpha);
}
