varying vec2 frag_tex_coord;

uniform sampler2D tex;

void main() {
    vec4 color = texture2D(tex, frag_tex_coord);

    if(color.a < 0.1)
        discard;

    float alpha = smoothstep(0.1, 0.11, color.a);
    float white = smoothstep(-0.15,-0.1, -color.a);

    gl_FragColor = vec4(white, white, 1.0, 0.75*alpha);
}
