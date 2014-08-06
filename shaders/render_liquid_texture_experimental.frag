#version 150 core
in vec2 frag_tex_coord;

uniform sampler2D tex;

out vec4 OutColor;

void main() {
    vec4 color = texture(tex, frag_tex_coord);

    if(color.a < 0.1)
        discard;

    float alpha = smoothstep(0.1, 0.11, color.a);

    vec4 base_color = vec4(26.0/255.0, 157.0/255.0, 150.0/255, 0.75*alpha);
    vec4 dark_color  = vec4(0.0/255, 58.0/255, 74.0/255, 1.0);

    OutColor = mix(mix(base_color, dark_color, smoothstep(0.7, 1.0, 1-frag_tex_coord.y)), vec4(1,1,1,1.0), smoothstep(-0.15,-0.1, -color.a));
}
