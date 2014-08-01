#version 150 core
in vec2 frag_tex_coord;

uniform sampler2D tex;

out vec4 OutColor;

void main() {
    vec4 tex_color = texture(tex, frag_tex_coord);
    float alpha = tex_color.a;

    if(alpha < 0.1)
        discard;

    float white = step(-0.15, -alpha);

    vec4 color = vec4(0.0,0.0,1.0,0.3);
    if(tex_color.r > 0.9 && tex_color.g > 0.9)
        color = vec4(0.5,0.6,0.8,0.9);

//    OutColor = color;    
    OutColor = vec4(tex_color.rg, 0.0, 1.0);
}
