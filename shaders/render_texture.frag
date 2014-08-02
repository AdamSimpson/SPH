#version 150 core
in vec2 frag_tex_coord;

uniform sampler2D tex;

out vec4 OutColor;

void main() {
    vec4 tex_color = texture(tex, frag_tex_coord);
    float alpha = tex_color.a;

    if(alpha < 0.20)
        discard;

    float white = step(-0.15, -alpha);
    float blue = step(0.1, alpha);
    vec4 color = vec4(white,white,blue,0.3);
    if(tex_color.r > 0.95 && tex_color.g < 0.05)
        color = vec4(0.0,0.1,1.0,0.4);

    OutColor = color;    
//    OutColor = vec4(tex_color.rg, 0.0, 1.0);
}
