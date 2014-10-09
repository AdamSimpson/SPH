#version 150 core

in vec2 frag_tex_coord;
in vec2 blur_tex_coords[14];

uniform sampler2D tex;
 
out vec4 color;

void main()
{
    float alpha = 0.0;
    alpha += texture(tex, blur_tex_coords[ 0]).a * 0.0044299121055113265;
    alpha += texture(tex, blur_tex_coords[ 1]).a * 0.00895781211794;
    alpha += texture(tex, blur_tex_coords[ 2]).a * 0.0215963866053;
    alpha += texture(tex, blur_tex_coords[ 3]).a * 0.0443683338718;
    alpha += texture(tex, blur_tex_coords[ 4]).a * 0.0776744219933;
    alpha += texture(tex, blur_tex_coords[ 5]).a * 0.115876621105;
    alpha += texture(tex, blur_tex_coords[ 6]).a * 0.147308056121;
    alpha += texture(tex, frag_tex_coord     ).a * 0.159576912161;
    alpha += texture(tex, blur_tex_coords[ 7]).a * 0.147308056121;
    alpha += texture(tex, blur_tex_coords[ 8]).a * 0.115876621105;
    alpha += texture(tex, blur_tex_coords[ 9]).a * 0.0776744219933;
    alpha += texture(tex, blur_tex_coords[10]).a * 0.0443683338718;
    alpha += texture(tex, blur_tex_coords[11]).a * 0.0215963866053;
    alpha += texture(tex, blur_tex_coords[12]).a * 0.00895781211794;
    alpha += texture(tex, blur_tex_coords[13]).a * 0.0044299121055113265;

    color = vec4(0.0, 0.0, 0.0, alpha);
}
