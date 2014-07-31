#version 150 core

varying vec2 frag_tex_coord;
varying vec2 blur_tex_coords[14];

uniform sampler2D tex;
 
void main()
{
    vec4 color = vec4(0.0);
    color += texture2D(tex, blur_tex_coords[ 0])*0.0044299121055113265;
    color += texture2D(tex, blur_tex_coords[ 1])*0.00895781211794;
    color += texture2D(tex, blur_tex_coords[ 2])*0.0215963866053;
    color += texture2D(tex, blur_tex_coords[ 3])*0.0443683338718;
    color += texture2D(tex, blur_tex_coords[ 4])*0.0776744219933;
    color += texture2D(tex, blur_tex_coords[ 5])*0.115876621105;
    color += texture2D(tex, blur_tex_coords[ 6])*0.147308056121;
    color += texture2D(tex, frag_tex_coord     )*0.159576912161;
    color += texture2D(tex, blur_tex_coords[ 7])*0.147308056121;
    color += texture2D(tex, blur_tex_coords[ 8])*0.115876621105;
    color += texture2D(tex, blur_tex_coords[ 9])*0.0776744219933;
    color += texture2D(tex, blur_tex_coords[10])*0.0443683338718;
    color += texture2D(tex, blur_tex_coords[11])*0.0215963866053;
    color += texture2D(tex, blur_tex_coords[12])*0.00895781211794;
    color += texture2D(tex, blur_tex_coords[13])*0.0044299121055113265;
    gl_FragColor = color;
}
