#version 150 core
in vec2 position;
in vec2 tex_coord;
 
out vec2 frag_tex_coord;
out vec2 blur_tex_coords[14];
 
void main()
{
    gl_Position = vec4(position, 0.0, 1.0);
    frag_tex_coord = tex_coord;
    blur_tex_coords[ 0] = frag_tex_coord + vec2(0.0, -0.028);
    blur_tex_coords[ 1] = frag_tex_coord + vec2(0.0, -0.024);
    blur_tex_coords[ 2] = frag_tex_coord + vec2(0.0, -0.020);
    blur_tex_coords[ 3] = frag_tex_coord + vec2(0.0, -0.016);
    blur_tex_coords[ 4] = frag_tex_coord + vec2(0.0, -0.012);
    blur_tex_coords[ 5] = frag_tex_coord + vec2(0.0, -0.008);
    blur_tex_coords[ 6] = frag_tex_coord + vec2(0.0, -0.004);
    blur_tex_coords[ 7] = frag_tex_coord + vec2(0.0, 0.004);
    blur_tex_coords[ 8] = frag_tex_coord + vec2(0.0, 0.008);
    blur_tex_coords[ 9] = frag_tex_coord + vec2(0.0, 0.012);
    blur_tex_coords[10] = frag_tex_coord + vec2(0.0, 0.016);
    blur_tex_coords[11] = frag_tex_coord + vec2(0.0, 0.020);
    blur_tex_coords[12] = frag_tex_coord + vec2(0.0, 0.024);
    blur_tex_coords[13] = frag_tex_coord + vec2(0.0, 0.028);
}
