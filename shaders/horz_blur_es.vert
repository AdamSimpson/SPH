attribute vec2 position;
attribute vec2 tex_coord;
 
varying vec2 frag_tex_coord;
varying vec2 blur_tex_coords[14];
 
void main()
{
    gl_Position = vec4(position, 0.0, 1.0);
    frag_tex_coord = tex_coord;
    float mult = 1.0;//radius of gaussian blur
    blur_tex_coords[ 0] = frag_tex_coord + mult*vec2(-0.028, 0.0);
    blur_tex_coords[ 1] = frag_tex_coord + mult*vec2(-0.024, 0.0);
    blur_tex_coords[ 2] = frag_tex_coord + mult*vec2(-0.020, 0.0);
    blur_tex_coords[ 3] = frag_tex_coord + mult*vec2(-0.016, 0.0);
    blur_tex_coords[ 4] = frag_tex_coord + mult*vec2(-0.012, 0.0);
    blur_tex_coords[ 5] = frag_tex_coord + mult*vec2(-0.008, 0.0);
    blur_tex_coords[ 6] = frag_tex_coord + mult*vec2(-0.004, 0.0);
    blur_tex_coords[ 7] = frag_tex_coord + mult*vec2( 0.004, 0.0);
    blur_tex_coords[ 8] = frag_tex_coord + mult*vec2( 0.008, 0.0);
    blur_tex_coords[ 9] = frag_tex_coord + mult*vec2( 0.012, 0.0);
    blur_tex_coords[10] = frag_tex_coord + mult*vec2( 0.016, 0.0);
    blur_tex_coords[11] = frag_tex_coord + mult*vec2( 0.020, 0.0);
    blur_tex_coords[12] = frag_tex_coord + mult*vec2( 0.024, 0.0);
    blur_tex_coords[13] = frag_tex_coord + mult*vec2( 0.028, 0.0);
}
