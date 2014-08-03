#version 150 core
out vec4 out_color;

void main() {
    // Convert from [0,1] to [-1, 1] coords
    vec2 local_frag_coord = (2.0 * gl_PointCoord) - 1.0;
    local_frag_coord.y = -local_frag_coord.y;

    // squared 2D distance from center of gl_point
    float rad_squared = dot(local_frag_coord, local_frag_coord);

    // Discard outside of circle
    if(rad_squared > 0.34)
        discard;
    
    // Alpha component, metallball "potential"
    // Simple form to avoid sqrt
    float intensity = 1.0 - 3.0*rad_squared;

    out_color = vec4(0,0,0, intensity);
}

