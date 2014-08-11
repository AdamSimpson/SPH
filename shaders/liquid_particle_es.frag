void main() {
    // Convert from [0,1] to [-1, 1] coords
    vec2 local_frag_coord = (2.0 * gl_PointCoord) - 1.0;
    local_frag_coord.y = -local_frag_coord.y;

    // squared 2D distance from center of gl_point
    float rad_squared = dot(local_frag_coord, local_frag_coord);

    // If outside of the 2D circle discard
    if(rad_squared > 0.34)
        discard;

    float intensity = 1.0 - 3.0*rad_squared;

    gl_FragColor = vec4(0.0, 0.0, 1.0, intensity);
}

