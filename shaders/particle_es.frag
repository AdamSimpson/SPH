precision highp float;

varying vec3 sphere_color;
varying vec2 circle_center;

uniform float radius_world;

const vec3 light_intensity = vec3(0.8, 0.8, 0.8);
const vec3 light_position = vec3(0.0, 0.0, 2.0);

void main() {
    // Convert from [0,1] to [-1, 1] coords
    vec2 local_frag_coord = (2.0 * gl_PointCoord) - 1.0;

    local_frag_coord.y = -local_frag_coord.y;

    // squared 2D distance from center of gl_point
    float rad_squared = dot(local_frag_coord, local_frag_coord);

    // If outside of the 2D circle discard
    if(rad_squared > 1.0)
        discard;

    // Calculate 3D normal
    //vec3 normal = normalize( vec3(local_frag_coord, sqrt(1.0 - rad_squared)));
    vec3 normal = normalize( vec3(local_frag_coord*local_frag_coord, 1.0 - rad_squared));

    // GL world coordinates
    vec3 frag_position = (normal * radius_world) + vec3(circle_center, 0.0);

    // Vector from frag to light
    vec3 frag_to_light = normalize(light_position - frag_position);

    // cosine of angle of incidence
    float cosAngleIncidence = dot(normal, frag_to_light);
    cosAngleIncidence = clamp(cosAngleIncidence, 0.0, 1.0);

    // diffuse lighting
    vec3 color = sphere_color * light_intensity * cosAngleIncidence;

    // ambient lighting
    color += sphere_color * 0.2;

    gl_FragColor = vec4(color, 1.0);
}
