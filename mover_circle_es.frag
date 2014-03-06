precision highp float;

varying vec2 frag_tex_coord;

uniform float radius;
uniform vec2 center;
uniform vec3 color;

const vec3 light_intensity = vec3(0.7, 0.7, 0.7);
const vec3 light_position = vec3(0.0, 0.0, 0.8);

void main() {
    // squared 2D distance from center of gl_point
    float rad_squared = dot(frag_tex_coord, frag_tex_coord);

    // If outside of the 2D circle discard
    if(rad_squared > 1.0)
        discard;

    // Calculate 3D normal
    vec3 normal = normalize( vec3(frag_tex_coord, sqrt(1.0 - rad_squared)));

    // GL world coordinates
    vec3 frag_position = (normal * radius) + vec3(center, 0.0);

    // Vector from frag to light
    vec3 frag_to_light = light_position - frag_position;

    // cosine of angle of incidence
    float cosAngleIncidence = clamp( dot(normal, frag_to_light) , 0, 1);

    // diffuse lighting
    vec3 sphere_color = color * light_intensity * cosAngleIncidence;

    // ambient lighting
    sphere_color += color * 0.3;

    // Specular lighting
//    sphere_color += vec3(0.7, 0.7, 0.7)*pow(cosAngleIncidence, 40.0);

    gl_FragColor = vec4(sphere_color, 1.0);
}

