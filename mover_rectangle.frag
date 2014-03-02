#version 150 core
in vec3 rec_color;
in vec2 local_frag_coord;

out vec4 out_color;

uniform vec2 center;
uniform vec3 color;

const vec3 light_intensity = vec3(0.7, 0.7, 0.7);
const vec3 light_position = vec3(0.0, 0.0, 0.8);

void main() {
    // Calculate 3D normal
    vec3 normal = normalize( vec3(local_frag_coord, sqrt(1.0 - local_frag_coord.y*local_frag_coord.y)));

    // GL world coordinates
    vec3 frag_position = vec3(gl_FragCoord.x/960.0 - 1.0, gl_FragCoord.y/540.0 - 1.0 ,0.0);

    // Vector from frag to light
    vec3 frag_to_light = light_position - frag_position;

    // cosine of angle of incidence
    float cosAngleIncidence = clamp( dot(normal, frag_to_light) , 0, 1);

    // diffuse lighting
    vec3 rec_color = color * light_intensity * cosAngleIncidence;

    // ambient lighting
    rec_color += color * 0.3;

    // Specular lighting
//    rec_color += vec3(0.7, 0.7, 0.7)*pow(cosAngleIncidence, 40.0);

    out_color = vec4(rec_color, 1.0);
}

