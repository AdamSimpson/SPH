precision highp float;
varying vec3 frag_color;
const vec2 center = vec2(0.5, 0.5);
const float radius = 0.5;

void main() {
    float distance_from_center = distance(center, gl_PointCoord);
    float in_radius = step(distance_from_center,radius);
    gl_FragColor = vec4(frag_color, in_radius);
}
