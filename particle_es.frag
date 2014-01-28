precision mediump float;
varying vec3 frag_color;
const mediump vec2 center = vec2(0.5, 0.5);
const mediump float radius = 0.5;

void main() {
    mediump float distance_from_center = distance(center, gl_PointCoord);
    lowp float in_radius = step(distance_from_center,radius);
    gl_FragColor = vec4(frag_color, in_radius);
}
