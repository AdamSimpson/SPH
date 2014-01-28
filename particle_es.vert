attribute vec2 position;
attribute vec3 color;
varying vec3 frag_color;
uniform float radius;

void main() {
   gl_Position = vec4(position, 0.0, 1.0);
   gl_PointSize = radius;
   frag_color = color;
}

