attribute vec2 position;
attribute vec3 color;
varying vec3 frag_color;
void main() {
   gl_Position = vec4(position, 0.0, 1.0);
   gl_PointSize = 5.0;
   frag_color = color;
}

