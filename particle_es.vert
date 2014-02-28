attribute vec2 position;
attribute vec3 color;
uniform float radius;

// These should not actually vary
varying vec3 sphere_color;
varying vec2 circle_center;

void main() {
   gl_Position = vec4(position, 0.0, 1.0);
   gl_PointSize = 2.0*radius;
   sphere_color = color;
   circle_center = position;
}

