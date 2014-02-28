#version 150 core
in vec2 position;
in vec3 color;
uniform float radius;

out vec3 sphere_color;
out vec2 circle_center;

void main() {
   gl_Position = vec4(position, 0.0, 1.0);
   gl_PointSize = 2.0*radius;

   sphere_color = color;
   circle_center = position;
}

