#version 150 core
in vec3 position;
in vec3 color;
uniform float diameter_pixels;

out vec3 sphere_color;
out vec3 sphere_center;

void main() {
   gl_Position = vec4(position, 1.0);
   gl_PointSize = diameter_pixels;

   sphere_color = color;
   sphere_center = position;
}

