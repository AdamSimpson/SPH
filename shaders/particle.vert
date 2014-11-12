#version 150 core
in vec3 position;
in vec3 color;
uniform float diameter_pixels;

out vec3 sphere_color;
out vec3 sphere_center;

uniform mat4 view;
uniform mat4 proj;

void main() {
   gl_Position = proj*view*vec4(position, 0.0);
   gl_PointSize = diameter_pixels;

   sphere_color = color;
   sphere_center = position;
}

