#version 150 core
in vec2 position;
in vec3 color;
uniform float diameter_pixels;

out vec3 sphere_color;
out vec2 circle_center;

void main() {
   gl_Position = vec4(position, 0.0, 1.0);
   gl_PointSize = diameter_pixels;

   sphere_color = color;
   circle_center = position;
}

