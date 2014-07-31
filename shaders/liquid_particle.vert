#version 150 core
in vec2 position;
uniform float diameter_pixels;

void main() {
   gl_Position = vec4(position, 0.0, 1.0);
   gl_PointSize = diameter_pixels;
}

