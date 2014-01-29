#version 150 core
in vec2 position;
in vec3 color;
uniform float radius;
out vec3 frag_color;

void main() {
   gl_Position = vec4(position, 0.0, 1.0);
   gl_PointSize = radius;
   frag_color = color;
}

