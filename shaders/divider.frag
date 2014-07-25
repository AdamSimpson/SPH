#version 150 core

in vec3 frag_color;
out vec4 OutColor;

void main() {
    OutColor = vec4(frag_color, 1.0);
}
