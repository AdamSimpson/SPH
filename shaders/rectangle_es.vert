attribute vec2 position;
uniform vec4 color;

void main() {
    gl_Position = vec4(position, 0.0, 1.0);
}

