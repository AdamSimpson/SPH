#version 150 core
in vec2 position;
in float particle_id;

uniform float diameter_pixels;
uniform sampler2D rg_tex;

out vec2 rg_color;

void main() {
   gl_Position = vec4(position, 0.0, 1.0);
   gl_PointSize = diameter_pixels;
   rg_color = texture(rg_tex, vec2(particle_id, 0.5)).rg;
}

