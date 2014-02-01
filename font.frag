#version 150 core
in vec2 frag_tex_coord;
in vec4 frag_color;

out vec4 outColor;
uniform sampler2D tex;

void main() {
    float alpha = texture(tex, frag_tex_coord).r;
    outColor = vec4(frag_color.rgb, alpha*frag_color.a);
}

