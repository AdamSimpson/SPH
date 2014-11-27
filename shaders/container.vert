#version 150 core
in vec3 position;

layout(std140) uniform GlobalMatrices
{
    mat4 worldToCameraMatrix;
    mat4 cameraToClip;
};

void main() {
    gl_Position = cameraToClip*worldToCameraMatrix*vec4(position, 1.0);
}
