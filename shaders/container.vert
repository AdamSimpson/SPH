#version 150 core
in vec3 position;
in vec3 normal;

out vec4 cameraSpaceFragPos;
out vec4 cameraSpaceNormal;

layout(std140) uniform GlobalMatrices
{
    mat4 worldToCameraMatrix;
    mat4 cameraToClip;
};

void main() {
    cameraSpaceFragPos = worldToCameraMatrix*vec4(position, 1.0);
    cameraSpaceNormal = worldToCameraMatrix*vec4(normal, 0.0);
    gl_Position = cameraToClip*cameraSpaceFragPos;
}
