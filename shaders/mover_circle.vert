#version 150 core

uniform vec3 center;
uniform vec3 color;

layout(std140) uniform GlobalMatrices
{
    mat4 worldToCameraMatrix;
    mat4 cameraToClip;
};

out VertexData
{
    vec3 cameraSpherePos;
    vec3 sphereColor;
} outData;


void main() {
   outData.cameraSpherePos = (worldToCameraMatrix*vec4(center, 1.0)).xyz;
   outData.sphereColor = color;
}

