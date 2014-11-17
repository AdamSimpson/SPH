#version 150 core

uniform vec3 center;
uniform vec3 color;

uniform mat4 view;
uniform mat4 proj;

out VertexData
{
    vec3 cameraSpherePos;
    vec3 sphereColor;
} outData;


void main() {
   outData.cameraSpherePos = (view*vec4(center, 1.0)).xyz;
   outData.sphereColor = color;
}

