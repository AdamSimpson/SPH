#version 150 core
in vec3 position;
in vec3 color;

uniform mat4 view;

out VertexData
{
    vec3 cameraSpherePos;
    vec3 sphereColor;
} outData;

void main() {
   outData.cameraSpherePos = (view*vec4(position, 1.0)).xyz;
   outData.sphereColor = color;
}

