#version 330

in FragData
{
    flat vec3 cameraSpherePos;
    flat vec3 sphereColor;
    smooth vec2 mapping;
};

out vec4 outputColor;

void main()
{
    outputColor = vec4(sphereColor, 1.0);
}

