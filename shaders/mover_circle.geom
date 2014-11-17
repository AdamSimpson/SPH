#version 330

layout(points) in;
layout(triangle_strip, max_vertices=4) out;

uniform mat4 proj;
uniform float sphereRadius;

in VertexData
{
    vec3 cameraSpherePos;
    vec3 sphereColor;
} vert[];

out FragData
{
    flat vec3 cameraSpherePos;
    flat vec3 sphereColor;
    smooth vec2 mapping;
};

const float g_boxCorrection = 1.5;

void main()
{
        vec4 cameraCornerPos;
        //Bottom-left
        mapping = vec2(-1.0, -1.0) * g_boxCorrection;
        cameraSpherePos = vec3(vert[0].cameraSpherePos);
        sphereColor = vert[0].sphereColor;
        cameraCornerPos = vec4(vert[0].cameraSpherePos, 1.0);
        cameraCornerPos.xy += vec2(-sphereRadius, -sphereRadius) * g_boxCorrection;
        gl_Position = proj * cameraCornerPos;
        gl_PrimitiveID = gl_PrimitiveIDIn;
        EmitVertex();

        //Top-left
        mapping = vec2(-1.0, 1.0) * g_boxCorrection;
        cameraSpherePos = vec3(vert[0].cameraSpherePos);
        sphereColor = vert[0].sphereColor;
        cameraCornerPos = vec4(vert[0].cameraSpherePos, 1.0);
        cameraCornerPos.xy += vec2(-sphereRadius, sphereRadius) * g_boxCorrection;
        gl_Position = proj * cameraCornerPos;
        gl_PrimitiveID = gl_PrimitiveIDIn;
        EmitVertex();

        //Bottom-right
        mapping = vec2(1.0, -1.0) * g_boxCorrection;
        cameraSpherePos = vec3(vert[0].cameraSpherePos);
        sphereColor = vert[0].sphereColor;
        cameraCornerPos = vec4(vert[0].cameraSpherePos, 1.0);
        cameraCornerPos.xy += vec2(sphereRadius, -sphereRadius) * g_boxCorrection;
        gl_Position = proj * cameraCornerPos;
        gl_PrimitiveID = gl_PrimitiveIDIn;
        EmitVertex();

        //Top-right
        mapping = vec2(1.0, 1.0) * g_boxCorrection;
        cameraSpherePos = vec3(vert[0].cameraSpherePos);
        sphereColor = vert[0].sphereColor;
        cameraCornerPos = vec4(vert[0].cameraSpherePos, 1.0);
        cameraCornerPos.xy += vec2(sphereRadius, sphereRadius) * g_boxCorrection;
        gl_Position = proj * cameraCornerPos;
        gl_PrimitiveID = gl_PrimitiveIDIn;
        EmitVertex();
}
