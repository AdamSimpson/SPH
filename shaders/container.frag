#version 330

uniform vec4 color;

in vec4 cameraSpaceFragPos;
in vec4 cameraSpaceNormal;

out vec4 OutColor;

layout(std140) uniform GlobalMatrices
{
    mat4 worldToCameraMatrix;
    mat4 cameraToClip;
};

struct PerLight
{
        vec4 cameraSpaceLightPos;
        vec4 lightIntensity;
};

struct Light
{
        vec4 ambientIntensity;
        float lightAttenuation;
        PerLight light;
} Lgt;

void main() {
    Lgt.ambientIntensity= vec4(0.2, 0.2, 0.2, 1.0);
    Lgt.light.cameraSpaceLightPos=worldToCameraMatrix*vec4(0.4, 0.0, -0.8, 1.0);
    Lgt.light.lightIntensity=vec4(0.6, 0.6, 0.6, 1.0);

    vec3 lightPos = Lgt.light.cameraSpaceLightPos.xyz;
    vec3 fragPos =  cameraSpaceFragPos.xyz;
    vec3 surfaceToLight = normalize(lightPos - fragPos);    

    float cosAngleIncidence = dot(normalize(cameraSpaceNormal.xyz), surfaceToLight);
    cosAngleIncidence = clamp(cosAngleIncidence, 0, 1);

    OutColor = (color * Lgt.light.lightIntensity * cosAngleIncidence) + (color * Lgt.ambientIntensity);
}
