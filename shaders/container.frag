#version 330

uniform vec4 color;

in vec4 cameraSpaceFragPos;
in vec4 cameraSpaceNormal;
in vec2 fragTexCoord;

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

float CalcAttenuation(in vec3 cameraSpacePosition,
                      in vec3 cameraSpaceLightPos,
                      out vec3 lightDirection)
{
        vec3 lightDifference =  cameraSpaceLightPos - cameraSpacePosition;
        float lightDistanceSqr = dot(lightDifference, lightDifference);
        lightDirection = lightDifference * inversesqrt(lightDistanceSqr);

        return (1 / ( 1.0 + 1.3 * lightDistanceSqr));
}

void main() {
    Lgt.ambientIntensity= vec4(0.1, 0.1, 0.1, 1.0);
    Lgt.light.cameraSpaceLightPos=worldToCameraMatrix*vec4(0.3, -0.1, -0.4, 1.0);
    Lgt.light.lightIntensity=vec4(0.8, 0.8, 0.8, 1.0);

    vec3 surfaceToLight = vec3(0.0);
    vec3 lightPos = Lgt.light.cameraSpaceLightPos.xyz;
    vec3 fragPos =  cameraSpaceFragPos.xyz;
    float attenIntensity = CalcAttenuation(fragPos, lightPos, surfaceToLight);

    float cosAngleIncidence = dot(normalize(cameraSpaceNormal.xyz), surfaceToLight);
    cosAngleIncidence = clamp(cosAngleIncidence, 0, 1);

    vec4 checker_color = color;

    int num_checks = 30;
    if ((int(floor(num_checks*fragTexCoord.x) + floor(num_checks*fragTexCoord.y)) & 1) == 0) {
        checker_color = vec4(0.3, 0.3, 0.3, 1.0);
    }

    OutColor = (checker_color * Lgt.light.lightIntensity * attenIntensity * cosAngleIncidence) + (checker_color * Lgt.ambientIntensity);
}
