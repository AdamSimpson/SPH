#version 330

layout(std140) uniform GlobalMatrices
{
    mat4 worldToCameraMatrix;
    mat4 cameraToClip;
};

layout(std140) uniform GlobalLight
{
    vec4 worldSpacePos;
    vec4 cameraSpacePos;
    vec4 intensity;
    vec4 ambientIntensity;
    float attenuation;
} LightData;

in FragData
{
    flat vec3 cameraSpherePos;
    flat vec3 sphereColor;
    smooth vec2 mapping;
};

out vec4 outputColor;

struct material_t
{
    vec4 diffuseColor;
    vec4 specularColor;
    float specularShininess;
} Mtl;

float CalcAttenuation(in vec3 cameraSpacePosition,
	out vec3 lightDirection)
{
	vec3 lightDifference =  LightData.cameraSpacePos.xyz - cameraSpacePosition;
	float lightDistanceSqr = dot(lightDifference, lightDifference);
	lightDirection = lightDifference * inversesqrt(lightDistanceSqr);
	
	return (1 / ( 1.0 + LightData.attenuation * lightDistanceSqr));
}

uniform float sphereRadius;

vec4 ComputeLighting(in vec3 cameraSpacePosition,
	in vec3 cameraSpaceNormal, in material_t material)
{
	vec3 lightDir;
	vec4 lightIntensity;
	if(LightData.cameraSpacePos.w == 0.0)
	{
		lightDir = vec3(LightData.cameraSpacePos);
		lightIntensity = LightData.intensity;
	}
	else
	{
		float atten = CalcAttenuation(cameraSpacePosition, lightDir);
		lightIntensity = atten * LightData.intensity;
	}
	
	vec3 surfaceNormal = normalize(cameraSpaceNormal);
	float cosAngIncidence = dot(surfaceNormal, lightDir);
	cosAngIncidence = cosAngIncidence < 0.0001 ? 0.0 : cosAngIncidence;
	
	vec3 viewDirection = normalize(-cameraSpacePosition);
	
	vec3 halfAngle = normalize(lightDir + viewDirection);
	float angleNormalHalf = acos(dot(halfAngle, surfaceNormal));
	float exponent = angleNormalHalf / material.specularShininess;
	exponent = -(exponent * exponent);
	float gaussianTerm = exp(exponent);

	gaussianTerm = cosAngIncidence != 0.0 ? gaussianTerm : 0.0;
	
	vec4 lighting = material.diffuseColor * lightIntensity * cosAngIncidence;
	lighting += material.specularColor * lightIntensity * gaussianTerm;
	
	return lighting;
}

void Impostor(out vec3 cameraPos, out vec3 cameraNormal)
{
	vec3 cameraPlanePos = vec3(mapping * sphereRadius, 0.0) + cameraSpherePos;
	vec3 rayDirection = normalize(cameraPlanePos);
	
	float B = 2.0 * dot(rayDirection, -cameraSpherePos);
	float C = dot(cameraSpherePos, cameraSpherePos) -
		(sphereRadius * sphereRadius);
	
	float det = (B * B) - (4 * C);
	if(det < 0.0)
		discard;
		
	float sqrtDet = sqrt(det);
	float posT = (-B + sqrtDet)/2;
	float negT = (-B - sqrtDet)/2;
	
	float intersectT = min(posT, negT);
	cameraPos = rayDirection * intersectT;
	cameraNormal = normalize(cameraPos - cameraSpherePos);
}

void main()
{
        Mtl.diffuseColor = vec4(sphereColor, 0.98);
        Mtl.specularColor = vec4(0.8, 0.8, 0.8, 0.98);
        Mtl.specularShininess = 0.1;

	vec3 cameraPos;
	vec3 cameraNormal;
	
	Impostor(cameraPos, cameraNormal);
	
	//Set the depth based on the new cameraPos.
	vec4 clipPos = cameraToClip * vec4(cameraPos, 1.0);
	float ndcDepth = clipPos.z / clipPos.w;
	gl_FragDepth = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;
	
	vec4 accumLighting = Mtl.diffuseColor * LightData.ambientIntensity;
	accumLighting += ComputeLighting(cameraPos, cameraNormal, Mtl);
	
	outputColor = accumLighting;
        outputColor.a = 0.98;
}

