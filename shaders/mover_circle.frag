#version 330

in FragData
{
    flat vec3 cameraSpherePos;
    flat vec3 sphereColor;
    smooth vec2 mapping;
};

out vec4 outputColor;

struct MaterialEntry
{
    vec4 diffuseColor;
    vec4 specularColor;
    float specularShininess;
} Mtl;

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
	
	return (1 / ( 1.0 + Lgt.lightAttenuation * lightDistanceSqr));
}

uniform mat4 proj;
uniform float sphereRadius;

vec4 ComputeLighting(in PerLight lightData, in vec3 cameraSpacePosition,
	in vec3 cameraSpaceNormal, in MaterialEntry material)
{
	vec3 lightDir;
	vec4 lightIntensity;
	if(lightData.cameraSpaceLightPos.w == 0.0)
	{
		lightDir = vec3(lightData.cameraSpaceLightPos);
		lightIntensity = lightData.lightIntensity;
	}
	else
	{
		float atten = CalcAttenuation(cameraSpacePosition,
			lightData.cameraSpaceLightPos.xyz, lightDir);
		lightIntensity = atten * lightData.lightIntensity;
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
        Mtl.diffuseColor = vec4(sphereColor, 1);
        Mtl.specularColor = vec4(0.8, 0.8, 0.8 ,0.8);
        Mtl.specularShininess = 0.1;

        Lgt.ambientIntensity= vec4(0.2, 0.2, 0.2, 1.0);
        Lgt.lightAttenuation = 1.0f / (25.0*25.0);
        Lgt.light.cameraSpaceLightPos=vec4(50.0, 50.0, 0.1, 1.0);
        Lgt.light.lightIntensity=vec4(0.6, 0.6, 0.6, 1.0);

	vec3 cameraPos;
	vec3 cameraNormal;
	
	Impostor(cameraPos, cameraNormal);
	
	//Set the depth based on the new cameraPos.
	vec4 clipPos = proj * vec4(cameraPos, 1.0);
	float ndcDepth = clipPos.z / clipPos.w;
	gl_FragDepth = ((gl_DepthRange.diff * ndcDepth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;
	
	vec4 accumLighting = Mtl.diffuseColor * Lgt.ambientIntensity;
		accumLighting += ComputeLighting(Lgt.light,
			cameraPos, cameraNormal, Mtl);
	
	outputColor = sqrt(accumLighting); //2.0 gamma correction
}

