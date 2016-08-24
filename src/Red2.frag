#version 330

uniform sampler2D diffuseTexture_CL;
uniform int normals;
uniform vec3 g_sunDir;

in  vec3 fragmentWorldPos;
in  vec3 fragmentNormal;
in  vec2 fragmentTexCoord;
in  vec3 fragmentView;

out vec4 fragColor;

void main(void)
{
	if (normals==1){
		fragColor.x = fragmentNormal.x;
		fragColor.y = fragmentNormal.y;
		fragColor.z = fragmentNormal.z;
		fragColor.w = 0;
	}
	else
		fragColor = texture(diffuseTexture_CL, fragmentTexCoord/18)*(max(dot(fragmentNormal,g_sunDir),0)+0.8f);
}