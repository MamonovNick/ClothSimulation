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
    vec3  l2   = normalize ( g_sunDir );
    vec3  v2   = normalize ( fragmentView );

    float nl    = dot ( fragmentNormal, l2 );
    float nv    = dot ( fragmentNormal, v2 );
    vec3  lProj = normalize ( l2 - fragmentNormal * nl );
    vec3  vProj = normalize ( v2 - fragmentNormal * nv );
    float cx    = max ( dot ( lProj, vProj ), 0.0 );

    float cosAlpha = nl > nv ? nl : nv;
    float cosBeta  = nl > nv ? nv : nl;
    float dx       = sqrt ( ( 1.0 - cosAlpha * cosAlpha ) * ( 1.0 - cosBeta * cosBeta ) ) / cosBeta;

	if (normals==1){
		fragColor.x = fragmentNormal.x;
		fragColor.y = fragmentNormal.y;
		fragColor.z = fragmentNormal.z;
		fragColor.w = 0;
	} 
	else
		fragColor = texture(diffuseTexture_CL, fragmentTexCoord/18)*max ( 0.0, nl ) * (1.8 + 0.00005 * cx * dx);
}