#version 330

uniform sampler2D diffuseTexture_PL;
uniform vec3 g_sunDir;
in vec2 fragmentTexCoord;
out vec4 fragColor;

void main(void)
{	
  fragColor = texture(diffuseTexture_PL, 4*fragmentTexCoord)*(max(g_sunDir.y,0)+1.0f);
}





