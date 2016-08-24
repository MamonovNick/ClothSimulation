#version 330

uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 objectMatrix;

uniform vec3 g_sunDir;

in vec4 vertex;
in vec3 normal;
in vec2 texCoord;

out vec3 fragmentWorldPos;
out vec3 fragmentNormal;
out vec2 fragmentTexCoord;

out vec3 fragmentSunView;
out vec3 fragmentNormalView;
out vec3 fragmentView;

void main(void)
{
  vec4 worldPos    = objectMatrix*vertex;
  vec4 viewPos     = modelViewMatrix*worldPos;
  
  fragmentView     = viewPos.xyz;
  fragmentWorldPos = worldPos.xyz;
  fragmentNormal   = normalize(mat3(objectMatrix)*normal);
  fragmentTexCoord = texCoord;

  fragmentSunView    = mat3(modelViewMatrix)*g_sunDir;
  fragmentNormalView = mat3(modelViewMatrix)*fragmentNormal;

  gl_Position        = projectionMatrix*viewPos;
}

