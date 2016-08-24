#pragma once

#include <GL/glew.h>

#include "GL/glus.h"
#include "../vsgl3/glHelper.h"

#include <memory>
#include <vector>


struct ClothMeshData
{
  ClothMeshData() : pinPong(false) {}

  std::shared_ptr<SimpleMesh>   pMesh; // mesh to draw points and lines - i.e. springs; 
  std::shared_ptr<SimpleMesh>   pTris; // auxilary index buffer and VAO to render triangles. Pease don't copy VBOs. :)

  std::vector<float4> vertPos0;        // per vertex position; use float4 for pos because SimpleMesh and glusShape use float4 vert pos.
  std::vector<float4> vertVel0;        // per vertex velocity; use float4 also to simplyfy things
	std::vector<int> vertTr;
	std::vector<float> vertTexCoord;
  
  std::vector<float4> vertPos1;        // a copy for pin pong
  std::vector<float4> vertVel1;        // a copy for pin pong

  bool pinPong;                        // (vertPos0 => vertPos1) or (vertPos1 => vertPos0)

  std::vector<float4> vertForces;      // it is possible allocate and free forces each sim step, but we prever to hold it here;
  std::vector<float>  vertMassInv;     // per vertex inverse mass;
  std::vector<float3> vertNormals;     // this is only for graphics; don't need it for sim;

  std::vector<int>    edgeIndices;     // indices of connected verts (with springs); double of connection number size;
  std::vector<float>  edgeHardness;    // Hooke coeff for springs 
  std::vector<float>  edgeInitialLen;  // 
  std::vector<float4> faceNormals;     // this is only for graphics; don't need it for sim;

  inline size_t vertexNumber()     const { return vertPos0.size(); }
  inline size_t connectionNumber() const { return edgeIndices.size() / 2; }

  void initGPUData(GLuint a_programId)   { pMesh->CreateGPUData(a_programId); pTris->CreateGPUData(a_programId); } // we have to call it before we are going to draw our cloth
  void updatePositionsGPU();
  void updateNormalsGPU();

  float4 g_wind;
	float ov_time;
	float ping_wind;
	bool p_wind;
};


ClothMeshData CreateCloth();

void SimStep(ClothMeshData* pMesh, float delta_t);
void RecalculateNormals(ClothMeshData* pMesh);
