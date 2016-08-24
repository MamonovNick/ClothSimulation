#include "ClothSim.h"
#include <cstdint>
#define HR 16
#define WR 20

#define initialHardnes 4.0f
#define G_Const 9.80665f
#define Res_Const 0.4f

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float dist(float4 a, float4 b)
{
	return sqrtf((b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y) + (b.z - a.z)*(b.z - a.z));
}

ClothMeshData CreateCloth()
{
  ClothMeshData mdata;
	float edgeLength = 0.05f;

	mdata.ov_time = 0;
	mdata.ping_wind = 0;
	mdata.p_wind = false;
  mdata.vertPos0.resize(HR*WR);
	mdata.vertVel0.resize(HR*WR);
	mdata.vertPos1.resize(HR*WR);
	mdata.vertVel1.resize(HR*WR);
	mdata.vertMassInv.resize(HR*WR);
	mdata.vertForces.resize(HR*WR);
	mdata.vertNormals.resize(HR*WR);
	mdata.faceNormals.resize(HR*WR*6);

	float4 init_v = float4(0, 0.5, 0, 1);
	float4 ex_v;

	for (int i = 0; i < HR; i++)
	{
		ex_v = init_v;
		for (int j = 0; j < WR; j++)
		{
			mdata.vertPos0[i*WR + j] = ex_v;
			mdata.vertVel0[i*WR + j] = float4(0, 0, 0, 0);
			mdata.vertMassInv[i*WR + j] = 0.009f;
			ex_v.x += edgeLength;
		}
		init_v.y -= edgeLength;
	}
	
	for (int i = 0; i < mdata.vertPos0.size(); i++)
	{
		for (int j = i + 1; j < mdata.vertPos0.size(); j++)
		{
			float4 vA = mdata.vertPos0[i];
			float4 vB = mdata.vertPos0[j];
			float distance = dist(vA, vB);
			if (distance < 2.0f*sqrtf(2.0f)*edgeLength)
			{
				float hardness = initialHardnes * (edgeLength / distance);
				mdata.edgeIndices.push_back(i);
				mdata.edgeIndices.push_back(j);
				mdata.edgeHardness.push_back(hardness);
				mdata.edgeInitialLen.push_back(distance);
			}
		}
	}

	for (int i = 0; i < HR - 1; i++)
	{
		for (int j = 0; j < WR - 1; j++)
		{
			mdata.vertTr.push_back(i * WR + j);
			mdata.vertTr.push_back(i * WR + j + 1);
			mdata.vertTr.push_back((i + 1)*WR + j);
			mdata.vertTr.push_back(i * WR + j + 1);
			mdata.vertTr.push_back((i + 1)*WR + j + 1);
			mdata.vertTr.push_back((i + 1)*WR + j);
		}
	}

	for (int i = 0; i < HR; i++)
	{
		for (int j = 0; j < WR; j++)
		{
			mdata.vertTexCoord.push_back(j);
			mdata.vertTexCoord.push_back(-i*WR/HR);
		}
	}

  mdata.vertPos1 = mdata.vertPos0;
  mdata.vertVel1 = mdata.vertVel0;

  mdata.g_wind = float4(0.1, 0, 0.1, 0);

  // you can use any intermediate mesh representation or load data to GPU (in VBOs) here immediately.                              <<===== !!!!!!!!!!!!!!!!!!

	mdata.pMesh = std::make_shared<SimpleMesh>();
	mdata.pTris = std::make_shared<SimpleMesh>();

	GLUSshape& shape = mdata.pMesh->m_glusShape;

	shape.numberVertices = mdata.vertPos0.size();
	shape.numberIndices = mdata.edgeIndices.size();

	shape.vertices = (GLUSfloat*)malloc(4 * shape.numberVertices * sizeof(GLUSfloat));
	shape.indices = (GLUSuint*)malloc(shape.numberIndices * sizeof(GLUSuint));

	memcpy(shape.vertices, &mdata.vertPos0[0], sizeof(float) * 4 * shape.numberVertices);
	memcpy(shape.indices, &mdata.edgeIndices[0], sizeof(int) * shape.numberIndices);

  // create graphics mesh; SimpleMesh uses GLUS Shape to store geometry; 
  // we copy data to GLUS Shape, and then these data will be copyed later from GLUS shape to GPU 
  //

	GLUSshape& shape_Tr = mdata.pTris->m_glusShape;

	shape_Tr.numberVertices = mdata.vertPos0.size();
	shape_Tr.numberIndices = mdata.vertTr.size();

	shape_Tr.vertices = (GLUSfloat*)malloc(4 * shape_Tr.numberVertices * sizeof(GLUSfloat));
	shape_Tr.indices = (GLUSuint*)malloc(shape_Tr.numberIndices * sizeof(GLUSuint));
	shape_Tr.texCoords = (GLUSfloat*)malloc(mdata.vertTexCoord.size()*sizeof(GLUSfloat));
	shape_Tr.normals = (GLUSfloat*)malloc(shape_Tr.numberVertices * 3 * sizeof(GLUSfloat));

	memcpy(shape_Tr.vertices, &mdata.vertPos0[0], sizeof(float) * 4 * shape_Tr.numberVertices);
	memcpy(shape_Tr.indices, &mdata.vertTr[0], sizeof(int) * shape_Tr.numberIndices);
	memcpy(shape_Tr.texCoords, &mdata.vertTexCoord[0], sizeof(float) * mdata.vertTexCoord.size());
	memcpy(shape_Tr.normals, &mdata.vertNormals[0], sizeof(float) * 3 * shape_Tr.numberVertices);

	
	//shape.texCoords = (GLUSfloat*)malloc(sizeof(int) * mdata.vertTexCoord.size());
	//memcpy(shape.texCoords, &mdata.vertTri[0], sizeof(int) * mdata.vertTri.size());
	
  // for tri mesh you will need normals, texCoords and different indices
  // 

  return mdata;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClothMeshData::updatePositionsGPU()
{
	if (pMesh == nullptr)
		return;

	glBindBuffer(GL_ARRAY_BUFFER, pMesh->m_vertexPosBufferObject); 
	CHECK_GL_ERRORS;
	glBufferData(GL_ARRAY_BUFFER, vertexNumber() * 4 * sizeof(GLfloat), (pinPong) ? ((GLfloat*)&(vertPos1[0])) : ((GLfloat*)&(vertPos0[0])), GL_STATIC_DRAW); 
	CHECK_GL_ERRORS;

	glBindBuffer(GL_ARRAY_BUFFER, pTris->m_vertexPosBufferObject);
	CHECK_GL_ERRORS;
	glBufferData(GL_ARRAY_BUFFER, vertexNumber() * 4 * sizeof(GLfloat), (GLfloat*)&vertPos0[0], GL_STATIC_DRAW);
	CHECK_GL_ERRORS;
}

void ClothMeshData::updateNormalsGPU()
{
  if (pMesh == nullptr || this->faceNormals.size() == 0)
    return;

	glBindBuffer(GL_ARRAY_BUFFER, pTris->m_vertexNormBufferObject);
	CHECK_GL_ERRORS;
	glBufferData(GL_ARRAY_BUFFER, vertexNumber() * 3 * sizeof(GLfloat), (GLfloat*)&vertNormals[0], GL_STATIC_DRAW);
	CHECK_GL_ERRORS;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SimStep(ClothMeshData* pMesh, float delta_t)
{
 
  // get in and out pointers
  //
  float4* inVertPos  = pMesh->pinPong ? &pMesh->vertPos1[0] : &pMesh->vertPos0[0];
  float4* inVertVel  = pMesh->pinPong ? &pMesh->vertVel1[0] : &pMesh->vertVel0[0];

  float4* outVertPos = pMesh->pinPong ? &pMesh->vertPos0[0] : &pMesh->vertPos1[0];
  float4* outVertVel = pMesh->pinPong ? &pMesh->vertVel0[0] : &pMesh->vertVel1[0];

  // accumulate forces first
  //

	float4 G_Vec = G_Const*float4(0, -1, 0, 0);
	
	for (size_t i = 0; i < pMesh->vertForces.size(); i++)
	{
		pMesh->vertForces[i] = pMesh->vertMassInv[i] * G_Vec - Res_Const * (inVertVel[i]) + pMesh->g_wind;
	}

	if (pMesh->p_wind)
	{
		if (pMesh->ping_wind <= 0.5)
		{
			pMesh->ping_wind += delta_t;
		}
		else
		{
			pMesh->p_wind = false;
			pMesh->ping_wind = 0;
			pMesh->ov_time = 0;
			pMesh->g_wind = float4(0, 0, 0.02, 0);
		}
	}
	else
	{
		if (pMesh->ov_time >= 10)
		{
			pMesh->g_wind = float4(0, 0, 0, 0);
			pMesh->p_wind = true;
		}
		else
		{
			
			for (int i = 5; i < HR; i++)
			{
				for (int j = 4; j < WR-4; j++)
				{
					pMesh->vertForces[i*WR + j] += float4(0.05, 0, 0.05, 0)*sin(pMesh->ov_time);
				}
			}
			for (int i = 0; i < HR; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					pMesh->vertForces[i*WR + j] += float4(0.09, 0, 0.09, 0)*cos(pMesh->ov_time);
				}
			}
			for (int i = 0; i < HR; i++)
			{
				for (int j = WR-4; j < WR; j++)
				{
					pMesh->vertForces[i*WR + j] += -float4(0.09, 0, 0.09, 0)*cos(pMesh->ov_time);
				}
			}
			pMesh->ov_time += delta_t;
		}
	}

  for (int connectId = 0; connectId < pMesh->connectionNumber(); connectId++)
  {
		int i = pMesh->edgeIndices[2*connectId];
		int j = pMesh->edgeIndices[2*connectId + 1];
		float dl = dist(inVertPos[i], inVertPos[j]) - pMesh->edgeInitialLen[connectId];
		float4 f_el = normalize(inVertPos[j] - inVertPos[i]) * pMesh->edgeHardness[connectId] * dl;
		pMesh->vertForces[i] = pMesh->vertForces[i] + f_el;
		pMesh->vertForces[j] = pMesh->vertForces[j] - f_el;
  }

  // update positions and velocity
  //

	for (size_t i = WR; i < pMesh->vertPos0.size(); i++) 
	{
		float4 a = pMesh->vertForces[i] / pMesh->vertMassInv[i];
		outVertVel[i] = inVertVel[i] + a * delta_t;
		outVertPos[i] = inVertPos[i] + inVertVel[i] * delta_t + a * delta_t * delta_t / 2;
	}

  pMesh->pinPong = !pMesh->pinPong; // swap pointers for next sim step
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float3 make3(float4 a)
{
	return float3(a.x, a.y, a.z);
}

void RecalculateNormals(ClothMeshData* pMesh)
{
	int k = 0;
	int g = 0;
	float4* inVertPos = pMesh->pinPong ? &pMesh->vertPos1[0] : &pMesh->vertPos0[0];

	pMesh->faceNormals[k] = cross(inVertPos[1] - inVertPos[0], inVertPos[WR] - inVertPos[0]);
	k++;
	pMesh->vertNormals[g] = make3(pMesh->faceNormals[k-1]);
	g++;

	for (int j = 1; j < WR-1; j++)
	{
		pMesh->faceNormals[k] = cross(inVertPos[j + 1] - inVertPos[j], inVertPos[WR+j] - inVertPos[j]);
		k++;
		pMesh->faceNormals[k] = cross(inVertPos[WR + j - 1] - inVertPos[j], inVertPos[j - 1] - inVertPos[j]);
		k++;
		pMesh->faceNormals[k] = cross(inVertPos[WR + j] - inVertPos[j], inVertPos[WR + j - 1] - inVertPos[j]);
		k++;
		pMesh->vertNormals[g] = make3((pMesh->faceNormals[k - 1] + pMesh->faceNormals[k - 1 -1] + pMesh->faceNormals[k - 1 -2])/3);
		g++;
	}

	pMesh->faceNormals[k] = cross(inVertPos[WR + WR - 1] - inVertPos[WR - 1], inVertPos[WR + WR - 2] - inVertPos[WR - 1]);
	k++;
	pMesh->faceNormals[k] = cross(inVertPos[WR + WR - 2] - inVertPos[WR - 1], inVertPos[WR - 2] - inVertPos[WR - 1]);
	k++;
	pMesh->vertNormals[g] = make3((pMesh->faceNormals[k - 1] + pMesh->faceNormals[k - 1 - 1]) / 2);
	g++;

	for (int i = 1; i < HR-1; i++)
	{
		for (int j = 0; j < WR; j++)
		{
			if (j == 0)
			{
				pMesh->faceNormals[k] = cross(inVertPos[(i - 1)*WR] - inVertPos[i*WR], inVertPos[(i - 1)*WR + 1] - inVertPos[i*WR]);
				k++;
				pMesh->faceNormals[k] = cross(inVertPos[(i - 1)*WR + 1] - inVertPos[i*WR], inVertPos[i*WR + 1] - inVertPos[i*WR]);
				k++;
				pMesh->faceNormals[k] = cross(inVertPos[i*WR + 1] - inVertPos[i*WR], inVertPos[(i + 1)*WR] - inVertPos[i*WR]);
				k++;
				pMesh->vertNormals[g] = make3((pMesh->faceNormals[k - 1] + pMesh->faceNormals[k - 1 - 1] + pMesh->faceNormals[k - 1 - 2]) / 3);
				g++;
				continue;
			}
			if (j == WR - 1)
			{
				pMesh->faceNormals[k] = cross(inVertPos[i*WR + j - 1] - inVertPos[i*WR + j], inVertPos[(i - 1)*WR + j] - inVertPos[i*WR + j]);
				k++;
				pMesh->faceNormals[k] = cross(inVertPos[(i + 1)*WR + j - 1] - inVertPos[i*WR + j], inVertPos[i*WR + j - 1] - inVertPos[i*WR + j]);
				k++;
				pMesh->faceNormals[k] = cross(inVertPos[(i + 1)*WR + j] - inVertPos[i*WR + j], inVertPos[(i + 1)*WR + j - 1] - inVertPos[i*WR + j]);
				k++;
				pMesh->vertNormals[g] = make3((pMesh->faceNormals[k - 1] + pMesh->faceNormals[k - 1 - 1] + pMesh->faceNormals[k - 1 - 2]) / 3);
				g++;
				continue;
			}
			pMesh->faceNormals[k] = cross(inVertPos[(i - 1)*WR + j] - inVertPos[i*WR + j], inVertPos[(i - 1)*WR + j + 1] - inVertPos[i*WR + j]);
			k++;
			pMesh->faceNormals[k] = cross(inVertPos[(i - 1)*WR + j + 1] - inVertPos[i*WR + j], inVertPos[i*WR + j + 1] - inVertPos[i*WR + j]);
			k++;
			pMesh->faceNormals[k] = cross(inVertPos[i*WR + j + 1] - inVertPos[i*WR + j], inVertPos[(i + 1)*WR + j] - inVertPos[i*WR + j]);
			k++;
			pMesh->faceNormals[k] = cross(inVertPos[i*WR + j - 1] - inVertPos[i*WR + j], inVertPos[(i - 1)*WR + j] - inVertPos[i*WR + j]);
			k++;
			pMesh->faceNormals[k] = cross(inVertPos[(i + 1)*WR + j - 1] - inVertPos[i*WR + j], inVertPos[i*WR + j - 1] - inVertPos[i*WR + j]);
			k++;
			pMesh->faceNormals[k] = cross(inVertPos[(i + 1)*WR + j] - inVertPos[i*WR + j], inVertPos[(i + 1)*WR + j - 1] - inVertPos[i*WR + j]);
			k++;
			pMesh->vertNormals[g] = make3((pMesh->faceNormals[k - 1] + pMesh->faceNormals[k - 1 - 1] + pMesh->faceNormals[k - 1 - 2] + pMesh->faceNormals[k - 1 - 3] + pMesh->faceNormals[k - 1 - 4] + pMesh->faceNormals[k - 1 - 5]) / 6);
			g++;
		}
	}

	pMesh->faceNormals[k] = cross(inVertPos[(HR - 2)*WR] - inVertPos[(HR - 1)*WR], inVertPos[(HR - 2)*WR + 1] - inVertPos[(HR - 1)*WR]);
	k++;
	pMesh->faceNormals[k] = cross(inVertPos[(HR - 2)*WR + 1] - inVertPos[(HR - 1)*WR], inVertPos[(HR - 1)*WR + 1] - inVertPos[(HR - 1)*WR]);
	k++;
	pMesh->vertNormals[g] = make3((pMesh->faceNormals[k - 1] + pMesh->faceNormals[k - 1 - 1]) / 2);
	g++;

	for (int j = 1; j < WR - 1; j++)
	{
		pMesh->faceNormals[k] = cross(inVertPos[(HR - 2)*WR + j] - inVertPos[(HR - 1)*WR + j], inVertPos[(HR - 2)*WR + j + 1] - inVertPos[(HR - 1)*WR + j]);
		k++;
		pMesh->faceNormals[k] = cross(inVertPos[(HR - 2)*WR + j + 1] - inVertPos[(HR - 1)*WR + j], inVertPos[(HR - 1)*WR + j + 1] - inVertPos[(HR - 1)*WR + j]);
		k++;
		pMesh->faceNormals[k] = cross(inVertPos[(HR - 1)*WR + j - 1] - inVertPos[(HR - 1)*WR + j], inVertPos[(HR - 2)*WR + j] - inVertPos[(HR - 1)*WR + j]);
		k++;
		pMesh->vertNormals[g] = make3((pMesh->faceNormals[k - 1] + pMesh->faceNormals[k - 1 - 1] + pMesh->faceNormals[k - 1 - 2]) / 3);
		g++;
	}

	pMesh->faceNormals[k] = cross(inVertPos[(HR - 1)*WR + WR - 2] - inVertPos[(HR - 1)*WR + WR - 1], inVertPos[(HR - 2)*WR + WR - 1] - inVertPos[(HR - 1)*WR + WR - 1]);
	k++;
	pMesh->vertNormals[g] = make3(pMesh->faceNormals[k - 1]);
	g++;

}