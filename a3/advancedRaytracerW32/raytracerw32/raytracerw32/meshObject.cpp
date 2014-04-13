#include <iostream>
#include <cstdlib>
#include <cmath>
#include "meshObject.h"

#define MAX_LINE_SIZE 256
// #define QUADMESH
// #define INCLUDE_TEXTURE
#define NORMALMESH
MeshObject::MeshObject(const char* filename)
{
    readObj(filename);

}

MeshObject::~MeshObject()
{
    //free(objectName);
}

/*
    Private stuff
*/

void MeshObject::readObj(const char* objName)
{
    FILE *obj = fopen(objName,"r");
    char buffer[MAX_LINE_SIZE];

    if (obj)
    {

        while (fgets(buffer,MAX_LINE_SIZE,obj))
            getNextElement(buffer);
    }
}

void MeshObject::getNextElement(const char* buffer)
{
    switch(buffer[0])
    {
        case PARSED_ELEMENT_VERTEX:
            switch(buffer[1])
            {
                case PARSED_ELEMENT_TEXTURE:
                    readTexCoord(buffer);
                    break;

                case PARSED_ELEMENT_NORMAL:
                    readNormal(buffer);
                    break;

                default:
                    readVertex(buffer);
            }
            break;

        case PARSED_ELEMENT_FACE:
#ifdef QUADMESH
            readFace4(buffer);
#endif
#ifdef INCLUDE_TEXTURE
			readFaceTexture(buffer);
#endif
#ifdef NORMALMESH
			readFace(buffer);
#endif
            break;

        //ignore
        case PARSED_ELEMENT_SMOOTH:
        case PARSED_ELEMENT_COMMENT:
        case PARSED_ELEMENT_MATERIAL:
        default:
            //printf("%s",buffer);
            break;
    }
}

void MeshObject::readVertex(const char* objString)
{
    Vector3f vec3f;
    sscanf(objString,"%*s %f %f %f\n",vec3f.coords,vec3f.coords+1,vec3f.coords+2);

    verts.push_back(vec3f);
}

void MeshObject::readNormal(const char* objString)
{
    Vector3f vec3f;
    sscanf(objString,"%*s %f %f %f\n",vec3f.coords,vec3f.coords+1,vec3f.coords+2);

    normals.push_back(vec3f);
}

void MeshObject::readTexCoord(const char* objString)
{
    Vector3f vec3f;
    sscanf(objString,"%*s %f %f %f\n",vec3f.coords,vec3f.coords+1,vec3f.coords+2);

    texCoords.push_back(vec3f);
}

void MeshObject::readFace(const char* objString)
{
    Face f;
    //this looks worse than I expected
    sscanf(objString,"f %d//%d %d//%d %d//%d\n",
        f.verts.coords,f.normals.coords,
        f.verts.coords+1,f.normals.coords+1,
        f.verts.coords+2,f.normals.coords+2);
    normalizeFaceIndexes(f);
    faces.push_back(f);
}

void MeshObject::readFaceTexture(const char* objString)
{
	Face f;
	Face f2;
	//this looks worse than I expected
	sscanf(objString, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n",
		f.verts.coords, f.texCoords.coords, f.normals.coords,
		f.verts.coords + 1, f.texCoords.coords + 1, f.normals.coords + 1,
		f.verts.coords + 2, f.texCoords.coords + 2, f.normals.coords + 2, 
		f2.verts.coords, f2.texCoords.coords, f2.normals.coords);
		f2.verts.coords[1] = f.verts.coords[0];
		f2.verts.coords[2] = f.verts.coords[2];
		f2.normals.coords[1] = f.normals.coords[0];
		f2.normals.coords[2] = f.normals.coords[2];
		f2.texCoords.coords[1] = f.texCoords.coords[0];
		f2.texCoords.coords[2] = f.texCoords.coords[2];
	normalizeFaceIndexes(f);
	normalizeFaceIndexes(f2);
	faces.push_back(f);
	faces.push_back(f2);
}


void MeshObject::readFace4(const char* objString)
{
	Face f;
	int a, b, c, d;
	//this looks worse than I expected
	int obj_len = strlen(objString);
	int slash_num = 0;
	for (size_t i = 0; i < obj_len; i++)
	{
		if (objString[i] == '/')
			slash_num++;
	}
	if (slash_num == 4)
	{
		Face f2;
		sscanf(objString, "f %d/%d %d/%d %d/%d %d/%d\n",
			f.verts.coords, f.normals.coords,
			f.verts.coords + 1, f.normals.coords + 1,
			f.verts.coords + 2, f.normals.coords + 2,
			f2.verts.coords, f2.normals.coords);
		f2.verts.coords[1] = f.verts.coords[0];
		f2.verts.coords[2] = f.verts.coords[2];
		f2.normals.coords[1] = f.normals.coords[0];
		f2.normals.coords[2] = f.normals.coords[2];
		normalizeFaceIndexes(f2);
		normalizeFaceIndexes(f);
		faces.push_back(f);
		faces.push_back(f2);
	}
	else
	{
		sscanf(objString, "f %d/%d %d/%d %d/%d\n",
			f.verts.coords, f.normals.coords,
			f.verts.coords + 1, f.normals.coords + 1,
			f.verts.coords + 2, f.normals.coords + 2);
		normalizeFaceIndexes(f);
		faces.push_back(f);
	}

}


void MeshObject::normalizeFaceIndexes(Face& f)
{
    for(int i=0; i < 3; ++i)
    {
        f.verts.coords[i] -= 1;
        f.normals.coords[i] -= 1;
        f.texCoords.coords[i] -= 1;
    }
}

void MeshObject::normalizeVectorCoord()
{

    double max_scale = 1e-5;

    for(int i=0; i < verts.size(); ++i)
    {
        max_scale = fmax(max_scale, fabs(verts[i].coords[0]));
        max_scale = fmax(max_scale, fabs(verts[i].coords[1]));
        max_scale = fmax(max_scale, fabs(verts[i].coords[2]));
    }
    std::cout<<"The object will be shrinked by factor "<<max_scale<<std::endl; 
    for(int i=0; i < verts.size(); ++i)
    {
        verts[i].coords[0] = verts[i].coords[0] / max_scale;
        verts[i].coords[1] = verts[i].coords[1] / max_scale;
        verts[i].coords[2] = verts[i].coords[2] / max_scale;
    }

}


