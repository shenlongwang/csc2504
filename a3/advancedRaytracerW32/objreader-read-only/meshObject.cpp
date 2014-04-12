#include "meshObject.h"

#define MAX_LINE_SIZE 256

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
            readFace(buffer);
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
    sscanf(objString,"f %d/%d/%d %d/%d/%d %d/%d/%d\n",
        f.verts.coords,f.normals.coords,f.texCoords.coords,
        f.verts.coords+1,f.normals.coords+1,f.texCoords.coords+1,
        f.verts.coords+2,f.normals.coords+2,f.texCoords.coords+2);

    normalizeFaceIndexes(f);
    faces.push_back(f);

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

