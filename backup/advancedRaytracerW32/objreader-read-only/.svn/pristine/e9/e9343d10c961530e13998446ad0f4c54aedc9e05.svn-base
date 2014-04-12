#include <stdio.h>
#include <stdlib.h>
#include "meshObject.h"
#include <vector>

using std::vector;

void dumpvectorfaces(vector<Face> vec)
{
    for (unsigned int i=0; i < vec.size(); i+=3)
    {
        printf("(%i,%i,%i)\n",vec[i].verts.coords[0],vec[i].verts.coords[1],vec[i].verts.coords[2]);
        printf("(%i,%i,%i)\n",vec[i].texCoords.coords[0],vec[i].texCoords.coords[1],vec[i].texCoords.coords[2]);
        printf("(%i,%i,%i)\n",vec[i].normals.coords[0],vec[i].normals.coords[1],vec[i].normals.coords[2]);
    }
}

void dumpvectorf(vector<Vector3f> vec)
{
    for (unsigned int i=0; i < vec.size(); ++i)
    {
        printf("(%f,%f,%f)\n",vec[i].coords[0],vec[i].coords[1],vec[i].coords[2]);
    }
}

int main(int argc,char **argv)
{
    if (argc < 2)
    {
        printf("Screw you I need argumments!\n");
        exit(1);
    }

    MeshObject mesh(argv[1]);

    printf("VERTEX\n");
    dumpvectorf(mesh.getVerts());
    printf("FACES\n");
    dumpvectorfaces(mesh.getFaces());
    printf("NORMALS\n");
    dumpvectorf(mesh.getNormals());
    printf("TEXCOORDS\n");
    dumpvectorf(mesh.getTexCoords());
}

