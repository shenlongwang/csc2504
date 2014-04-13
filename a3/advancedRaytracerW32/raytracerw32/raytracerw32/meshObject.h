//
//  TODO
//

/*

    Break stuff in separate files.
    Find a way to denote vertex groups.
    Check if input was succesfull?

*/
#ifndef __MESH_OBJECT_H__
#define __MESH_OBJECT_H__

#include <cstdio>
#include <cstdlib>
#include <vector>

using std::vector;

//structure to hold vertex data
typedef struct _VECTOR3F {
    float coords[3];
}Vector3f;

typedef struct _VECTOR3I {
    int coords[3];
}Vector3i;

typedef struct _FACE {
    Vector3i verts;
    Vector3i normals;
    Vector3i texCoords;
}Face;

enum PARSED_ELEMENT {
    PARSED_ELEMENT_COMMENT = '#',
    PARSED_ELEMENT_VERTEX = 'v',
    PARSED_ELEMENT_NORMAL = 'n',
    PARSED_ELEMENT_TEXTURE = 't',
    PARSED_ELEMENT_FACE = 'f',
    PARSED_ELEMENT_GROUP = 'g',
    PARSED_ELEMENT_OBJECT = 'o',
    PARSED_ELEMENT_SMOOTH = 's',
    PARSED_ELEMENT_MATERIAL = 'm'
};

class MeshObject{

   public:
       MeshObject(const char* filename);
       ~MeshObject();

       inline vector<Vector3f>& getVerts() { return verts; }
       inline vector<Face>& getFaces() { return faces; }
       inline vector<Vector3f>& getNormals() { return normals; }
       inline vector<Vector3f>& getTexCoords() { return texCoords; }

       inline int getFaceSides() { return faceSides; }
       void normalizeVectorCoord();

   private:

       void readVertex(const char* objString);
       void readNormal(const char* objString);
       void readTexCoord(const char* objString);
       void readFace(const char* objString);
	   void readFace4(const char* objString);
	   void readFaceTexture(const char* objString);

       void getNextElement(const char* buffer);

       void readObj(const char* filename);
       void normalizeFaceIndexes(Face& f);
       //
       //  DATA
       //

       //all the raw vertex in order
       vector<Vector3f> verts;

       //all the texture coordinates in order
       vector<Vector3f> texCoords;

       //all the normals in order
       vector<Vector3f> normals;

       //the face Indexes
       vector<Face> faces;

       //number representing how much sides the faces have for this object
       int faceSides;

};

#endif /* __MESH_OBJECT_H__ */

