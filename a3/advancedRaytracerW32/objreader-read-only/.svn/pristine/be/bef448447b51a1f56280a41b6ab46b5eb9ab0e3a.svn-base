#ifndef __MESH_H__
#define __MESH_H__

class Mesh
{
    public:
        Mesh();
        Mesh(const char* filename);
        virtual ~Mesh();

        inline MeshObject* getObjects() { return objects; }

        inline int getSize() { return size; }

    private:
        int size;
        MeshObject *objects;
};

#endif /* __MESH_H__ */

