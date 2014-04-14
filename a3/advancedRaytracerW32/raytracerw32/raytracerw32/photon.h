#include "util.h"
#include <vector>

using namespace std;
enum PhotonType { OFF, DIFFUSE, SPECULAR, ABSORPTION, VOLUME };

struct DistanceIndex
{
	double distance;
	int index;
};

struct KDNode{
	size_t index;
	KDNode *left, *right;
};

struct Photon
{
	Photon(){}
	Photon(Point3D p, Colour c, Vector3D d, short f) : pos(p), power(c), dir(d), flag(f){ type = OFF; specular_history = 0; diffuse_history = 0; volume_history = 0; }
	Point3D pos;
	Colour power;
	Vector3D dir;
	short flag;
	PhotonType type;
	unsigned int specular_history;
	unsigned int diffuse_history;
	unsigned int volume_history;
};

class PhotonMap
{
public:
	PhotonMap(){};
	~PhotonMap(){ delete tree_mem; };
	int getSize(){ return photon_map.size(); };
	void addNewPhoton(Point3D pos, Colour col, Vector3D dir, short flag);
	void addExistPhoton(Photon exist_photon);
	void initialPhotonMap(int photon_num);
	Photon getPhoton(int ind){ if (ind < photon_map.size()) return photon_map[ind]; };
	void displayPhotoMap();
	Colour findKNN(Point3D point, int num);
	Colour findKNNMedia(Point3D eye, Vector3D ray, int num);
	bool saveToFile(int width, int height, Point3D eye, Vector3D view,
		Vector3D up, double fov, char* fileName);
	// function related to KDtree
	void initKDtree(int k);
	inline double dist(KDNode*a, Point3D b, int dim);
	inline void swap(KDNode *a, KDNode *b);
	void printKDNode(KDNode *root, int depth);
	KDNode* FindMedian(KDNode *startNode, KDNode *endNode, int ind);
	KDNode* MakeKDTree(KDNode *treeNode, int nodeSize, int ind, int dim);
	void FindNearest(KDNode *rootNode, KDNode *queryNode, int ind,
		int dim, KDNode **bestNodes, double *bestDist);
	void pushNodes(KDNode **resultTree, double *dist, int dim);
	void FindKNearest(KDNode *rootNode, Point3D queryNode, int ind,
		int dim, int k, KDNode **bestNodes, double *bestDist);
	Colour findKNNColorWithKDTree(Point3D point, int num);
private:
	vector<Photon> photon_map;
	KDNode *tree_root;
	KDNode *tree_mem;
	KDNode **bestNodes;
	double *bestDist;
};