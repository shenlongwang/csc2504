#include "util.h"
#include <vector>

using namespace std;
enum PhotonType { OFF, DIFFUSE, SPECULAR, ABSORPTION };

struct DistanceIndex
{
	double distance;
	int index;
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
	~PhotonMap(){};
	int getSize(){ return photon_map.size(); };
	void addNewPhoton(Point3D pos, Colour col, Vector3D dir, short flag);
	void addExistPhoton(Photon exist_photon);
	void initialPhotonMap(int photon_num);
	Photon getPhoton(int ind){ if (ind < photon_map.size()) return photon_map[ind]; };
	void displayPhotoMap();
	Colour findKNN(Point3D point, int num);
	bool saveToFile(int width, int height, Point3D eye, Vector3D view,
		Vector3D up, double fov, char* fileName);
	//bool balanceMap();
	//bool renderMap();
private:
	vector<Photon> photon_map;
};

