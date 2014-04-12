#include "photon.h"
#include "bmp_io.h"
#include <algorithm>    // std::nth_element, std::random_shuffle


bool compare(DistanceIndex disInd1, DistanceIndex disInd2)
{
	return disInd1.distance < disInd2.distance;
}

void PhotonMap::displayPhotoMap(){
	for (int i = 0; i < photon_map.size(); ++i)
	{
		std::cout << i << "th photon, pos: " << photon_map[i].pos << ", col: " << photon_map[i].power << ", dir: " << photon_map[i].dir << std::endl;
	}
}


void PhotonMap::initialPhotonMap(int photon_num){

	Point3D pos(0.0, 0.0, 0.0);
	Colour col(0.5, 0.5, 0.5);
	Vector3D dir(0.0, 0.0, 1.0);
	short flag = 0;
	for (int i = 0; i < photon_num; i++)
	{
		addNewPhoton(pos, col, dir, flag);
	}
	
}

void PhotonMap::addNewPhoton(Point3D pos, Colour col, Vector3D dir, short flag){
	Photon new_photon (pos, col, dir, flag);
	photon_map.push_back(new_photon);
}

void PhotonMap::addExistPhoton(Photon exist_photon){
	photon_map.push_back(exist_photon);
}



bool PhotonMap::saveToFile(int width, int height, Point3D eye, Vector3D view,
	Vector3D up, double fov, char* fileName){
	return false;
}

Colour PhotonMap::findKNN(Point3D point, int num)
{ // find knn
	vector<DistanceIndex> distance;
	Vector3D dist_vector;
	DistanceIndex disind_temp;
	for (int i = 0; i < photon_map.size(); i++)
	{
		dist_vector = point - photon_map[i].pos;
		disind_temp.distance = pow(dist_vector[0], 2) + pow(dist_vector[1], 2) + pow(dist_vector[2], 2);
		disind_temp.index = i;
		distance.push_back(disind_temp);
	}
	std::nth_element(distance.begin(), distance.begin() + num, distance.end(), compare);
	// vector<Photon> knnMap;
	int photon_index;
	Colour mean_color;
	double weight = 0;
	double weight_sum = 0;
	double knn_radius = sqrt(distance[num - 1].distance);
	for (int i = 0; i < num; i++)
	{
		photon_index = distance[i].index;
		// std::cout << distance[i].distance << std::endl;
		weight = 1 - distance[i].distance / (knn_radius * 2);
		mean_color = mean_color + weight * photon_map[photon_index].power;
		weight_sum = weight_sum + weight;
		// knnMap.push_back(photon_map[photon_index]);
	}
	// std::cout << mean_color << std::endl;
	mean_color = (1 / (knn_radius*knn_radius+0.05)) * mean_color;
	// std::cout << mean_color << std::endl;
	return mean_color;

}