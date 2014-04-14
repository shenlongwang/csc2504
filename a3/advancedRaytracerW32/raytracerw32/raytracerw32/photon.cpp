#include "photon.h"
#include "bmp_io.h"
#include <algorithm>    // std::nth_element, std::random_shuffle


bool compare(DistanceIndex disInd1, DistanceIndex disInd2)
{
	return disInd1.distance < disInd2.distance;
}

bool compareLarge(DistanceIndex disInd1, DistanceIndex disInd2)
{
	return disInd1.distance > disInd2.distance;
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


Colour PhotonMap::findKNNMedia(Point3D eye, Vector3D ray, int num)
{ // find knn
	vector<DistanceIndex> distance;
	Vector3D dist_vector;
	DistanceIndex disind_temp;
	for (int i = 0; i < photon_map.size(); i++)
	{
		dist_vector = photon_map[i].pos - eye;
		dist_vector.normalize();
		disind_temp.distance = pow(dist_vector.dot(ray), 8);
		disind_temp.index = i;
		distance.push_back(disind_temp);
		// std::cout << disind_temp.distance << std::endl;
	}
	std::nth_element(distance.begin(), distance.begin() + num, distance.end(), compareLarge);
	// vector<Photon> knnMap;
	int photon_index;
	Colour mean_color(0.0, 0.0, 0.0);
	double weight = 0;
	double weight_sum = 0;
	// double knn_radius = pow((distance[num - 1].distance - 0.96) / 0.04, 4);
	double knn_radius = pow((distance[num - 1].distance - 0.985) / 0.015, 4);
	std::cout << knn_radius << std::endl;
	if (double(distance[num - 1].distance < 0.985))
		return mean_color;
	//for (int i = 0; i < num; i++)
	//{
	//	photon_index = distance[i].index;
	//	// weight = pow((distance[i].distance - 0.99)/0.01, 4);
	//	// std::cout << distance[i].distance << std::endl;
	//	// std::cout << weight << std::endl;
	//	mean_color = mean_color + photon_map[photon_index].power;
	//	weight_sum = weight_sum + 1;
	//	// std::cout << distance[i].distance << std::endl;
	//	// knnMap.push_back(photon_map[photon_index]);
	//}
	//// std::cout << mean_color << std::endl;
	//mean_color = (1 / weight_sum + DBL_EPSILON) * (knn_radius)* mean_color;
	mean_color = (knn_radius)*photon_map[distance[num - 1].index].power;
	std::cout << mean_color << std::endl;
	return mean_color;

}


inline double PhotonMap::dist(KDNode *a, Point3D b, int dim)
{
	// Calculate the distance between two Nodes;
	double t, d = 0;
	while (dim--)
	{
		t = photon_map[a->index].pos[dim] - b[dim];
		d += t * t;
	}
	return d;
}

inline void PhotonMap::swap(KDNode *a, KDNode *b)
{
	// Swap the data of the two nodes;
	int temp;
	temp = a->index;
	a->index = b->index;
	b->index = temp;
}

void PhotonMap::printKDNode(KDNode *root, int depth)
{
	std::cout << depth << "th Degree Nodes:" << std::endl;
	std::cout << "(" << photon_map[root->index].pos[0] << ", " << photon_map[root->index].pos[1] << ", " << photon_map[root->index].pos[2] << ")" << std::endl;
	if (depth < 0)
		return;
	if (root->left != NULL)
		printKDNode(root->left, depth + 1);
	if (root->right != NULL)
		printKDNode(root->right, depth + 1);
}

KDNode* PhotonMap::FindMedian(KDNode *startNode, KDNode *endNode, int ind)
{
	// Find the median of a set of nodes, sorted according to INDth value. 
	if (endNode <= startNode) return NULL;
	if (endNode == startNode + 1) return startNode;
	KDNode *iterateNode, *storeNode, *medianNode = startNode + (endNode - startNode) / 2;
	double pivot;
	while (true)
	{
		pivot = photon_map[medianNode->index].pos[ind];
		swap(medianNode, endNode - 1);
		storeNode = startNode;
		// Put all the nodes smaller than the median before median
		for (iterateNode = startNode; iterateNode < endNode; iterateNode++)
		{
			if (photon_map[iterateNode->index].pos[ind] < pivot)
			{
				if (iterateNode != storeNode)
					swap(iterateNode, storeNode);
				storeNode++;
			}
		}
		swap(storeNode, endNode - 1);

		if (photon_map[storeNode->index].pos[ind] == photon_map[medianNode->index].pos[ind])
			return medianNode;

		if (storeNode > medianNode)
			endNode = storeNode;
		else
			startNode = storeNode;
	}
}

KDNode* PhotonMap::MakeKDTree(KDNode *treeNode, int nodeSize, int ind, int dim)
{
	KDNode *node;
	if (!nodeSize) return NULL;
	if ((node = FindMedian(treeNode, treeNode + nodeSize, ind)))
	{
		ind = (ind + 1) % dim; // Update index
		node->left = MakeKDTree(treeNode, node - treeNode, ind, dim);
		node->right = MakeKDTree(node + 1, treeNode + nodeSize - (node + 1), ind, dim);
	}
	return node;
}

void PhotonMap::initKDtree(int k)
{
	int num_pts = photon_map.size();
	tree_mem = new KDNode[num_pts];
	for (size_t i = 0; i < num_pts; i++)
	{
		tree_mem[i].index = i;
	}
	std::cout << "Num of Points" << num_pts << std::endl;
	std::cout << "Making KDtree Now...." << std::endl;
	tree_root = MakeKDTree(tree_mem, num_pts, 0, 3);
	bestDist = new double[k];
	bestNodes = new KDNode*[k];
	for (size_t i = 0; i < k; i++)
	{
		KDNode temp_pts;
		temp_pts.index = 0;
		bestNodes[i] = &temp_pts;
		bestDist[i] = 1e6;
	}
}


inline void PhotonMap::pushNodes(KDNode **resultTree, double *dist, int dim)
{
	double temp;
	KDNode* tempNode;
	for (size_t i = dim - 1; i > 0; i--)
	{
		if (!resultTree[i - 1] || dist[i] < dist[i - 1])
		{
			tempNode = resultTree[i - 1];
			resultTree[i - 1] = resultTree[i];
			resultTree[i] = tempNode;
			temp = dist[i - 1];
			dist[i - 1] = dist[i];
			dist[i] = temp;
		}
		else
			return;
	}
}

void PhotonMap::FindKNearest(KDNode *rootNode, Point3D queryNode, int ind,
	int dim, int k, KDNode **bestNodes, double *bestDist)
{
	double distRoot, distRootInd, distSqRootInd;
	if (!rootNode) return;
	distRoot = dist(rootNode, queryNode, dim);
	distRootInd = photon_map[rootNode->index].pos[ind] - queryNode[ind];
	distSqRootInd = distRootInd * distRootInd;

	if (!bestNodes[k - 1] || distRoot < bestDist[k - 1])
	{
		bestNodes[k - 1] = rootNode;
		bestDist[k - 1] = distRoot;
		pushNodes(bestNodes, bestDist, k);
	}

	if (!bestDist[k - 1]) return; // Root is exact equal to query

	if (ind++ >= dim) ind = 0; // Recursively change index

	FindKNearest(distRootInd > 0 ? rootNode->left : rootNode->right,
		queryNode, ind, dim, k, bestNodes, bestDist);
	if (distSqRootInd >= bestDist[k - 1]) return; // No need to search other parts
	FindKNearest(distRootInd > 0 ? rootNode->right : rootNode->left,
		queryNode, ind, dim, k, bestNodes, bestDist);

}



Colour PhotonMap::findKNNColorWithKDTree(Point3D point, int num)
{ // find knn

	// printKDNode(tree_root, 4);
	for (size_t i = 0; i < num; i++)
		bestDist[i] = 1e6;
	FindKNearest(tree_root, point, 0, 3, num, bestNodes, bestDist);
	
	 // print the nearest neighbours
	// std::cout << "Query Point is" << std::endl;
	// std::cout << point << std::endl;
	// std::cout << "K Nearest Neighbour is" << std::endl;
	//for (size_t i = 0; i < num; i++)
	//{
	//	std::cout << photon_map[bestNodes[i]->index].pos << std::endl;
	//	std::cout << bestDist[i] << std::endl;
	//}

	size_t photon_index;
	Colour mean_color;
	double weight = 0;
	double weight_sum = 0;
	double knn_radius = sqrt(bestDist[num - 1]);
	for (int i = 0; i < num; i++)
	{
		photon_index = bestNodes[i]->index;
		weight = 1 - bestDist[i] / (pow(knn_radius, 2) + DBL_EPSILON);
		mean_color = mean_color + weight * photon_map[photon_index].power;
		weight_sum = weight_sum + weight;
	}
	mean_color = (1 / (knn_radius*knn_radius + 0.05)) * mean_color;
	// std::cout << mean_color << std::endl;
	return mean_color;

}