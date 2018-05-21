#include "Sphere.h"

#include "../util/Common.h"
#include "../util/Draw_header.h"
#include <Eigen/Dense>
#include <iostream>
Sphere::Sphere(){

}
void Sphere::init(int numSpheres, double**spheres, Tri mTriangles){

	std::cout << "\nSphere:: Initialize\n" << std::endl;
	m_spheres.reserve(numSpheres);
	//get coordinates and radius
	for (int i = 0; i < numSpheres; i++){
		sphere s;
		s.x[0] = spheres[i][0];
		s.x[1] = spheres[i][1];
		s.x[2] = spheres[i][2];
		s.x[3] = spheres[i][3];
		s.x[3] = s.x[3] * s.x[3];
		s.tag = int(spheres[i][4]);
		m_spheres.push_back(s);
	}

	//build kd tree
	m_KdTree.BuildTree(numSpheres, spheres, 3);

	if (false){
		SphereInfo(0);
	}
	//get overlapping spheres 
    m_inside = new int[numSpheres];
	for (int i = 0; i < numSpheres; i++){
		setOverlap(i);
	}

	//compute all intersection pairs 
	for (int i = 0; i < numSpheres; i++){
		setSeeds(i, true,false);
	}


	//checking sampling conditions (no sphere contain the center of other spheres)
	/*for (int i = 0; i <int(m_spheres.size()); i++){
		for (int j = 0; j < int(m_spheres[i].overlap.size()); j++){
			int sp = m_spheres[i].overlap[j];
			if (Dist(m_spheres[i].x[0], m_spheres[i].x[1], m_spheres[i].x[2], m_spheres[sp].x[0], m_spheres[sp].x[1], m_spheres[sp].x[2])
				< std::min(m_spheres[i].x[3], m_spheres[sp].x[3])){
				std::cout << "\n sampling condition vilated:: sphere [" << i << "] contains center of sphere [" << sp << "] " << std::endl;
				Draw_m_sphere(i);
				Draw_m_sphere(sp);
			}

		}
	}*/


	//get the triangle associated with each sphere in brute force way 	
	for (int i = 0; i < numSpheres; i++){
		int myT;
		double close_dd = SIZE_MAX;
		for (size_t t = 0; t < mTriangles.ids.size(); t++){
			int id0 = mTriangles.ids[t][0];
			int id1 = mTriangles.ids[t][1];
			int id2 = mTriangles.ids[t][2];
			double x2, y2, z2;
			double dd = PointTriangleDistance(mTriangles.coord[id0][0], mTriangles.coord[id0][1], mTriangles.coord[id0][2],
				                              mTriangles.coord[id1][0], mTriangles.coord[id1][1], mTriangles.coord[id1][2],
				                              mTriangles.coord[id2][0], mTriangles.coord[id2][1], mTriangles.coord[id2][2],
				                              spheres[i][0], spheres[i][1], spheres[i][2], x2, y2, z2);
			if (dd < close_dd){
				myT = t;
				close_dd = dd;
			}
		}
		m_spheres[i].tri = myT;

		
		if (false){
			DrawOneTri(myT, mTriangles.coord[mTriangles.ids[myT][0]][0], mTriangles.coord[mTriangles.ids[myT][0]][1], mTriangles.coord[mTriangles.ids[myT][0]][2],
				            mTriangles.coord[mTriangles.ids[myT][1]][0], mTriangles.coord[mTriangles.ids[myT][1]][1], mTriangles.coord[mTriangles.ids[myT][1]][2],
				            mTriangles.coord[mTriangles.ids[myT][2]][0], mTriangles.coord[mTriangles.ids[myT][2]][1], mTriangles.coord[mTriangles.ids[myT][2]][2]);
		}
	}		
}
Sphere::~Sphere(){
	delete[] m_inside;
}


void Sphere::setOverlap(int id){
	//compute the overlapping spheres with sphere id and store them in m_spheres
	if (id >= int(m_spheres.size())){
		std::cout << "Error Sphere::getSeeds:: trying to access out-of-range!!!" << std::endl;
		return;
	}
	if (m_spheres[id].overlap.size() > 0){
		m_spheres[id].overlap.clear();
	}

	//get all spheres whose centers lies inside a sphere centered at i
	//with radius 10 times the i's radius
	int  numInside = 0;
	m_KdTree.rangeQuery(m_spheres[id].x[0], m_spheres[id].x[1], m_spheres[id].x[2], 15*15 * m_spheres[id].x[3], m_inside, numInside);

	//for each sphere you get, test to check if it actually touches i
	for (int j = 0; j < numInside; j++){
		if (m_inside[j] == id){ continue; }
		if (Dist(m_spheres[m_inside[j]].x[0], m_spheres[m_inside[j]].x[1], m_spheres[m_inside[j]].x[2],
			m_spheres[id].x[0], m_spheres[id].x[1], m_spheres[id].x[2]) <
			m_spheres[m_inside[j]].x[3] + m_spheres[id].x[3] + 2.0*sqrt(m_spheres[id].x[3] * m_spheres[m_inside[j]].x[3]) + _tol_sq){

			m_spheres[id].overlap.push_back(m_inside[j]);
		}
	}
}
bool Sphere::isSliverCandidate(double xx, double yy, double zz, double r_2, int skip1, int skip2, std::vector<int>overlap){
	//check if sphere at xx,yy,zz with radius r_2 will create a sliver with 
	//spheres in overlap while discarding skip1 and skip2
	/*sphere mySphere;
	mySphere.x[0] = xx;
	mySphere.x[1] = yy;
	mySphere.x[2] = zz;
	mySphere.x[3] = r_2;
	
	m_skip.resize(4);
	m_skip[0] = skip1;
	m_skip[1] = skip2;

	for (int i = 0; i<int(overlap.size())-1; i++){
		int sp1 = overlap[i];
		
		if (sp1 == skip1 || sp1 == skip2){ continue; }
		m_skip[2] = sp1;

		for (int j = i + 1; j<int(overlap.size()); j++){
			int sp2 = overlap[j];
			if (sp2 == skip1 || sp2 == skip2){ continue; }
			seed pair1, pair2;
			m_skip[3] = sp2;
			if (insect(m_spheres[sp1], m_spheres[sp2], mySphere, pair1, pair2)){
				
				bool pair1_covered = isCovered(pair1, overlap);
				bool pair2_covered = isCovered(pair2, overlap);


				if ((pair1_covered && !pair2_covered) || (pair2_covered && !pair1_covered)){
					return true;
				}
			}

		}

	}
	m_skip.clear();
	return false;*/


	//use skip1 as a subsitiution for (xx,yy,zz) by changing its info and don't treat it 
	//as a skip sphere anymore
	//do triple loop over all spheres in overlap and detect all intersection pairs 
	//in order to account for skip1 existance, add it to overlap (the caller won't see this)
	//do the coverage check appropaiately using either overlapping spheres or this overlap vector
	//note that 
	sphere mySphere = m_spheres[skip1];
	
	m_spheres[skip1].x[0] = xx;
	m_spheres[skip1].x[1] = yy;
	m_spheres[skip1].x[2] = zz;
	m_spheres[skip1].x[3] = r_2;

	m_skip.resize(4);
	m_skip[0] = skip2;
	
	std::vector<int>skip1_vec(skip1);
	if (std::find(overlap.begin(), overlap.end(), skip1) == overlap.end()){
		
		overlap.push_back(skip1);
	}

	for (int k = 0; k<int(overlap.size()) - 2; k++){
		int sp1 = overlap[k];
		if (sp1 == skip2){ continue; }
		m_skip[1] = sp1;

		for (int i = k + 1; i<int(overlap.size()) - 1; i++){
			int sp2 = overlap[i];
			if (sp2 == skip2){ continue; }
			m_skip[2] = sp2;

			for (int j = k + 1; j<int(overlap.size()); j++){
				int sp3 = overlap[j];
				if (sp3 == skip2){ continue; }
				m_skip[3] = sp3;
				
				seed pair1, pair2;

				if (insect(m_spheres[sp1], m_spheres[sp2], m_spheres[sp3], pair1, pair2)){

					bool pair1_covered, pair2_covered;
					if (sp1 == skip1 || sp2 == skip1 || sp3 == skip1){
						pair1_covered = isCovered(pair1, overlap);
						pair2_covered = isCovered(pair2, overlap);
					}
					else{
						pair1_covered = isCovered(pair1, m_spheres[sp1].overlap) || isCovered(pair1, skip1_vec);
						pair2_covered = isCovered(pair2, m_spheres[sp1].overlap) || isCovered(pair1, skip1_vec);
					}


					if ((pair1_covered && !pair2_covered) || (pair2_covered && !pair1_covered)){
						m_spheres[skip1] = mySphere;						
						return true;
					}
				}

			}

		}
	}
	
	m_spheres[skip1] = mySphere;
	m_skip.clear();
	return false;
}
bool Sphere::addSphere(double x, double y, double z, double r_2, int tri){
	//Add a new sphere and 
	//compute its overlapping sphere and store them 
	//and its intersection pairs and compute them 
	//return true at success
	sphere s;
	s.x[0] = x;
	s.x[1] = y;
	s.x[2] = z;
	s.x[3] = r_2;
	s.tri = tri;
	m_spheres.push_back(s);

	//setOverlap(m_spheres.size() - 1);

	setSeeds(m_spheres.size() - 1, false,false);
	

	return true;
}
bool Sphere::setSeeds(int id, bool withOrder, bool checkFirst)
{
	//compute and set the intersection pairs that live over the sphere id
	//withOrder should be true if you are computing the intersection pair
	//for all the spheres (in m_spheres) by iterating over them 
	//this prevent re-computing the same intersection pair as a new pair 
	//otherwise, when compute the intersection pair of a new sphere 
	//withOrder should be false which is the default
	///if checkFirst, then will check if the three spheres have a common intersection pair before adding it
	//if so, the intersection won't be added 
	if (m_spheres[id].overlap.empty()){
		std::cout << " id=" << id << std::endl;
		std::cout << "Error Sphere::setSeeds:: should populate m_spheres[id].overlap first (call setOverlap)" << std::endl;
	}
	m_skip.resize(3);
	m_skip[0] = id;

	for (int i = 0; i <int( m_spheres[id].overlap.size()) - 1; i++){
		int id1 = m_spheres[id].overlap[i];
		m_skip[1] = id1;
		for (int j = i+1; j < int(m_spheres[id].overlap.size()); j++){
			int id2 = m_spheres[id].overlap[j];
			m_skip[2] = id2;
			if (withOrder && (id > id1 || id > id2)){
				//id should be the smallest index 
				continue;
			}

			seed pair1, pair2;
			if (insect(m_spheres[id], m_spheres[id1], m_spheres[id2], pair1, pair2)){

				bool pair1_covered = isCovered(pair1, m_spheres[id].overlap);
				bool pair2_covered = isCovered(pair2, m_spheres[id].overlap);

				if (pair1_covered && pair2_covered){ continue; }

				if ((pair1_covered && !pair2_covered) || (pair2_covered && !pair1_covered)){
					std::cout <<"Error Sphere::setSeeds:: Sliver found!!!" << std::endl;
				}

				if (checkFirst && areHavingCommonSeed(id, id1, id2)){
					continue;
				}

				//spheres that created this pair 
				pair1.spheres[0] = id;
				pair1.spheres[1] = id1;
				pair1.spheres[2] = id2;
				pair2.spheres[0] = id;
				pair2.spheres[1] = id1;
				pair2.spheres[2] = id2;

				//seeds sibling relations 
				pair2.sibling = m_seeds.size();
				pair1.sibling = m_seeds.size() + 1;


				//add the seeds to seed list (m_seeds)
				m_seeds.push_back(pair1);
				m_seeds.push_back(pair2);


				//add the seeds (id) to the spheres 
				m_spheres[id].seeds.push_back(m_seeds.size() - 1);
				m_spheres[id].seeds.push_back(m_seeds.size() - 2);

				m_spheres[id1].seeds.push_back(m_seeds.size() - 1);
				m_spheres[id1].seeds.push_back(m_seeds.size() - 2);

				m_spheres[id2].seeds.push_back(m_seeds.size() - 1);
				m_spheres[id2].seeds.push_back(m_seeds.size() - 2);			

			}
		}
	}
	m_skip.clear();
	return true;
}
bool Sphere::areHavingCommonSeed(int sp1, int sp2, int sp3){
	//check if the three sphere are having a common intersection pair 

	for (size_t i = 0; i < m_spheres[sp1].seeds.size(); i++){
		//for each seed in sp1 seed list 
		int seed1 = m_spheres[sp1].seeds[i];

		if (std::find(m_spheres[sp2].seeds.begin(), m_spheres[sp2].seeds.end(), seed1) != m_spheres[sp2].seeds.end()&& 
			std::find(m_spheres[sp3].seeds.begin(), m_spheres[sp3].seeds.end(), seed1) != m_spheres[sp3].seeds.end()){
			//if seed1 is in the seed list of sp2 and sp3 
			return true;
		}
	}
	return false;
}
void Sphere::re_setSeeds(std::vector<int>sphere_list){
	//for each sphere in sphere_list, recompute it seeds 
	//check if the seed exist first before adding it 
	for (size_t i = 0; i < sphere_list.size(); i++){		
		setSeeds(sphere_list[i], false, true);
	}
}

bool Sphere::isCovered(seed pair, std::vector<int>overlap){
	//check if this vertex is covered by other spheres 
	//and skip spheres stored in m_skip
	//used in conjunction  with Sphere::setSeeds
	if (m_skip.empty()){
		std::cout <<"Error Sphere::isCovered:: m_skip should containt three spheres" << std::endl;
	}
	
	//for (int i = 0; i < int(m_spheres[sphere_id].overlap.size()); i++){
	//	int sp = m_spheres[sphere_id].overlap[i];
	for (int i = 0; i < int(overlap.size()); i++){
			int sp = overlap[i];
		if (std::find(m_skip.begin(), m_skip.end(), sp) != m_skip.end()){ continue; }	
		double d = Dist(pair.x[0], pair.x[1], pair.x[2], m_spheres[sp].x[0], m_spheres[sp].x[1], m_spheres[sp].x[2]);
		if (d < m_spheres[sp].x[3] + _tol_sq){
			return true;
		}
	}

	return false;


}


bool Sphere::removeTwoAddOne(int id1_remove, int id2_remove, double x, double y, double z, double r_2, int tri)
{
	//remove two shperes and add new one
	//remove id1_remove and id2_remove and add a sphere centered at (x,y,z,r_2)
	//that lives on triangle tri

	//this means we replace the information of id1 with the new spheres 
	//and then remove id2 
	//this has lower cost than removing two spheres and adding one since we are using vector 

	//1) remove any seed that lives on id1_remove 
	std::vector<int> removeList(m_spheres[id1_remove].seeds);
	removeSeedList(removeList);

	if (!m_spheres[id1_remove].seeds.empty()){
		std::cout << " Error Sphere::removeTwoAddOne delete all seeds is not done correctly (1)" << std::endl;
	}

	//2) remove id1_remove from the overlap list of other spheres 
	for (int i = 0; i < int(m_spheres[id1_remove].overlap.size()); i++){
		int sp1 = m_spheres[id1_remove].overlap[i];
		if (!removeFromVector(id1_remove, m_spheres[sp1].overlap)){
			std::cout << " Error Sphere::removeTwoAddOne sphere-sphere overlap relation is invalid (1)" << std::endl;
		}		
	}
	m_spheres[id1_remove].overlap.clear();


	//3) load the new sample in id1_remove position 
	m_spheres[id1_remove].x[0] = x;
	m_spheres[id1_remove].x[1] = y;
	m_spheres[id1_remove].x[2] = z;
	m_spheres[id1_remove].x[3] = r_2;
	m_spheres[id1_remove].tri = tri;
	m_spheres[id1_remove].tag = 0;


	//brute froce overlap (can not use kd-tree as the coordinates changed)
	//make sure to also add id1_remove to overlapping spheres with it
	for (int i = 0; i < int(m_spheres.size()); i++){
		if (i == id1_remove || i == id2_remove){ continue; }
		if (Dist(m_spheres[i].x[0], m_spheres[i].x[1], m_spheres[i].x[2], x, y, z) <
			m_spheres[i].x[3] + r_2 + 2.0*sqrt(r_2 * m_spheres[i].x[3]) + _tol_sq){

			m_spheres[id1_remove].overlap.push_back(i);
			m_spheres[i].overlap.push_back(id1_remove);
		}
	}
	//setOverlap(id1_remove);
	setSeeds(id1_remove,false,false);


	//Delete id2_remove
	//4) remove any seed that lives on id2_remove 
	std::vector<int> removeList1(m_spheres[id2_remove].seeds);
	removeSeedList(removeList1);
	if (!m_spheres[id2_remove].seeds.empty()){
		std::cout << " Error Sphere::removeTwoAddOne delete all seeds is not done correctly (2)" << std::endl;
	}

	//5) remove id2_remove from the overlap list of the other spheres 
	for (int i = 0; i < int(m_spheres[id2_remove].overlap.size()); i++){
		int sp1 = m_spheres[id2_remove].overlap[i];
		if (!removeFromVector(id2_remove, m_spheres[sp1].overlap)){
			std::cout << " Error Sphere::removeTwoAddOne sphere-sphere overlap relation is invalid (2)" << std::endl;
		}
	}
	m_spheres[id2_remove].overlap.clear();

	//6)one-way swap of id2_remove with back() and then pop_back last element 
	changeMyName_sphere(m_spheres.size() - 1, id2_remove);
	m_spheres[id2_remove] = m_spheres.back();
	m_spheres.pop_back();



	return true;

}
void Sphere::changeMyName_seed(int source, int target)
{
	//loop over the sphere that made up the seed source, change its name to target
	for (int i = 0; i < 3; i++){
		int sp = m_seeds[source].spheres[i];
		for (int j = 0; j < int(m_spheres[sp].seeds.size()); j++){
			if (m_spheres[sp].seeds[j] == source){
				m_spheres[sp].seeds[j] = target;
				break;
			}
		}
	}

}
void Sphere::changeMyName_sphere(int source, int target)
{
	//loop over all seeds that live on sphere source and change its name to target
	//loop over all spheres that overlap with source and change its name to target
	for (int i = 0; i <int(m_spheres[source].seeds.size()); i++){
		int seed = m_spheres[source].seeds[i];
		for (int j = 0; j < 3; j++){
			if (m_seeds[seed].spheres[j] == source){
				m_seeds[seed].spheres[j] = target;
				break;
			}
		}
	}


	for (int i = 0; i < int(m_spheres[source].overlap.size()); i++){
		int sp = m_spheres[source].overlap[i];
		for (int j = 0; j < int(m_spheres[sp].overlap.size()); j++){
			if (m_spheres[sp].overlap[j] == source){
				m_spheres[sp].overlap[j] = target;
				break;
			}
		}
	}
}

bool Sphere::getOverlap(int id, std::vector<int>&overlap){
	//get the overlapping spheres with id 
	//return false if failed 
	if (id >= int(m_spheres.size())){
		std::cout << "Error Sphere::getOverlap:: trying to access out-of-range!!!" << std::endl;
		return false;
	}
	if (overlap.size() != 0){
		overlap.clear();
	}
	for (int i = 0; i < int(m_spheres[id].overlap.size()); i++){		
		overlap.push_back(m_spheres[id].overlap[i]);
		
	}

	return true;
}
bool Sphere::getSeeds(int id, std::vector<int>&seeds){
	//get the seeds that live on the sphere id 
	//return false if fail 
	if (id >= int(m_spheres.size())){
		std::cout << "Error Sphere::getSeeds:: trying to access out-of-range!!!" << std::endl;
		return false;
	}
	if (seeds.size() != 0){
		seeds.clear();
	}
	for (int i = 0; i < int(m_spheres[id].seeds.size()); i++){
		if (m_spheres[id].seeds[i] < 0){ continue; }

		seeds.push_back(m_spheres[id].seeds[i]);
	}
	

	return true;
}
bool Sphere::insect(sphere s1, sphere s2, sphere s3, seed&pair1, seed&pair2){
	//get the intersection of s1, s2, s3
	//store it in pair1 and pair2 

	Eigen::Vector3d s1_s = Eigen::Vector3d(s1.x[0], s1.x[1], s1.x[2]) - Eigen::Vector3d(s3.x[0], s3.x[1], s3.x[2]);
	Eigen::Vector3d s2_s = Eigen::Vector3d(s2.x[0], s2.x[1], s2.x[2]) - Eigen::Vector3d(s3.x[0], s3.x[1], s3.x[2]);

	Eigen::Vector3d v1 = (s1_s).normalized();
	Eigen::Vector3d v2 = (s2_s).normalized();
	Eigen::Vector3d n = v1.cross(v2);
	n = n.normalized();

	double len1 = (s1_s).norm();
	double len2 = (s2_s).norm();

	double x1 = (len1*len1 + s3.x[3] - s1.x[3]) / (2 * len1);
	double x2 = (len2*len2 + s3.x[3] - s2.x[3]) / (2 * len2);

	Eigen::Vector2d d;
	d << (Eigen::Vector3d(s3.x[0], s3.x[1], s3.x[2]) + x1 * v1).dot(v1),
		 (Eigen::Vector3d(s3.x[0], s3.x[1], s3.x[2]) + x2 * v2).dot(v2);

	Eigen::Matrix2d A; A << v1(0, 0), v1(1, 0), v2(0, 0), v2(1, 0);
	d = A.inverse() * d;
	Eigen::Vector3d y;
	y << d(0, 0), d(1, 0), 0;

	double b = 2 * (y - Eigen::Vector3d(s1.x[0], s1.x[1], s1.x[2])).dot(n);
	double c = (y - Eigen::Vector3d(s1.x[0], s1.x[1], s1.x[2])).squaredNorm() - s1.x[3];
	double disc = b*b - 4 * c;

	if (disc < 0){
		return false;
	}

	disc = sqrt(disc);
	double t1 = (-b + disc) / 2;
	double t2 = (-b - disc) / 2;
	Eigen::Vector3d q1 = y + t1 * n;
	Eigen::Vector3d q2 = y + t2 * n;
	pair1.x[0] = q1[0];
	pair1.x[1] = q1[1];
	pair1.x[2] = q1[2];

	pair2.x[0] = q2[0];
	pair2.x[1] = q2[1];
	pair2.x[2] = q2[2];
	
	return true;
}

bool Sphere::removeSeedList(std::vector<int>&removeList)
{
	std::sort(removeList.begin(), removeList.end());

	//remove a list of seeds
	for (int s = removeList.size() - 1; s >= 0; s--){
		int my_seed = removeList[s];

		//remove it from the spheres 
		for (int sp = 0; sp < 3; sp++){
			int my_sphere = m_seeds[my_seed].spheres[sp];

			if (!removeFromVector(my_seed, m_spheres[my_sphere].seeds)){
				std::cout << " Error Sphere::removeSeedList seed-sphere relation is invalid" << std::endl;
			}
		}

		//one-way swap with the end element and then pop back
		changeMyName_seed(m_seeds.size() - 1, my_seed);
		m_seeds[my_seed] = m_seeds.back();
		m_seeds.pop_back();

	}

	return true;
}

void Sphere::Draw_all_m_spheres()
{
	double **sp = new double*[m_spheres.size()];
	for (int i = 0; i < int(m_spheres.size()); i++){
		sp[i] = new double[4];
		sp[i][0] = m_spheres[i].x[0];
		sp[i][1] = m_spheres[i].x[1];
		sp[i][2] = m_spheres[i].x[2];
		sp[i][3] = sqrt(m_spheres[i].x[3]);
	}
	DrawManySpheres("debug_out/all_spheres.obj", m_spheres.size(), sp, true);

	for (int i = 0; i <int(m_spheres.size()); i++){
		delete[] sp[i];
	}
	delete[]sp;

}
void Sphere::Draw_all_m_seeds()
{
	double **sp = new double*[m_seeds.size()];
	double r = sqrt(m_spheres[0].x[3])/10.0;
	for (int i = 0; i < int(m_seeds.size()); i++){
		sp[i] = new double[4];
		sp[i][0] = m_seeds[i].x[0];
		sp[i][1] = m_seeds[i].x[1];
		sp[i][2] = m_seeds[i].x[2];
		sp[i][3] = r;
	}
	DrawManySpheres("debug_out/all_seeds.obj", m_seeds.size(), sp, true);

	for (int i = 0; i < int(m_seeds.size()); i++){
		delete[] sp[i];
	}
	delete[]sp;

}
void Sphere::Draw_m_sphere_overlap(int id)
{
	for (int i = 0; i < int(m_spheres[id].overlap.size()); i++){
		Draw_m_sphere(m_spheres[id].overlap[i]);
	}
}
void Sphere::Draw_m_sphere(int i)
{
	DrawOneSphere(i, m_spheres[i].x[0], m_spheres[i].x[1], m_spheres[i].x[2], sqrt(m_spheres[i].x[3]), 3);
}
void Sphere::Draw_m_seed(int i)
{
	DrawOneSphere(-i, m_seeds[i].x[0], m_seeds[i].x[1], m_seeds[i].x[2], 0.05*sqrt(m_spheres[0].x[3]), 2);
}
void Sphere::Draw_m_sphere_seeds(int id)
{
	for (int i = 0; i<int(m_spheres[id].seeds.size()); i++){
		Draw_m_seed(m_spheres[id].seeds[i]);
	}
}
void Sphere::printSpheres_csv(std::string filename){

	std::ofstream file(filename, std::ios::out);
	file.precision(30);
	file << "x coord, y coord, z coord, radius" << std::endl;
	for (int i = 0; i < int(m_spheres.size()); i++){
		file << m_spheres[i].x[0] << ", " << m_spheres[i].x[1] << ", " << m_spheres[i].x[2] << ", " << sqrt(m_spheres[i].x[3]) << std::endl;
	}
	file.close();

}
void Sphere::printSeeds_csv(std::string filename){

	std::ofstream file (filename, std::ios::out);
	file.precision(30);
	file << "x1 coord, x2 coord, x3 coord" << std::endl;
	for (int i = 0; i < int(m_seeds.size()); i++){
		file << m_seeds[i].x[0] << ", " << m_seeds[i].x[1] << ", " << m_seeds[i].x[2] << std::endl;
	}
	file.close();
}
void Sphere::vcSurface(){
	std::ofstream file("simp_surface.obj", std::ios::out);
	file.precision(30);
	for (size_t i = 0; i < m_spheres.size(); i++){
		file << "v " << m_spheres[i].x[0] << " " << m_spheres[i].x[1] << " " << m_spheres[i].x[2] << std::endl;
	}

	for (size_t s = 0; s < m_seeds.size(); s++){
		int sb = m_seeds[s].sibling;
		if (s < sb){
			//only draw the triangle once
			file << "f " << m_seeds[s].spheres[0] + 1 << "	" << m_seeds[s].spheres[1] + 1 << "	" << m_seeds[s].spheres[2] + 1 << std::endl;
		}
	}
	file.close();
	
}