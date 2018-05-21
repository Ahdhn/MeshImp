#ifndef __SPHERE__
#define __SPHERE__
#include <vector>
#include "../util/KdTree.h"
#include <iostream>

class Sphere
{
public:
	Sphere();	
	~Sphere();

	//initialize the data structure
	void init(int numSpheres, double**spheres, Tri mTriangles);

	//Add a new sphere and 
	//compute its overlapping sphere and store them 
	//and its intersection pairs and compute them 
	//return true at success
	bool addSphere(double x, double y, double z, double r_2, int tri);

	//remove the sphere by removing its intersection pairs
	//delete it from all the sphere that overlap with it 
	//returen true at success 
	//bool removeSphere(int id);


	//remove two shperes and add new one
	//remove id1_remove and id2_remove and add a sphere centered at (x,y,z,r_2)
	//that lives on triangle tri
	bool removeTwoAddOne(int id1_remove, int id2_remove, double x, double y, double z, double r_2, int tri);
	

	
	//get the coordinates of sphere id
	//return false if it does not exist
	bool getCoordinates(int id, double&x, double&y, double&z, double&r_2, int&tri){
		if (id >= int(m_spheres.size()) || id<0){ return false; }
		x = m_spheres[id].x[0];
		y = m_spheres[id].x[1];
		z = m_spheres[id].x[2];
		r_2 = m_spheres[id].x[3];
		tri = m_spheres[id].tri;		
		return true;
	};

	int getTag(int id){
		return m_spheres[id].tag;
	}

	//get the coordinates of sphere id
	//return false if it does not exist
	bool getSeedCoordinates(int id, double&x, double&y, double&z){
		x = m_seeds[id].x[0];
		y = m_seeds[id].x[1];
		z = m_seeds[id].x[2];
		return true;
	
	}
	//get the overlapping spheres with id 
	//return false if failed 
	bool getOverlap(int id, std::vector<int>&overlap);

	//get the seeds that live on the sphere id 
	//return false if fail 
	bool getSeeds(int id, std::vector<int>&seeds);

	//return total number of spheres stored currently
	int getNumSphere(){ return m_spheres.size(); }

	void printSpheres_csv(std::string filename);
	void printSeeds_csv(std::string filename);


	bool isSliverCandidate(double xx, double yy, double zz, double r_2, int skip1, int skip2, std::vector<int>overlap);

	//for each sphere in sphere_list, recompute it seeds 
	//check if the seed exist first before adding it 
	void re_setSeeds(std::vector<int>sphere_list);


	//output surface
	void vcSurface();

	//display info
	void SphereInfo(int id){
		std::cout << "\n Sphere[" << id << "] at (" << m_spheres[id].x[0] << ", " << m_spheres[id].x[1] << ", " << m_spheres[id].x[2] << ")" << std::endl;
		std::cout << "       radius = " << m_spheres[id].x[3] << " on triangle[" << m_spheres[id].tri << "] with tag[" << m_spheres[id].tag << "]" << std::endl;
		std::cout << "       overlap with(";
		for (size_t i = 0; i < m_spheres[id].overlap.size(); i++){
			std::cout << m_spheres[id].overlap[i];
			if (i == m_spheres[id].overlap.size() - 1){
				std::cout << ")";
			}
			else{
				std::cout << ", ";
			}
		}
		std::cout << std::endl;
	}


	//drawing routines 
	void Draw_m_sphere(int i);
	void Draw_m_sphere_seeds(int i);
	void Draw_m_seed(int i);
	void Draw_m_sphere_overlap(int i);
	void Draw_all_m_spheres();
	void Draw_all_m_seeds();

private:

	//Seed here means intersection points 
	/************************************/

	KdTree m_KdTree; 

	struct seed {
		double x[3]; //x,y,z
		int sibling; //its sibling (seeds are created in pairs)
		int spheres[3];//the three spheres created this seed 		
	};
	std::vector<seed> m_seeds;

	struct sphere{
		double x[4];//x,y,z,r
		int tri;//the trie where the center of the sphere lives 
		std::vector<int> overlap;//the spheres that this sphere overlap with 
		std::vector<int> seeds;//the seeds this sphere produced
		int tag;
	};
	std::vector<sphere> m_spheres;
	int*m_inside;

	void setOverlap(int id);
	
	std::vector<int>m_skip;

	//get the intersection of s1, s2, s3
	//store it in pair1 and pair2 
	bool insect(sphere s1, sphere s2, sphere s3, seed&pair1, seed&pair2);

	//compute and set the intersection pairs that live over the sphere id
	//withOrder should be true if you are computing the intersection pair
	//for all the spheres (in m_spheres) by iterating over them 
	//this prevent re-computing the same intersection pair as a new pair 
	//otherwise, when compute the intersection pair of a new sphere 
	//withOrder should be false which is the default
	bool setSeeds(int id, bool withOrder, bool checkFirst);
	
	//remove a list of seeds 
	bool removeSeedList(std::vector<int>&removeList);

	
	bool isCovered(seed pair1, std::vector<int>overlap);

	//check if the three sphere are having a common intersection pair 
	bool areHavingCommonSeed(int sp1, int sp2, int sp3);

	

	void changeMyName_seed(int source, int target);
	void changeMyName_sphere(int source, int target);
	
	
};



#endif /*__SPHERE__*/