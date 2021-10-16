#pragma once

#include <cmath>
#include <vector>

#define sqr(a) (a)*(a)
#define cube(a) (a)*(a)*(a)

struct Object {

	int size;
	int iterator;

	std::vector <double> mass;

	std::vector <double> x;
	std::vector <double> y;
	std::vector <double> z;

	std::vector <double> vx;
	std::vector <double> vy;
	std::vector <double> vz;

	std::vector <double> fx;
	std::vector <double> fy;
	std::vector <double> fz;
};

void init_Object(Object& n, const int num_objects){
	n.size = num_objects;
	n.iterator = 0;
}

void add_Object(Object& n, double mass, double x, double y, double z){
	n.mass.push_back(mass);
	
	n.x.push_back(x);	
	n.y.push_back(y);
	n.z.push_back(z);

	n.vx.push_back(0);
	n.vy.push_back(0);
	n.vz.push_back(0);

	n.fx.push_back(0);
	n.fy.push_back(0);
	n.fz.push_back(0);
}

inline void reset_forces(Object& n){
	std::fill(n.fx.begin(), n.fx.end(), 0);
	std::fill(n.fy.begin(), n.fy.end(), 0);
        std::fill(n.fz.begin(), n.fz.end(), 0);
}

void adjust_for_boundary(Object& n, double size_enclosure){
	for(int i = 0; i < n.size; i++){
		if(n.x[i] < 0){
			n.x[i] = 0;
			n.vx[i] *= -1;
		}
		if(n.x[i] > size_enclosure){
			n.x[i] = size_enclosure;
			n.vx[i] *= -1;
		}
		if(n.y[i] < 0){
			n.y[i] = 0;
			n.vy[i] *= -1;
		}
		if(n.y[i] > size_enclosure){
                        n.y[i] = size_enclosure;
                        n.vy[i] *= -1;
		}
		if(n.z[i] < 0){
			n.z[i] = 0;
			n.vz[i] *= -1;
		}
		if(n.z[i] > size_enclosure){
			n.z[i] = size_enclosure;
			n.vz[i] *= -1;
		}
	}
}

// j merges into i (j deleted)
void merge_objects(Object& n, int i, int j){
	// Merge attributes
	n.mass[i] += n.mass[j];
	n.vx[i] += n.vx[j];
    n.vy[i] += n.vy[j];
    n.vz[i] += n.vz[j];
	
	// Delete second object

	n.mass.erase(n.mass.begin()+j);
	n.x.erase(n.x.begin()+j);
	n.y.erase(n.y.begin()+j);
	n.z.erase(n.z.begin()+j);
	n.vx.erase(n.vx.begin()+j);
	n.vy.erase(n.vy.begin()+j);
	n.vz.erase(n.vz.begin()+j);
	n.fx.erase(n.fx.begin()+j);
	n.fy.erase(n.fy.begin()+j);
	n.fz.erase(n.fz.begin()+j);

	n.size--;
}



inline double dst_sqr(Object& n, int i1, int i2) {
	return sqr(n.x[i1] - n.x[i2]) + sqr(n.y[i1] - n.y[i2]) + sqr(n.z[i1] - n.z[i2]);
}

inline double dst_cube(Object& n, int i1, int i2) {
	return std::pow(dst_sqr(n, i1, i2), 1.5);
}

