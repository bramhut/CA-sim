#pragma once

#include <cmath>

#define sqr(a) (a)*(a)
#define cube(a) (a)*(a)*(a)

// Struct containing all relevant information about a dimension

class Object {
public:
	double mass;

	// x,y,z
	double p[3];
	double f[3] = { 0,0,0 };
	double v[3] = { 0,0,0 };

	// Set to true to remove this object (used in collision checking)
	bool removeFlag = false;
	
	Object(const double mass = 0, const double x = 0, const double y = 0, const double z = 0) : mass(mass), p{ x,y,z } {}
};

inline double dst_sqr(const Object &a, const Object &b) {
	return sqr(a.p[0] - b.p[0]) + sqr(a.p[1] - b.p[1]) + sqr(a.p[2] - b.p[2]);
}

inline double dst_cube(const Object& a, const Object& b) {
	double dst = std::sqrt(dst_sqr(a, b));
	return cube(dst);
}

