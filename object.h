#pragma once

#include <cmath>

#define sqr(a) (a)*(a)
#define cube(a) (a)*(a)*(a)

struct Object {
	double mass;

	double x;
	double y;
	double z;

	double vx = 0;
	double vy = 0;
	double vz = 0;

	double fx = 0;
	double fy = 0;
	double fz = 0;

	Object(const double mass = 0, const double x = 0, const double y = 0, const double z = 0) : mass(mass), x(x), y(y), z(z) {}
};

inline double dst_sqr(Object &a, Object &b) {
	return sqr(a.x - b.x) + sqr(a.y - b.y) + sqr(a.z - b.z);
}

inline double dst_cube(Object& a, Object& b) {
	return std::pow(dst_sqr(a,b), 1.5);
}

