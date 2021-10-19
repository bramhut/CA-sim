#pragma once

#include <cmath>
#include <vector>

#define sqr(a) (a)*(a)
#define cube(a) (a)*(a)*(a)

struct Object;
inline double dst_sqr(Object* n, size_t i1, size_t i2);
inline double dst_cube(Object* n, size_t i1, size_t i2);

struct Object {

	size_t size;

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

	// Constructor
	Object(const size_t size, const uint64_t seed, const double size_enclosure) : size(size),
		mass(size),
		x(size),
		y(size),
		z(size),
		vx(size),
		vy(size),
		vz(size),
		fx(size),
		fy(size),
		fz(size)
	{
		// Initialize the RNG
		std::mt19937_64 gen(seed);
		std::uniform_real_distribution<> uniform_distr(0, size_enclosure);
		std::normal_distribution<double> normal_distr(1E21, 1E15);

		// Add the required amount of objects
		for (size_t i = 0; i < size; ++i) {
			x[i] = uniform_distr(gen);
			y[i] = uniform_distr(gen);
			z[i] = uniform_distr(gen);
			mass[i] = normal_distr(gen);
		}
	}

	// Reset the forces to zero
	inline void reset_forces() {
		std::fill(fx.begin(), fx.end(), 0);
		std::fill(fy.begin(), fy.end(), 0);
		std::fill(fz.begin(), fz.end(), 0);
	}

	inline void adjust_for_boundary(const double size_enclosure) {
		for (size_t i = 0; i < size; i++) {
			if (x[i] < 0) {
				x[i] = 0;
				vx[i] *= -1;
			}
			if (x[i] > size_enclosure) {
				x[i] = size_enclosure;
				vx[i] *= -1;
			}
			if (y[i] < 0) {
				y[i] = 0;
				vy[i] *= -1;
			}
			if (y[i] > size_enclosure) {
				y[i] = size_enclosure;
				vy[i] *= -1;
			}
			if (z[i] < 0) {
				z[i] = 0;
				vz[i] *= -1;
			}
			if (z[i] > size_enclosure) {
				z[i] = size_enclosure;
				vz[i] *= -1;
			}
		}
	}

	// j merges into i (j deleted)
	inline void merge_objects(size_t i, size_t j) {
		// Merge attributes
		mass[i] += mass[j];
		vx[i] += vx[j];
		vy[i] += vy[j];
		vz[i] += vz[j];

		// Delete second object

		mass.erase(mass.begin() + j);
		x.erase(x.begin() + j);
		y.erase(y.begin() + j);
		z.erase(z.begin() + j);
		vx.erase(vx.begin() + j);
		vy.erase(vy.begin() + j);
		vz.erase(vz.begin() + j);
		fx.erase(fx.begin() + j);
		fy.erase(fy.begin() + j);
		fz.erase(fz.begin() + j);

		size--;
	}

	inline void check_collisions(size_t& i) {
		size_t iterator = 0;
		while (iterator < i)
		{
			// Merge when distance is less than 1
			if (dst_sqr(this, i, iterator) < 1)
			{
				// Collision detected, merge iterator object into i
				merge_objects(i, iterator);

				// Decrement i, as i is now one index lower
				i--;
#ifndef NDEBUG
				std::printf("Two bodies collided. New size: %.2E\n", object.mass[i]);
#endif
			}
			else
			{
				iterator++;
			}
		}
	}

};

inline double dst_sqr(Object* n, size_t i1, size_t i2) {
	return sqr(n->x[i1] - n->x[i2]) + sqr(n->y[i1] - n->y[i2]) + sqr(n->z[i1] - n->z[i2]);
}

inline double dst_cube(Object* n, size_t i1, size_t i2) {
	return std::pow(dst_sqr(n, i1, i2), 1.5);
}
