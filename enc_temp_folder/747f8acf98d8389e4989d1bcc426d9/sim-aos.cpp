#include "sim-aos.h"

int main()
{
	uint64_t seed = 31728674;
	const double size_enclosure = 1000000;
	const double num_objects = 100000;
	const double time_step = 0.01;
	// Check the input parameters



	// Initialize the RNG
	std::mt19937_64 gen(seed);
	std::uniform_real_distribution<> uniform_distr(0, size_enclosure);
	std::normal_distribution<double> normal_distr(1E21, 1E15);
	
	// Create the necessary amount of objects and store them in a vector of class Object
	
	std::vector<Object> objects(num_objects);
	auto rnd_object = [&gen, &uniform_distr, &normal_distr] {
		return Object(normal_distr(gen), uniform_distr(gen), uniform_distr(gen), uniform_distr(gen));
	};
	std::generate(objects.begin(), objects.end(), rnd_object);
	
	/*std::vector<Object> objects;
	for (auto i = 0; i < num_objects; i++) {
		double x = uniform_distr(gen);
		double y = uniform_distr(gen);
		double z = uniform_distr(gen);
		objects.push_back({normal_distr(gen),x,y,z });
	}*/
	

	// Forces
	objects[0].fx = objects[0].fy = objects[0].fz = 0;	// Reset the force on the first object, others are reset in the loop
	size_t objectsSize = objects.size();
	for (unsigned int i = 0; i < objectsSize; i++) {
		for (unsigned int j = i + 1; j < objectsSize; j++) {
			if (i == 0) { 
				// Reset all forces to zero
				objects[j].fx = objects[j].fy = objects[j].fz = 0;
			}
			double massGravDist = objects[i].mass * objects[j].mass * G / dst_cube(objects[i], objects[j]);
			double fx = massGravDist * (objects[i].x - objects[j].x);
			double fy = massGravDist * (objects[i].y - objects[j].y);
			double fz = massGravDist * (objects[i].z - objects[j].z);
			objects[i].fx += fx; 
			objects[j].fx -= fx;
			objects[i].fy += fy;
			objects[j].fy -= fy;
			objects[i].fz += fz;
			objects[j].fz -= fz;
		}
		// All forces on objects[i] are now computed
		// Calculate the velocity change

		// F=ma -> a=F/m
		// dv=a*dt -> dv=F/m*dt

		objects[i].vx += objects[i].fx / objects[i].mass * time_step;
		objects[i].vy += objects[i].fy / objects[i].mass * time_step;
		objects[i].vz += objects[i].fz / objects[i].mass * time_step;

		// Update the position of the object
		objects[i].x += objects[i].vx * time_step;
		objects[i].y += objects[i].vy * time_step;
		objects[i].z += objects[i].vz * time_step;
	}
	
	// Printing
	std::printf("\t  x\t\t  y\t\t  z\n");
	for (const auto& i : objects) {
		static unsigned int j = 0;
		std::printf("%04d: f: %.2E \t%.2E \t%.2E\n", j, i.fx, i.fy, i.fz);
		std::printf("%04d: p: %.2E \t%.2E \t%.2E\n",j, i.x, i.y, i.z);
		std::printf("%04d: v: %.2E \t%.2E \t%.2E\n",j, i.vx, i.vy, i.vz);
		j++;
	}
	
	return 0;
}
