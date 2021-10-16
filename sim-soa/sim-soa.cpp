#include "sim-soa.h"

int main(int argc, char** argv) {
	auto t1 = std::chrono::high_resolution_clock::now();	// Start measuring the execution time
	// CONTROL VARIABLES
	// uint64_t seed = 31728674;
	// const double size_enclosure = 1000000;
	// const double num_objects = 2;
	// const double time_step = 0.01;
	const double merge_distance = 500;

	// Check the input parameters
	const char *arguments[5] = {"num_objects", "num_iterations",
								"random_seed", "size_enclosure", "time_step"};
	std::cout << "sim-soa invoked with " << argc - 1 << " parameters."
			  << "\n"
			  << "Arguments:\n";

	// Iterate for every argument needed
	for (int i = 1; i < 6; i++) {
		// Only assign variables that exist, variables that don't exist get an ?
		if (argc > i) {
			std::cout << " " << arguments[i - 1] << ": " << argv[i] << "\n";
		} else {
			std::cout << " " << arguments[i - 1] << ": ?"
					  << "\n";
		}
	}

	if (argc != 6) {
		std::cerr << "Error: Wrong number of parameters";
		return -1;
	}

	const int num_objects = std::stoi(argv[1]);
	const int num_iterations = std::stoi(argv[2]);
	const uint64_t seed = std::stoull(argv[3]);
	const double size_enclosure = std::stod(argv[4]);
	const double time_step = std::stod(argv[5]);

	const unsigned int num_iterations_unsigned = (unsigned int)num_iterations;

	if (num_objects < 0) {
		std::cerr << "Error: Invalid number of objects";
		return -1;
	}
	if (num_iterations < 0) {
		std::cerr << "Error: Invalid number of iterations";
		return -1;
	}
	if ((int)seed < 0) {
		std::cerr << "Error: Invalid seed";
		return -1;
	}
	if (size_enclosure < 0) {
		std::cerr << "Error: Invalid box size";
		return -1;
	}
	if (time_step < 0) {
		std::cerr << "Error: Invalid time increment";
		return -1;
	}

	// Initialize the RNG
	std::mt19937_64 gen(seed);
	std::uniform_real_distribution<> uniform_distr(0, size_enclosure);
	std::normal_distribution<double> normal_distr(1E21, 1E15);

	// Generate Object; assign random values to mass and position
	struct Object object;

	init_Object(object, num_objects);
	
	for(int i = 0; i < num_objects; i++){
		add_Object(object, normal_distr(gen), uniform_distr(gen), uniform_distr(gen), uniform_distr(gen));
	}

	std::ofstream initial;
	initial.open("init_config.txt");
	initial << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << "\n";

	std::ofstream final;
	final.open("final_config.txt");
	final << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << "\n";



	// Time loop
	for (size_t iteration = 0; iteration < num_iterations_unsigned; iteration++) {

		// Reset all forces to zero
		reset_forces(object);

		// Calculate the force, change in velocity and position
		for (int i = 0; i < object.size; i++) {

			for(int j = i+1; j < object.size; j++) {
			
				double massGravDist = object.mass[i] * object.mass[j] * G / dst_cube(object, i, j);
				double fx = massGravDist * (object.x[j] - object.x[i]);
				double fy = massGravDist * (object.y[j] - object.y[i]);
				double fz = massGravDist * (object.z[j] - object.z[i]);
			
				object.fx[i] += fx;
				object.fx[j] -= fx;
				object.fy[i] += fy;
				object.fy[j] -= fy;
				object.fz[i] += fz;
				object.fz[j] -= fz;
			}

			// All forces on objects[i] are now computed, calculate the velocity change
			// F=ma -> a=F/m
			// dv=a*dt -> dv=F/m*dt
			
			object.vx[i] += object.fx[i] / object.mass[i] * time_step;
			object.vy[i] += object.fy[i] / object.mass[i] * time_step;
			object.vz[i] += object.fz[i] / object.mass[i] * time_step;
			
			// Update the position of the object
			
			object.x[i] += object.vx[i] * time_step;
			object.y[i] += object.vy[i] * time_step;
			object.z[i] += object.vz[i] * time_step;

			// If objects are outside of boundary, set them to the perimeter
			
			adjust_for_boundary(object, size_enclosure);

			// Check for collisions (for all objects j < i)
			
			int iterator = 0;
			while(iterator < i){
				if (dst_sqr(object, i, iterator) < sqr(merge_distance)){
					// Collision detected, merge iterator object into i

					merge_objects(object, i, iterator);
					
					// Decrement i, as i is now one index lower
					i--;

					std::printf("Two bodies collided. New size: %.2E\n", object.mass[i]);

				}else{
					iterator++;
				}
			}
		}


		if (iteration == 0) {
            // Store all initial values in file
			
			for(int i = 0; i < object.size; i ++){
				initial << std::fixed << std::setprecision(3) << object.x[i] << " " << object.y[i] << " " << object.z[i] << " " << object.vx[i] << " " << object.vy[i] << " " << object.vz[i] << " " << object.mass[i] << "\n";
			}

        } else if (iteration == num_iterations_unsigned - 1) {
			// Store all final values in file

			for(int i = 0; i < object.size; i ++){
                final << std::fixed << std::setprecision(3) << object.x[i] << " " << object.y[i] << " " << object.z[i] << " " << object.vx[i] << " " << object.vy[i] << " " << object.vz[i] << " " << object.mass[i] << "\n";
    		}
		
		}

		#ifndef NDEBUG
            std::printf("it %d\t  x\t\t  y\t\t  z\n", (int)iteration);
    	    unsigned int j = 0;
            for (int i = 0; i < object.size; i++) {
                std::printf("%04d: f: %.2E \t%.2E \t%.2E\n", j, object.fx[i], object.fy[i], object.fz[i]);
                std::printf("%04d: p: %.2E \t%.2E \t%.2E\n", j, object.x[i], object.y[i], object.z[i]);
                std::printf("%04d: v: %.2E \t%.2E \t%.2E\n\n", j, object.vx[i], object.vy[i], object.vz[i]);
                j++;
            }
            if (object.size > 1){
                std::printf("Distance (0-1) %.2E\n", std::sqrt(dst_sqr(object, 0, 1)));
			}
		#endif



	
	}// End time loop

	std::printf("Position of final object: %.2E, %.2E, %.2E", object.x[0], object.y[0], object.z[0]);

	// Measure execution time and print it
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> exec_ms = t2 - t1;
	std::printf("Total execution time was %.2f ms.\n", exec_ms.count());
	
	return 0;
}
