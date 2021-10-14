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
	struct Object object = new struct Object();
	init_Object(object, num_objects);
	
	for(int i = 0; i < num_objects; i++){
		add_Object(n, normal_distr(gen), uniform_distr(gen), uniform_distr(gen), uniform_distr(gen));
	}


	// Time loop
		// Reset all forces to zero
		//
		// note: use function "reset_forces" from object.h

		// Calculate the force, change in velocity and position

			// All forces on objects[i] are now computed, calculate the velocity change
			// F=ma -> a=F/m
			// dv=a*dt -> dv=F/m*dt

			// Update the position of the object

			// Check for boundary bounce
			//
			// note: use function "adjust_for_boundary" from object.h		


			// Check for collisions (for all objects j < i)




	// Measure execution time and print it
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> exec_ms = t2 - t1;
	std::printf("Total execution time was %.2f ms.", exec_ms.count());
	
	return 0;
}
