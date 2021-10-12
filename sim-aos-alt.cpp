#include "sim-aos.h"

void recurse_objects(std::vector<Object> &objects, int i, int j);
void set_to_zero(std::vector<Object> &objects, int i);
void run_timeloop(std::vector<Object> &objects);


int main(int argc, char** argv) {
	auto t1 = std::chrono::high_resolution_clock::now();	// Start measuring the execution time
	// CONTROL VARIABLES
	// uint64_t seed = 31728674;
	// const double size_enclosure = 1000000;
	// const double num_objects = 2;
	// const double time_step = 0.01;

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

    double merge_distance = 10;

	const int num_objects = std::stoi(argv[1]);
	const int num_iterations = std::stoi(argv[2]);
	const uint64_t seed = std::stoull(argv[3]);
	double size_enclosure = std::stod(argv[4]);
	double time_step = std::stod(argv[5]);

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
	
	// Create the necessary amount of objects and store them in a vector of class Object
	std::vector<Object> objects(num_objects);
	auto rnd_object = [&gen, &uniform_distr, &normal_distr] {
		return Object(normal_distr(gen), uniform_distr(gen), uniform_distr(gen), uniform_distr(gen));
	};
	std::generate(objects.begin(), objects.end(), rnd_object);

    // Print initial position of first object
    std::printf("Initial position of the first is: %f, %f, %f\n", objects[0].x, objects[0].y, objects[0].z);


    // Run the timeloop the designated number of times
    for(int i = 0; i < num_iterations; i ++){
	    run_timeloop(objects);
    }

	// Measure execution time and print it
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> exec_ms = t2 - t1;
    std::printf("Final position of the first is: %f, %f, %f\n", objects[0].x, objects[0].y, objects[0].z);
	std::printf("Total execution time was %.2f ms.\n", exec_ms.count());
	
	return 0;
}


// Updates each object's position, detects and resolves collisions
void recurse_objects(std::vector<Object> &objects, int i, int j){

    // Calculating forces
    double massGravDist = objects[i].mass * objects[j].mass * G / dst_cube(objects[i], objects[j]);
	double fx = massGravDist * (objects[j].x - objects[i].x);
	double fy = massGravDist * (objects[j].y - objects[i].y);
	double fz = massGravDist * (objects[j].z - objects[i].z);
	objects[i].fx += fx;
	objects[j].fx -= fx;
	objects[i].fy += fy;
	objects[j].fy -= fy;
	objects[i].fz += fz;
	objects[j].fz -= fz;
    
    // Recursing
    if(j+1 < objects.size()){
        recurse_objects(objects, i, j+1);
    }else if(i+2 < objects.size()){
        recurse_objects(objects, i+1, i+2);
    }
    

    // Calculating velocity for i and j
    for(int k = j; k == i || k == j; k += (i-j)){// k = i, j
        if(k < objects.size() && objects[k].flag == false){// Reduces unnecessary calculations
            // v = a*t = F/m*t
            
            objects[k].vx += objects[k].fx / objects[k].mass * time_step;
	        objects[k].vy += objects[k].fy / objects[k].mass * time_step;
	        objects[k].vz += objects[k].fz / objects[k].mass * time_step;
    
            // Calculating position, boundary control
            objects[k].x += objects[k].vx * time_step;
            if (objects[k].x < 0) { objects[k].x = 0; objects[k].vx *= -1; }
            if (objects[k].x > size_enclosure) { objects[k].x = size_enclosure; objects[k].vx *= -1; }	    
        
            objects[k].y += objects[k].vy * time_step;
            if (objects[k].y < 0) { objects[k].y = 0; objects[k].vy *= -1; }
            if (objects[k].y > size_enclosure) { objects[k].y = size_enclosure; objects[k].vy *= -1; }

            objects[k].z += objects[k].vz * time_step;
            if (objects[k].z < 0) { objects[k].z = 0; objects[k].vz *= -1; }
            if (objects[k].z > size_enclosure) { objects[k].z = size_enclosure; objects[k].vz *= -1; }

            objects[k].flag = true;
        }
    }

    // Detecting collision
    if (i < objects.size() && j < objects.size() && dst_sqr(objects[i], objects[j]) < sqr(merge_distance)){
        // Collision detected merge object j into i
		objects[j].mass += objects[i].mass;
		objects[j].vx += objects[i].vx;
		objects[j].vy += objects[i].vy;
		objects[j].vz += objects[i].vz;
		//std::printf("Two bodies collided. New mass: %.2E\n", objects[i].mass);

		// Delete second object
		objects.erase(objects.begin() + i);    
    }

}

// Sets all objects' forces and velocities to 0
void set_to_zero(std::vector<Object> &objects, int i){
    // Setting forces and velocities to 0
    objects[i].fx = 0; objects[i].fy = 0; objects[i].fz = 0;
    objects[i].flag = false;

    // Recursing
    if(i+1 < objects.size()){
        set_to_zero(objects, i+1);
    }
}


// Runs a time loop
void run_timeloop(std::vector<Object> &objects){
    if(objects.size() > 1){
        set_to_zero(objects, 0);
        recurse_objects(objects, 0, 1);
    }
}

