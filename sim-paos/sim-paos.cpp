#include "sim-paos.h"

// CONTROL VARIABLES
const double merge_distance = 1;
bool en_benchmark = false;

// Checks for collisions between object i and objects 0 to i - 1
void checkCollisions(std::vector<Object>& objects, size_t& i) {
    auto it = objects.begin();
    while (it != objects.begin() + i) {  // for all objects j < i
        if (dst_sqr(objects[i], *it) < sqr(merge_distance)) {
            // Collision detected, merge object j (it) into i
            objects[i].mass += it->mass;
            objects[i].vx += it->vx;
            objects[i].vy += it->vy;
            objects[i].vz += it->vz;

            // Printing (only in debug)
#ifndef NDEBUG
            std::printf("Two bodies collided. New mass: %.2E\n", objects[i].mass);
#endif

            // Delete second object
            it = objects.erase(it);
            i--;  // Decrement i, as we just deleted a entry j < i
        }
        else {
            ++it;
        }
    }
}

int main(int argc, char** argv) {
    auto t1 = std::chrono::high_resolution_clock::now();  // Start measuring the execution time

    // Check the input parameters
    const char* arguments[5] = {"num_objects", "num_iterations",
                                "random_seed", "size_enclosure", "time_step"};

    // If an optional argument is provided and it matches en_benchmark, enable benchmark mode
    if (argc == 7) {
        std::string arg6(argv[6]);
        if (arg6 == "en_benchmark") {
            en_benchmark = true;
        }
    }
    if (!en_benchmark) {
        std::cout << "sim-paos invoked with " << argc - 1 << " parameters."
            << "\n"
            << "Arguments:\n";
    }

    // Iterate for every argument needed
    if (!en_benchmark) {
        for (int i = 1; i < 6; i++) {
            // Only assign variables that exist, variables that don't exist get an ?
            if (argc > i) {
                std::cout << " " << arguments[i - 1] << ": " << argv[i] << "\n";
            }
            else {
                std::cout << " " << arguments[i - 1] << ": ?"
                    << "\n";
            }
        }
    }

    // Check the parameter count
    if (argc < 6 || argc > 7) {
        std::cerr << "Error: Wrong number of parameters\n";
        return -1;
    }

    const int num_objects = std::stoi(argv[1]);
    const int num_iterations = std::stoi(argv[2]);
    const uint64_t seed = std::stoull(argv[3]);
    const double size_enclosure = std::stod(argv[4]);
    const double time_step = std::stod(argv[5]);

    if (num_objects < 0) {
        std::cerr << "Error: Invalid number of objects\n";
        return -2;
    }
    if (num_iterations < 0) {
        std::cerr << "Error: Invalid number of iterations\n";
        return -2;
    }

    // Seed is already an unsigned 64 bit integer, no need to check for validity

    if (size_enclosure < 0) {
        std::cerr << "Error: Invalid box size\n";
        return -2;
    }
    if (time_step < 0) {
        std::cerr << "Error: Invalid time increment\n";
        return -2;
    }

    // Initialize the RNG
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<> uniform_distr(0, size_enclosure);
    std::normal_distribution<double> normal_distr(1E21, 1E15);

    // Create the necessary amount of objects and store them in a vector of class Object
    std::vector<Object> objects(num_objects);
    auto rnd_object = [&gen, &uniform_distr, &normal_distr] {
        double x = uniform_distr(gen);
        double y = uniform_distr(gen);
        double z = uniform_distr(gen);
        double m = normal_distr(gen);
        return Object(m, x, y, z);
    };
    std::generate(objects.begin(), objects.end(), rnd_object);

    // Check for collisions before starting
    for (size_t i = 0; i < objects.size(); i++) {
        checkCollisions(objects, i);
    }

    // Print the initial config
    std::ofstream initial;
    initial.open("init_config.txt");
    initial << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << objects.size() << "\n";
    for (const auto& i : objects) {
        initial << std::fixed << std::setprecision(3) << i.x << " " << i.y << " " << i.z << " " << i.vx << " " << i.vy << " " << i.vz << " " << i.mass << "\n";
    }
    initial.close();

    // Time loop
    for (size_t iteration = 0; iteration < (unsigned) num_iterations; iteration++) {  

        // Calculate the force, change in velocity and position for every object
        for (size_t i = 0; i < objects.size(); i++) {
            const size_t objectsSize = objects.size();
            for (size_t j = i + 1; j < objectsSize; j++) {
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
            }

            // All forces on objects[i] are now computed, calculate the velocity change
            // F=ma -> a=F/m
            // dv=a*dt -> dv=F/m*dt
            objects[i].vx += objects[i].fx / objects[i].mass * time_step;
            objects[i].vy += objects[i].fy / objects[i].mass * time_step;
            objects[i].vz += objects[i].fz / objects[i].mass * time_step;

            // Reset all forces to zero
            objects[i].fx = objects[i].fy = objects[i].fz = 0;

            // Update the position of the object
            objects[i].x += objects[i].vx * time_step;
            objects[i].y += objects[i].vy * time_step;
            objects[i].z += objects[i].vz * time_step;

            // Check for boundary bounce
            if (objects[i].x > size_enclosure) {
                objects[i].x = size_enclosure;
                objects[i].vx *= -1;
            }
            if (objects[i].y > size_enclosure) {
                objects[i].y = size_enclosure;
                objects[i].vy *= -1;
            }
            if (objects[i].z > size_enclosure) {
                objects[i].z = size_enclosure;
                objects[i].vz *= -1;
            }

            if (objects[i].x < 0) {
                objects[i].x = 0;
                objects[i].vx *= -1;
            }
            if (objects[i].y < 0) {
                objects[i].y = 0;
                objects[i].vy *= -1;
            }
            if (objects[i].z < 0) {
                objects[i].z = 0;
                objects[i].vz *= -1;
            }

            // Check for collisions (for all objects j < i)
            checkCollisions(objects, i);
        }

        // Printing (only in debug)
#ifndef NDEBUG
        std::printf("it %d\t  x\t\t  y\t\t  z\n", (int)iteration);
        unsigned int j = 0;
        for (const auto& i : objects) {
            std::printf("%04d: f: %.2E \t%.2E \t%.2E\n", j, i.fx, i.fy, i.fz);
            std::printf("%04d: p: %.2E \t%.2E \t%.2E\n", j, i.x, i.y, i.z);
            std::printf("%04d: v: %.2E \t%.2E \t%.2E\n\n", j, i.vx, i.vy, i.vz);
            j++;
        }
        if (objects.size() > 1)
            std::printf("Distance (0-1) %.2E\n", std::sqrt(dst_sqr(objects[0], objects[1])));
#endif
    }  // END OF TIME LOOP

    // Printing final config
    std::ofstream final;
    final.open("final_config.txt");
    final << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << objects.size() << "\n";
    for (const auto& i : objects) {
        final << std::fixed << std::setprecision(3) << i.x << " " << i.y << " " << i.z << " " << i.vx << " " << i.vy << " " << i.vz << " " << i.mass << "\n";
    }
    final.close();

    // Measure execution time and print it
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_ms = t2 - t1;
    if (en_benchmark) {
        std::printf("%f", exec_ms.count());
    }
    else {
        std::printf("Total execution time was %f ms.\n", exec_ms.count());
    }


    return 0;
}
