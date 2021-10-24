#include "sim-soa.h"

bool en_benchmark = false;

int main(int argc, char** argv) {
    auto t1 = std::chrono::high_resolution_clock::now();  // Start measuring the execution time

    // Check the input parameters
    const char* arguments[5] = {"num_objects", "num_iterations",
                                "random_seed", "size_enclosure", "time_step"};
    if (argc == 7) {
        std::string arg6(argv[6]);
        if (arg6 == "en_benchmark") {
            en_benchmark = true;
        }
    }
    if (!en_benchmark) {
        std::cout << "sim-soa invoked with " << argc - 1 << " parameters."
                  << "\n"
                  << "Arguments:\n";
    }

    // Iterate for every argument needed
    if (!en_benchmark) {
        for (int i = 1; i < 6; i++) {
            // Only assign variables that exist, variables that don't exist get an ?
            if (argc > i) {
                std::cout << " " << arguments[i - 1] << ": " << argv[i] << "\n";
            } else {
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

    // Generate Object; assign random values to mass and position
    Object object((size_t)num_objects, seed, size_enclosure);

    // Check for collisions before starting
    for (size_t i = 0; i < object.size; i++) {
        object.check_collisions(i);
    }

    // Print the initial config
    std::ofstream initial;
    initial.open("init_config.txt");
    initial << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << object.size << "\n";
    for (size_t i = 0; i < object.size; ++i) {
        initial << std::fixed << std::setprecision(3) << object.x[i] << " " << object.y[i] << " " << object.z[i] << " " << object.vx[i] << " " << object.vy[i] << " " << object.vz[i] << " " << object.mass[i] << "\n";
    }
    initial.close();

    // Time loop
    for (size_t iteration = 0; iteration < (unsigned)num_iterations; iteration++) {
        // Reset all forces to zero
        object.reset_forces();

        // Calculate the force, change in velocity and position
        for (size_t i = 0; i < object.size; i++) {
            for (size_t j = i + 1; j < object.size; j++) {
                double massGravDist = object.mass[i] * object.mass[j] * G / dst_cube(&object, i, j);
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

            object.adjust_for_boundary(size_enclosure, i);

            // Check for collisions (for all objects j < i)

            object.check_collisions(i);
        }

        // Printing (only in debug)
#ifndef NDEBUG
        std::printf("it %d\t  x\t\t  y\t\t  z\n", (int)iteration);
        unsigned int j = 0;
        for (size_t i = 0; i < object.size; i++) {
            std::printf("%04d: f: %.2E \t%.2E \t%.2E\n", j, object.fx[i], object.fy[i], object.fz[i]);
            std::printf("%04d: p: %.2E \t%.2E \t%.2E\n", j, object.x[i], object.y[i], object.z[i]);
            std::printf("%04d: v: %.2E \t%.2E \t%.2E\n\n", j, object.vx[i], object.vy[i], object.vz[i]);
            j++;
        }
        if (object.size > 1) {
            std::printf("Distance (0-1) %.2E\n", std::sqrt(dst_sqr(&object, 0, 1)));
        }
#endif

    }  // End time loop

    // Printing final config
    std::ofstream final;
    final.open("final_config.txt");
    final << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << object.size << "\n";
    for (size_t i = 0; i < object.size; ++i) {
        final << std::fixed << std::setprecision(3) << object.x[i] << " " << object.y[i] << " " << object.z[i] << " " << object.vx[i] << " " << object.vy[i] << " " << object.vz[i] << " " << object.mass[i] << "\n";
    }
    final.close();

    // Measure execution time and print it
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> exec_ms = t2 - t1;
    if (en_benchmark) {
        std::printf("%f", exec_ms.count());
    } else {
        std::printf("Total execution time was %f ms.\n", exec_ms.count());
    }

    return 0;
}
