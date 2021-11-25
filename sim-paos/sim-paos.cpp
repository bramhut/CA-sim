#include "sim-paos.h"

/*


Include pseudocode in the report!


*/

// GLOBAL VARIABLES

// SETTINGS
int num_objects;
int num_iterations;
uint64_t seed;
double size_enclosure;
double time_step;
bool en_benchmark = false;

// OBJECTS VECTOR
std::vector<Object> objects;

// FUNCTIONS

// Watch class used for easy benchmarking
watch collisionWatch, updateObjWatch, totalWatch;

// Compare the pairs as if they were executed in sequential order (i first then j)
inline bool operator<(const Pair& p1, const Pair& p2) {
    return (p1.j - p1.i * num_objects) < (p2.j - p2.i * num_objects);
}

// Checks for collisions between object i and objects 0 to i - 1
void checkCollisions() {
    collisionWatch.start();

    std::set<Pair> toRemove;
    int objectsSize = (int) objects.size();

#pragma omp parallel for schedule(guided)
    for (int i = 0; i < objectsSize; i++) {
        for (int j = i - 1; j >= 0; j--) {
            if (dst_sqr(objects[i], objects[j]) < 1) {
#pragma omp critical
                toRemove.emplace(i, j);
            }
        }
    }

    // All collisions have been detected, now merge the collided objects
    //std::printf("New collision check\n");
    bool needRemoval = !toRemove.empty();
    while (!toRemove.empty()) {
        // Remove the last element from the set
        auto it = std::prev(toRemove.end());
        auto i = (*it).i;
        auto j = (*it).j;
        toRemove.erase(it);

        if (objects[j].removeFlag) {
            continue;
        }
        objects[j].removeFlag = true;

        objects[i].mass += objects[j].mass;
        objects[i].v[0] += objects[j].v[0];
        objects[i].v[1] += objects[j].v[1];
        objects[i].v[2] += objects[j].v[2];

        //std::printf("\tMerging %zi -> %zi\n", j,i);
    }

    // Remove elements if necessary
    if (needRemoval) {
        objects.erase(
            std::remove_if(
                objects.begin(),
                objects.end(),
                [&](const Object obj) -> bool {
                    return obj.removeFlag;
                }),
            objects.end());
    }
    collisionWatch.stop();
}

void updateObjects() {
    updateObjWatch.start();
    auto objectsSize = objects.size();
    //std::printf("Updating dim %i, objects size %zi\n", dim, objectsSize);
    for (size_t i = 0; i < objectsSize; ++i) {
        for (size_t j = i + 1; j < objectsSize; j++) {
            double mgd = objects[i].mass * objects[j].mass * G / dst_cube(objects[i], objects[j]);
            for (size_t dim = 0; dim < 3; dim++) {
                double f = mgd * (objects[j].p[dim] - objects[i].p[dim]);
                objects[i].f[dim] += f;
                objects[j].f[dim] -= f;
            }
        }
        // For all dimensions
        for (size_t dim = 0; dim < 3; dim++) {
            // Calculate velocity
            objects[i].v[dim] += objects[i].f[dim] / objects[i].mass * time_step;

            // Reset the force to zero
            objects[i].f[dim] = 0;

            // Update the position of the object
            objects[i].p[dim] += objects[i].v[dim] * time_step;

            // Check for boundary bounce
            if (objects[i].p[dim] > size_enclosure) {
                objects[i].p[dim] = size_enclosure;
                objects[i].v[dim] *= -1;
            }

            if (objects[i].p[dim] < 0) {
                objects[i].p[dim] = 0;
                objects[i].v[dim] *= -1;
            }
        }
    }
    updateObjWatch.stop();
}

int main(int argc, char** argv) {
    totalWatch.start();  // Start measuring the execution time

    //omp_set_num_threads(8);

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

    num_objects = std::stoi(argv[1]);
    num_iterations = std::stoi(argv[2]);
    seed = std::stoull(argv[3]);
    size_enclosure = std::stod(argv[4]);
    time_step = std::stod(argv[5]);

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
    objects.resize(num_objects);
    auto rnd_object = [&gen, &uniform_distr, &normal_distr] {
        double x = uniform_distr(gen);
        double y = uniform_distr(gen);
        double z = uniform_distr(gen);
        double m = normal_distr(gen);
        return Object(m, x, y, z);
    };
    std::generate(objects.begin(), objects.end(), rnd_object);

    // Check for collisions before starting
    checkCollisions();

    // Print the initial config
    std::ofstream initial;
    initial.open("init_config.txt");
    initial << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << objects.size() << "\n";
    for (const auto& i : objects) {
        initial << std::fixed << std::setprecision(3) << i.p[0] << " " << i.p[1] << " " << i.p[2] << " " << i.v[0] << " " << i.v[1] << " " << i.v[2] << " " << i.mass << "\n";
    }
    initial.close();

    // Time loop
    for (size_t iteration = 0; iteration < (unsigned)num_iterations; iteration++) {
        updateObjects();

        // Check for collisions (for all objects j < i)
        checkCollisions();

    }  // END OF TIME LOOP

    // Printing final config
    std::ofstream final;
    final.open("final_config.txt");
    final << std::fixed << std::setprecision(3) << size_enclosure << " " << time_step << " " << objects.size() << "\n";
    for (const auto& i : objects) {
        final << std::fixed << std::setprecision(3) << i.p[0] << " " << i.p[1] << " " << i.p[2] << " " << i.v[0] << " " << i.v[1] << " " << i.v[2] << " " << i.mass << "\n";
    }
    final.close();

    // Measure execution time and print it
    totalWatch.stop();
    if (en_benchmark) {
        std::printf("%f", totalWatch.getCount() / 1000000.0);
    } else {
        double updateObjRel = (double)updateObjWatch.getCount() / totalWatch.getCount() * 100;
        double collisionTimeRel = (double)collisionWatch.getCount() / totalWatch.getCount() * 100;
        std::printf("Total execution time: %.1fms: UpdateObjTime %.1fms (%.1f%%), CollisionTime %.1fms (%.1f%%), Others (%.1f%%)\n",
                    totalWatch.getCount() / 1000000.0,
                    updateObjWatch.getCount() / 1000000.0,
                    updateObjRel,
                    collisionWatch.getCount() / 1000000.0,
                    collisionTimeRel,
                    100.0 - updateObjRel - collisionTimeRel);
    }

    return 0;
}
