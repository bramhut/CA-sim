#pragma once

#include <chrono>

// Watch class used for easy benchmarking
class watch {
    std::chrono::steady_clock::time_point t1;
    std::chrono::steady_clock::time_point t2;
    uint64_t count = 0;
public:
    watch() {
        start();
    }
    uint64_t stop() {
        t2 = std::chrono::steady_clock::now();
        count += std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        return count;
    }
    void start() {
        t1 = std::chrono::steady_clock::now();
    }
    uint64_t getCount() {
        return count;
    }
};