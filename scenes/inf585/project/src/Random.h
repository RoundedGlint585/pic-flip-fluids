//
// Created by roundedglint585 on 3/6/21.
//

#ifndef PROJECT_RANDOM_H
#define PROJECT_RANDOM_H

#include <random>
#include "Image.h"

namespace randomGenerator {
    static thread_local std::mt19937 generator;

    inline float floatRandom(const float min, const float max) {
        std::uniform_real_distribution<float> distribution(min, max);
        return distribution(generator);
    }

}
#endif //PROJECT_RANDOM_H
