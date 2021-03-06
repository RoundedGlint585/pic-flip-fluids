//
// Created by roundedglint585 on 3/6/21.
//

#ifndef PROJECT_PARTICLES_H
#define PROJECT_PARTICLES_H
#include "MACGrid.h"

class Particles {
public:
    Particles(size_t particlesPerCellCount, const MACGrid &grid);

    void step(float dt);
private:
    MACGrid grid;
    size_t particlesCount;
    std::vector<vcl::vec2> positions, velocities;

    void toGrid();
    void fromGrid();
};


#endif //PROJECT_PARTICLES_H
