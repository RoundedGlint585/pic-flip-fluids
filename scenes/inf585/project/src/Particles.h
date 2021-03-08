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
    size_t getParticlesCount();
    std::vector<vcl::vec2> getParticlePositions();
    std::vector<vcl::vec2> getParticleVelocities();

private:
    MACGrid grid;
    size_t particlesCount;
    std::vector<vcl::vec2> positions, velocities;

    void toGrid(); // update grid velocity based on particles around
    void updateExternalForces(float dt); // for now only gravity
    void fromGrid();
    void addPointToInterpolation(vcl::grid_2D<float> &field, vcl::grid_2D<float> &weight, float value, barycentricCoordinate xCoord, barycentricCoordinate yCoord) const;
    void moveParticles(float dt);
};


#endif //PROJECT_PARTICLES_H
