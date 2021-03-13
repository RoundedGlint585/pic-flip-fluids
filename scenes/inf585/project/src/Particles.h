//
// Created by roundedglint585 on 3/6/21.
//

#ifndef PROJECT_PARTICLES_H
#define PROJECT_PARTICLES_H

#include "MACGrid.h"

class Particles {
public:

    Particles(size_t particlesPerCellCount, const MACGrid &grid_);

    void step(float dt);

    std::vector<vcl::vec2> getParticlePositions();

    std::vector<vcl::vec2> getParticleVelocities();

private:
    MACGrid grid;
    std::vector<vcl::vec2> positions, velocities; // positions and velocities

    void addPointToInterpolation(vcl::grid_2D<float> &field, vcl::grid_2D<float> &weight, float value,
                                 barycentricCoords xCoord, barycentricCoords yCoord) const;

    vcl::vec2 clampPosAccordingToGrid(const vcl::grid_2D<float> &grid, const vcl::vec2 &pos) const;

    void toGrid();

    void fromGrid();

    void moveParticles(float dt);

};


#endif //PROJECT_PARTICLES_H
