//
// Created by roundedglint585 on 3/3/21.
//

#ifndef PROJECT_SIMULATION_H
#define PROJECT_SIMULATION_H
#include "vcl/vcl.hpp"
#include "Particles.h"

// SPH Particle
/*struct particle_element
{
    vcl::vec3 p; // Position
    vcl::vec3 v; // Speed
    vcl::vec3 f; // Force

    float rho;      // density at this particle position
    float pressure; // pressure at this particle position

    particle_element() : p{0,0,0},v{0,0,0},f{0,0,0},rho(0),pressure(0) {}
};*/

void simulate(Particles particles, float dt);



#endif //PROJECT_SIMULATION_H
