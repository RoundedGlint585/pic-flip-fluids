//
// Created by roundedglint585 on 3/6/21.
//

#include "Particles.h"
#include "Random.h"

Particles::Particles(size_t particlesPerCellCount, const MACGrid& grid) : grid(grid) {
    float cellSize = grid.getCellSize();
    for(size_t i = 0; i < grid.getYCellNumber(); i++){
        for(size_t j = 0; j < grid.getXCellNumber(); j++){
            //randomly generate particles per cell
            float xMin = j*cellSize;
            float xMax = static_cast<float>(j+1)*cellSize;
            float yMin = i * cellSize;
            float yMax = static_cast<float>(i+1) * cellSize;
            for(size_t k = 0; k < particlesPerCellCount; k++){
                positions.emplace_back(randomGenerator::floatRandom(xMin, xMax), randomGenerator::floatRandom(yMin, yMax));
                velocities.emplace_back(0.f, 0.f);
            }
        }
    }
}
