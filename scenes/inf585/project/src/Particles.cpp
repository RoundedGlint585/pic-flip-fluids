//
// Created by roundedglint585 on 3/6/21.
//

#include "Particles.h"
#include "Random.h"

Particles::Particles(size_t particlesPerCellCount, const MACGrid &grid) : grid(grid) {
    float cellSize = grid.getCellSize();
    for (size_t i = 0; i < grid.getYCellNumber(); i++) {
        for (size_t j = 0; j < grid.getXCellNumber(); j++) {
            //randomly generate particles per cell
            float xMin = j * cellSize;
            float xMax = static_cast<float>(j + 1) * cellSize;
            float yMin = i * cellSize;
            float yMax = static_cast<float>(i + 1) * cellSize;
            for (size_t k = 0; k < particlesPerCellCount; k++) {
                positions.emplace_back(randomGenerator::floatRandom(xMin, xMax),
                                       randomGenerator::floatRandom(yMin, yMax));
                velocities.emplace_back(0.f, 0.f);
            }
        }
    }
    this->grid.getCellTypes().fill(cellType::FLUID_CELL);
}

void Particles::toGrid() {
    //to check indexing, looks problematic
    auto staggeredHorizontal = grid.getStaggeredHorizontal();
    staggeredHorizontal.fill(0);

    vcl::grid_2D<float> weights(grid.getYCellNumber() + 1, grid.getXCellNumber() + 1);
    weights.fill(0);
    for (size_t i = 0; i < positions.size(); i++) {
        auto &[x,y] = positions[i];
        auto xCoord = grid.barycentricOnX(x);
        auto yCoord = grid.barycentricOnY(y);
        addPointToInterpolation(staggeredHorizontal, weights, velocities[i][0], xCoord, yCoord);
    }
    for(size_t i = 0; i < staggeredHorizontal.dimension[0]; i++){
        for(size_t j = 0; j < staggeredHorizontal.dimension[1]; j++){
            if(staggeredHorizontal(i, j) > 0){
                staggeredHorizontal(i, j) /= weights(i, j);
            }
        }
    }

    auto staggeredVertical = grid.getStaggeredVertical();
    staggeredVertical.fill(0);

    weights.fill(0);

    for (size_t i = 0; i < positions.size(); i++) {
        auto &[x,y] = positions[i];
        auto xCoord = grid.barycentricOnX(x);
        auto yCoord = grid.barycentricOnY(y);
        addPointToInterpolation(staggeredHorizontal, weights, velocities[i][1], xCoord, yCoord);
    }
    for(size_t i = 0; i < staggeredVertical.dimension[0]; i++){
        for(size_t j = 0; j < staggeredVertical.dimension[1]; j++){
            if(staggeredVertical(i, j) > 0){
                staggeredVertical(i, j) /= weights(i, j);
            }
        }
    }
}

void Particles::addPointToInterpolation(vcl::grid_2D<float> &field, vcl::grid_2D<float> &weight, float value, barycentricCoordinate xCoord, barycentricCoordinate yCoord) const {
    auto &[xIndex, xOffset] = xCoord;
    auto &[yIndex, yOffset] = yCoord;

    float coef = (1 - xOffset) * (1 - yOffset);
    field(yIndex, xIndex) += coef * value;
    weight(yIndex, xIndex) += coef;

    coef = xOffset * (1 - yOffset);
    field(yIndex + 1, xIndex) += coef * value;
    weight(yIndex + 1, xIndex) += coef;

    coef = (1 - xOffset) * yOffset;
    field(yIndex, xIndex+1) += coef * value;
    weight(yIndex, xIndex + 1) += coef;

    coef = xOffset * yOffset;
    field(yIndex + 1, xIndex+1) += coef * value;
    weight(yIndex + 1, xIndex + 1) += coef;

}

void Particles::updateExternalForces(float dt) {
    const float g = 9.8;
    grid.getStaggeredHorizontal() -= dt * g;
}
