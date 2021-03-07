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
    auto u = grid.getU();
    u.fill(0);

    vcl::grid_2D<float> weights(grid.getYCellNumber() + 1, grid.getXCellNumber() + 1);
    weights.fill(0);
    for (size_t i = 0; i < positions.size(); i++) {
        auto &[x,y] = positions[i];
        auto xCoord = grid.barycentricOnX(x);
        auto yCoord = grid.barycentricOnY(y);
        addPointToInterpolation(u, weights, velocities[i][0], xCoord, yCoord);
    }
    for(size_t i = 0; i < u.dimension[0]; i++){
        for(size_t j = 0; j < u.dimension[1]; j++){
            if(u(i, j) > 0){
                u(i, j) /= weights(i, j);
            }
        }
    }

    auto v = grid.getV();
    v.fill(0);

    weights.fill(0);

    for (size_t i = 0; i < positions.size(); i++) {
        auto &[x,y] = positions[i];
        auto xCoord = grid.barycentricOnX(x);
        auto yCoord = grid.barycentricOnY(y);
        addPointToInterpolation(v, weights, velocities[i][1], xCoord, yCoord);
    }
    for(size_t i = 0; i < v.dimension[0]; i++){
        for(size_t j = 0; j < v.dimension[1]; j++){
            if(v(i, j) > 0){
                v(i, j) /= weights(i, j);
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
    grid.getV() -= dt * g;
}

void Particles::fromGrid() {
    for(size_t i = 0; i < positions.size(); i++){
        float u = vcl::interpolation_bilinear(grid.getU(), positions[i][0], positions[i][1]);
        float v = vcl::interpolation_bilinear(grid.getV(), positions[i][0], positions[i][1]);
        velocities[i] = {u, v}; //PIC step
    }
}

void Particles::step(float dt) {

}
