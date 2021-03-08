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
    this->particlesCount = particlesPerCellCount * grid.getXCellNumber() * grid.getYCellNumber();
}

size_t Particles::getParticlesCount() {
    return particlesCount;
}

std::vector<vcl::vec2> Particles::getParticlePositions() {
    return positions;
}

std::vector<vcl::vec2> Particles::getParticleVelocities() {
    return velocities;
}

void Particles::toGrid() {
    //to check indexing, looks problematic
    auto u = grid.getU();
    u.fill(0);

    vcl::grid_2D<float> weights(grid.getYCellNumber() + 1, grid.getXCellNumber() + 1);
    weights.fill(0);
    for (size_t i = 0; i < positions.size(); i++) {
        auto &[x, y] = positions[i];
        auto xCoord = grid.barycentricOnX(x);
        auto yCoord = grid.barycentricOnY(y);
        addPointToInterpolation(u, weights, velocities[i][0], xCoord, yCoord);
    }
    for (size_t i = 0; i < u.dimension[0]; i++) {
        for (size_t j = 0; j < u.dimension[1]; j++) {
            if (u(i, j) > 0) {
                u(i, j) /= weights(i, j);
            }
        }
    }

    auto v = grid.getV();
    v.fill(0);

    weights.fill(0);

    for (size_t i = 0; i < positions.size(); i++) {
        auto &[x, y] = positions[i];
        auto xCoord = grid.barycentricOnX(x);
        auto yCoord = grid.barycentricOnY(y);
        addPointToInterpolation(v, weights, velocities[i][1], xCoord, yCoord);
    }
    for (size_t i = 0; i < v.dimension[0]; i++) {
        for (size_t j = 0; j < v.dimension[1]; j++) {
            if (v(i, j) > 0) {
                v(i, j) /= weights(i, j);
            }
        }
    }
}

void Particles::addPointToInterpolation(vcl::grid_2D<float> &field, vcl::grid_2D<float> &weight, float value,
                                        barycentricCoordinate xCoord, barycentricCoordinate yCoord) const {
    auto &[xIndex, xOffset] = xCoord;
    auto &[yIndex, yOffset] = yCoord;

    float coef = (1 - xOffset) * (1 - yOffset);
    field(yIndex, xIndex) += coef * value;
    weight(yIndex, xIndex) += coef;

    coef = xOffset * (1 - yOffset);
    field(yIndex + 1, xIndex) += coef * value;
    weight(yIndex + 1, xIndex) += coef;

    coef = (1 - xOffset) * yOffset;
    field(yIndex, xIndex + 1) += coef * value;
    weight(yIndex, xIndex + 1) += coef;

    coef = xOffset * yOffset;
    field(yIndex + 1, xIndex + 1) += coef * value;
    weight(yIndex + 1, xIndex + 1) += coef;

}

void Particles::updateExternalForces(float dt) {
    const float g = 9.8;
    std::cout << "test: " << dt * g << std::endl;
    grid.getV() -= dt * g;
}

void Particles::fromGrid() {
    for (size_t i = 0; i < positions.size(); i++) {
        vcl::vec2 positionIndexed = clampPosAccordingToGrid(grid.getU(), positions[i]);
        float u = vcl::interpolation_bilinear(grid.getU(), positionIndexed[0], positionIndexed[1]);
        positionIndexed = clampPosAccordingToGrid(grid.getV(), positions[i]);
        float v = vcl::interpolation_bilinear(grid.getV(), positionIndexed[0], positionIndexed[1]);
        velocities[i] = {u, v}; //PIC step
    }
}

void Particles::step(float dt) {
    moveParticles(dt);
    toGrid();
    updateExternalForces(dt);
    grid.updateDistanceField();
    grid.interpolateVelocityWithFastSweep();
    grid.updateBoundaries();
    grid.divFreeField();
    fromGrid();
}

void Particles::moveParticles(float dt) {
    for (auto &position : positions) {
        vcl::vec2 positionIndexedU = clampPosAccordingToGrid(grid.getU(), position);
        vcl::vec2 positionIndexedV = clampPosAccordingToGrid(grid.getV(), position);
        float u = vcl::interpolation_bilinear(grid.getU(), positionIndexedU[0], positionIndexedU[1]);
        float v = vcl::interpolation_bilinear(grid.getV(), positionIndexedV[0], positionIndexedV[1]);
        vcl::vec2 firstStep = position + 0.5f * dt * vcl::vec2{u, v};
        positionIndexedU = clampPosAccordingToGrid(grid.getU(), firstStep);
        positionIndexedV = clampPosAccordingToGrid(grid.getV(), firstStep);
        u = vcl::interpolation_bilinear(grid.getU(), positionIndexedU[0], positionIndexedU[1]);
        v = vcl::interpolation_bilinear(grid.getV(), positionIndexedV[0], positionIndexedV[1]);
        position += dt * vcl::vec2{u, v};
    }
}

vcl::vec2 Particles::clampPosAccordingToGrid(const vcl::grid_2D<float> &grid, const vcl::vec2 &pos) const{
    float cellSize = this->grid.getCellSize();
    return {std::clamp(pos[0]/cellSize, 0.f, static_cast<float>(grid.dimension[0]-1.001)), // 1.001 is sort of workaround
            std::clamp(pos[1]/cellSize, 0.f, static_cast<float>(grid.dimension[1]-1.001))};

}
