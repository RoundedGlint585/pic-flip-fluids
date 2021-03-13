//
// Created by roundedglint585 on 3/6/21.
//

#include "Particles.h"
#include "Random.h"


Particles::Particles(size_t particlesPerCellCount, const MACGrid &grid) : grid(grid) {
    float cellSize = grid.cellSize;
    this->grid.cellTypes.fill(cellType::EMPTY_CELL);

    for (size_t j = 10; j < grid.yCellCount - 5; j++) { //sry for this way of initialization
        for (size_t i = 2; i < grid.xCellCount / 2; i++) {
            float xMin = i * cellSize;
            float xMax = static_cast<float>(i + 1) * cellSize;
            float yMin = j * cellSize;
            float yMax = static_cast<float>(j + 1) * cellSize;
            for (size_t k = 0; k < particlesPerCellCount; k++) {
                positions.emplace_back(randomGenerator::floatRandom(xMin, xMax),
                                       randomGenerator::floatRandom(yMin, yMax));
                velocities.emplace_back(0.f, 0.f);
            }
            this->grid.cellTypes(i, j) = FLUID_CELL;
        }
    }
    this->grid.updateBoundaries();
}

void Particles::toGrid() {

    auto &u = grid.u;
    u.fill(0);

    vcl::grid_2D<float> weights(grid.cellTypes.dimension.x + 1, grid.cellTypes.dimension.y + 1);
    weights.fill(0);
    for (size_t i = 0; i < positions.size(); i++) {
        auto &[x, y] = positions[i];
        auto xCoord = grid.barycentricOnX(x);
        auto yCoord = grid.barycentricOnY(y);
        addPointToInterpolation(u, weights, velocities[i][0], xCoord, yCoord);
    }
    for (size_t i = 0; i < u.dimension[0]; i++) {
        for (size_t j = 0; j < u.dimension[1]; j++) {
            if (weights(i, j) != 0) {
                u(i, j) /= weights(i, j);
            }
        }
    }

    auto &v = grid.v;
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
            if (weights(i, j) != 0) {
                v(i, j) /= weights(i, j);
            }
        }
    }

    grid.cellTypes.fill(EMPTY_CELL);
    for (auto &[x, y]: positions) {
        auto xCoord = grid.barycentricOnX(x);
        auto yCoord = grid.barycentricOnY(y);
        grid.cellTypes(xCoord.index, yCoord.index) = FLUID_CELL;
    }

}

void Particles::fromGrid() {
    for (size_t i = 0; i < positions.size(); i++) {
        vcl::vec2 positionIndexed = clampPosAccordingToGrid(grid.u, positions[i]);
        float u = vcl::interpolation_bilinear(grid.u, positionIndexed.x, positionIndexed.y);
        float du = vcl::interpolation_bilinear(grid.du, positionIndexed.x, positionIndexed.y);
        positionIndexed = clampPosAccordingToGrid(grid.v, positions[i]);
        float v = vcl::interpolation_bilinear(grid.v, positionIndexed.x, positionIndexed.y);
        float dv = vcl::interpolation_bilinear(grid.dv, positionIndexed.x, positionIndexed.y);
        velocities[i] = 1 * vcl::vec2{u, v}; //PIC step
        velocities[i] += 1 * vcl::vec2{du, dv};
    }
}

void Particles::moveParticles(float dt) {
    for (auto &position : positions) {
        vcl::vec2 positionIndexedU = clampPosAccordingToGrid(grid.u, position);
        vcl::vec2 positionIndexedV = clampPosAccordingToGrid(grid.v, position);
        float u = vcl::interpolation_bilinear(grid.u, positionIndexedU.x, positionIndexedU.y);
        float v = vcl::interpolation_bilinear(grid.v, positionIndexedV.x, positionIndexedV.y);
        vcl::vec2 firstStep = position + 0.5f * dt * vcl::vec2{u, v};
        positionIndexedU = clampPosAccordingToGrid(grid.u, firstStep);
        positionIndexedV = clampPosAccordingToGrid(grid.v, firstStep);
        u = vcl::interpolation_bilinear(grid.u, positionIndexedU.x, positionIndexedU.y);
        v = vcl::interpolation_bilinear(grid.v, positionIndexedV.x, positionIndexedV.y);
        position += dt * vcl::vec2{u, v};
    }
}

void Particles::addPointToInterpolation(vcl::grid_2D<float> &field, vcl::grid_2D<float> &weight, float value,
                                        barycentricCoords xCoord, barycentricCoords yCoord) const {
    auto &[xIndex, xOffset] = xCoord;
    auto &[yIndex, yOffset] = yCoord;

    float coef = (1 - xOffset) * (1 - yOffset);
    field(xIndex, yIndex) += coef * value;
    weight(xIndex, yIndex) += coef;

    coef = xOffset * (1 - yOffset);
    field(xIndex + 1, yIndex) += coef * value;
    weight(xIndex + 1, yIndex) += coef;

    coef = (1 - xOffset) * yOffset;
    field(xIndex, yIndex + 1) += coef * value;
    weight(xIndex, yIndex + 1) += coef;

    coef = xOffset * yOffset;
    field(xIndex + 1, yIndex + 1) += coef * value;
    weight(xIndex + 1, yIndex + 1) += coef;
}

vcl::vec2 Particles::clampPosAccordingToGrid(const vcl::grid_2D<float> &grid, const vcl::vec2 &pos) const {
    return {std::clamp(pos[0] / this->grid.cellSize, 0.001f,
                       static_cast<float>(grid.dimension[0] - 1.001)), // 1.001 is sort of workaround
            std::clamp(pos[1] / this->grid.cellSize, 0.001f, static_cast<float>(grid.dimension[1] - 1.001))};

}

void Particles::step(float dt) {
    moveParticles(dt);
    toGrid();
    grid.update(dt);
    fromGrid();
}

std::vector<vcl::vec2> Particles::getParticlePositions() {
    return positions;
}

std::vector<vcl::vec2> Particles::getParticleVelocities() {
    return velocities;
}
