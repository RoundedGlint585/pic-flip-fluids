//
// Created by roundedglint585 on 3/6/21.
//

#include "MACGrid.h"
#include "boundaries.h"

#include <cmath>

MACGrid::MACGrid(size_t xCellNumber, size_t yCellNumber, float cellSize) : xCellNumber(xCellNumber),
                                                                           yCellNumber(yCellNumber),
                                                                           cellSize(cellSize) {
    pressure = vcl::grid_2D<float>(xCellNumber, yCellNumber);
    pressure.fill(0);
    density = vcl::grid_2D<float>(xCellNumber, yCellNumber);
    density.fill(0);
    cellTypes = vcl::grid_2D<cellType>(xCellNumber, yCellNumber);
    cellTypes.fill(cellType::EMPTY_CELL);
    distanceField = vcl::grid_2D<float>(xCellNumber, yCellNumber);
    distanceField.fill(0);

    v = vcl::grid_2D<float>(xCellNumber, yCellNumber + 1);
    u = vcl::grid_2D<float>(xCellNumber + 1, yCellNumber);

}

size_t MACGrid::getXCellNumber() const {
    return xCellNumber;
}

size_t MACGrid::getYCellNumber() const {
    return yCellNumber;
}

float MACGrid::getCellSize() const {
    return cellSize;
}

barycentricCoordinate MACGrid::barycentricOnX(float x) const {
    float cellsCoord = x / cellSize;
    int cellIndex = static_cast<int>(cellsCoord);
    if (cellIndex < 0) {
        return {0, 0.f};
    } else if (cellIndex > xCellNumber - 2) {
        return {xCellNumber - 2, 1.f};
    } else {
        return {static_cast<size_t>(cellsCoord), cellsCoord - std::floor(cellsCoord)};
    }
}

barycentricCoordinate MACGrid::barycentricOnY(float y) const {
    float cellsCoord = y / cellSize;
    int cellIndex = static_cast<int>(cellsCoord);
    if (cellIndex < 0) {
        return {0, 0.f};
    } else if (cellIndex > yCellNumber - 2) {
        return {yCellNumber - 2, 1.f};
    } else {
        return {static_cast<size_t>(cellsCoord), cellsCoord - std::floor(cellsCoord)};
    }
}

vcl::grid_2D<float> &MACGrid::getU() {
    return u;
}

vcl::grid_2D<float>& MACGrid::getV() {
    return v;
}

vcl::grid_2D<float>& MACGrid::getDensity() {
    return density;
}

vcl::grid_2D<cellType> &MACGrid::getCellTypes() {
    return cellTypes;
}

void MACGrid::updateDistanceField() {
    //expect that cell types already updated
    distanceField.fill(std::numeric_limits<float>::max());//fill with inf
    for (size_t i = 0; i < cellTypes.dimension[0]; i++) {
        for (size_t j = 0; j < cellTypes.dimension[1]; j++) {
            if (cellTypes(i, j) == cellType::FLUID_CELL) {
                distanceField(i, j) = 0.0;
            }
        }
    }

    //fast sweeping to update field
    auto updateFunc = [this](int i, int j, int di, int dj) {
        if (cellTypes(i, j) == cellType::FLUID_CELL) {
            float a = distanceField(i - di, j);
            float b = distanceField(i, j - dj);
            float initDistance = std::min(a, b) + 1;
            if (initDistance > std::max(a, b)) {
                initDistance = (a + b + std::sqrt(2.f - std::pow<float>(a - b, 2))) / 2;
            }
            distanceField(i, j) = std::min(distanceField(i, j), initDistance);
        }
    };
    const size_t iterationCount = 2;
    for (size_t iter = 0; iter < iterationCount; iter++) {
        performSweep(1, distanceField.dimension[0], 1, distanceField.dimension[1], updateFunc);
        performSweep(1, distanceField.dimension[0], distanceField.dimension[1] - 2, 0, updateFunc);
        performSweep(distanceField.dimension[0] - 2, 0, 1, distanceField.dimension[1], updateFunc);
        performSweep(distanceField.dimension[0] - 2, 0, distanceField.dimension[1] - 2, 0, updateFunc);
    }
}

void MACGrid::interpolateVelocityWithFastSweep() {
    sweepHorizontal(2);
    sweepVertical(2);
}

void MACGrid::sweepHorizontal(size_t iterationCount) {

    auto updateFunc = [this](int i, int j, int di, int dj) {
        if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i - 1, j) == cellType::EMPTY_CELL) {
            float distanceDeltaI = distanceField(i, j) - distanceField(i - di, j);
            if (distanceDeltaI < 0) {
                return;
            }
            float distanceDeltaJ = distanceField(i, j) - distanceField(i, j - dj);
            if (distanceDeltaJ < 0) {
                return;
            }
            float coeff =
                    (distanceDeltaI + distanceDeltaJ) < std::numeric_limits<float>::epsilon() ? 0.5 : distanceDeltaI /
                                                                                                      (distanceDeltaI +
                                                                                                       distanceDeltaJ);
            u(i, j) =
                    coeff * u(i - di, j) + (1 - coeff) * u(i, j - dj);
        }
    };

    for (size_t k = 0; k < iterationCount; k++) {
        performSweep(1, distanceField.dimension[0], 1, distanceField.dimension[1], updateFunc);
        performSweep(1, distanceField.dimension[0], distanceField.dimension[1] - 2, 0, updateFunc);
        performSweep(distanceField.dimension[0] - 2, 0, 1, distanceField.dimension[1], updateFunc);
        performSweep(distanceField.dimension[0] - 2, 0, distanceField.dimension[1] - 2, 0, updateFunc);
    }
    for (size_t i = 0; i < u.dimension[0]; i++) {
        u(i, 0u) = u(i, 1u);
        u(i, u.dimension[1] - 1) = u(i, u.dimension[1] - 2);
    }
    for (size_t i = 0; i < u.dimension[1]; i++) {
        u(0u, i) = u(1u, i);
        u(u.dimension[0] - 1, i) = u(u.dimension[0] - 2, i);
    }
}

void MACGrid::sweepVertical(size_t iterationCount) {

    auto updateFunc = [this](int i, int j, int di, int dj) {
        if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i, j - 1) == cellType::EMPTY_CELL) {
            float distanceDeltaI = distanceField(i, j) - distanceField(i, j - dj);
            if (distanceDeltaI < 0) {
                return;
            }
            float distanceDeltaJ = distanceField(i, j) - distanceField(i - di, j);
            if (distanceDeltaJ < 0) {
                return;
            }
            float coeff =
                    (distanceDeltaI + distanceDeltaJ) < std::numeric_limits<float>::epsilon() ? 0.5 : distanceDeltaI /
                                                                                                      (distanceDeltaI +
                                                                                                       distanceDeltaJ);
            v(i, j) =
                    coeff * v(i - di, j) + (1 - coeff) * v(i, j - dj);
        }
    };

    for (size_t k = 0; k < iterationCount; k++) {
        performSweep(1, distanceField.dimension[0], 1, distanceField.dimension[1], updateFunc);
        performSweep(1, distanceField.dimension[0], distanceField.dimension[1] - 2, 0, updateFunc);
        performSweep(distanceField.dimension[0] - 2, 0, 1, distanceField.dimension[1], updateFunc);
        performSweep(distanceField.dimension[0] - 2, 0, distanceField.dimension[1] - 2, 0, updateFunc);
    }

    for (size_t i = 0; i < v.dimension[0]; i++) {
        v(i, 0u) = v(i, 1u);
        v(i, v.dimension[1] - 1) = v(i, v.dimension[1] - 2);
    }
    for (size_t i = 0; i < u.dimension[1]; i++) {
        v(0u, i) = v(1u, i);
        v(v.dimension[0] - 1, i) = v(v.dimension[0] - 2, i);
    }
}

void MACGrid::updateBoundaries() {
    for (size_t i = 0; i < cellTypes.dimension[0]; i++) {
        cellTypes(i, 0u) = cellType::SOLID_CELL;
        cellTypes(i, cellTypes.dimension[1] - 1) = cellType::SOLID_CELL;
    }
    for (size_t i = 0; i < cellTypes.dimension[1]; i++) {
        cellTypes(0u, i) = cellType::SOLID_CELL;
        cellTypes(cellTypes.dimension[0] - 1, i) = cellType::SOLID_CELL;
    }

    for (size_t i = 0; i < u.dimension[1]; i++) {
        u(0u, i) = 0;
        u(1u, i) = 0;
        u(u.dimension[0]-1, i) = 0;
        u(u.dimension[0]-2, i) = 0;
    }
    for (size_t i = 0; i < v.dimension[0]; i++) {
        v(i, 0u) = 0;
        v(i, 1u) = 0;
        v(i, v.dimension[1] - 1) = 0;
        v(i, v.dimension[1] - 2) = 0;
    }
}

vcl::grid_2D<float> MACGrid::getDivergence() const {
    auto div = vcl::grid_2D<float>(yCellNumber, xCellNumber);
    div.fill(0);
    for (size_t j = 0; j < div.dimension[1]; j++) {
        for (size_t i = 0; i < div.dimension[0]; i++) {
            if (cellTypes(i, j) == cellType::FLUID_CELL) {
                div(i, j) = u(i + 1, j) - u(i, j) + v(i, j + 1) - v(i, j);
            }
        }
    }
    return div;
}

void
MACGrid::performSweep(int fromX, int toX, int fromY, int toY, const std::function<void(int, int, int, int)> &function) {
    int di = fromX <= toX ? 1 : -1;
    int dj = fromY <= toY ? 1 : -1;
    for (int j = fromY; j != toY; j += dj) {
        for (int i = fromX; i != toX; i += di) {
            function(i, j, di, dj);
        }
    }
}

void MACGrid::divFreeField() {
    auto div = getDivergence();
    // Gauss Seidel
    vcl::grid_2D<float> q = vcl::grid_2D<float>(div.dimension[0], div.dimension[1]);
    for (size_t k_iter = 0; k_iter < 20; ++k_iter) {
        for (size_t x = 1; x < div.dimension[0] - 1; ++x) {
            for (size_t y = 1; y < div.dimension[1] - 1; ++y) {
                q(x, y) = (q(x + 1, y) + q(x - 1, y) + q(x, y + 1) + q(x, y - 1) - div(x, y)) / 4.0f;
            }
        }
        set_boundary(q);
    }

    for (size_t x = 1; x < u.dimension[0] - 2; ++x) {
        for (size_t y = 1; y < u.dimension[1] - 1; ++y) {
            u(x, y) -= (q(x + 1, y) - q(x - 1, y));
        }
    }
    for (size_t x = 1; x < v.dimension[0] - 1; ++x) {
        for (size_t y = 1; y < v.dimension[1] - 2; ++y) {
            v(x, y) = v(x, y) - (q(x, y + 1) - q(x, y - 1));
        }
    }
    set_boundary(u);
    set_boundary(v);
}

barycentricCoordinate MACGrid::barycentricOnXUnsafe(float x) const {
    float cellsCoord = x / cellSize;
    int cellIndex = static_cast<int>(cellsCoord);
    return {static_cast<size_t>(cellsCoord), cellsCoord - std::floor(cellsCoord)};
}

barycentricCoordinate MACGrid::barycentricOnYUnsafe(float y) const{
    float cellsCoord = y / cellSize;
    int cellIndex = static_cast<int>(cellsCoord);
    return {static_cast<size_t>(cellsCoord), cellsCoord - std::floor(cellsCoord)};
}

