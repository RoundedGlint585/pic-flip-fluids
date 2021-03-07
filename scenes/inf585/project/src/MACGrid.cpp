//
// Created by roundedglint585 on 3/6/21.
//

#include "MACGrid.h"

#include <cmath>

MACGrid::MACGrid(size_t xCellNumber, size_t yCellNumber, float cellSize) : xCellNumber(xCellNumber),
                                                                           yCellNumber(yCellNumber),
                                                                           cellSize(cellSize) {
    pressure = vcl::grid_2D<float>(yCellNumber, xCellNumber);
    pressure.fill(0);
    density = vcl::grid_2D<float>(yCellNumber, xCellNumber);
    density.fill(0);
    cellTypes = vcl::grid_2D<cellType>(yCellNumber, xCellNumber);
    cellTypes.fill(cellType::EMPTY_CELL);

    staggeredHorizontal = vcl::grid_2D<float>(yCellNumber + 1, xCellNumber);
    staggeredVertical = vcl::grid_2D<float>(yCellNumber, xCellNumber + 1);

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
    float cellsCoord = x / cellSize - 0.5 * cellSize;
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
    float cellsCoord = y / cellSize - 0.5 * cellSize;
    int cellIndex = static_cast<int>(cellsCoord);
    if (cellIndex < 0) {
        return {0, 0.f};
    } else if (cellIndex > yCellNumber - 2) {
        return {yCellNumber - 2, 1.f};
    } else {
        return {static_cast<size_t>(cellsCoord), cellsCoord - std::floor(cellsCoord)};
    }
}

vcl::grid_2D<float> &MACGrid::getStaggeredHorizontal() {
    return staggeredHorizontal;
}

vcl::grid_2D<float> &MACGrid::getStaggeredVertical() {
    return staggeredVertical;
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
    auto shortestDistance = [](float a, float b) -> float {
        float initDistance = std::min(a, b) + 1;
        if (initDistance > std::max(a, b)) {
            return (a + b + std::sqrt(2.f - std::pow<float>(a - b, 2))) / 2;
        } else {
            return initDistance;
        }
    };
    //fast sweeping to update field

    const size_t iterationCount = 2;
    for (size_t iter = 0; iter < iterationCount; iter++) {
        for (size_t i = 1; i < distanceField.dimension[0]; i++) {
            for (size_t j = 1; j < distanceField.dimension[1]; j++) {
                if (cellTypes(i, j) == cellType::FLUID_CELL) {
                    distanceField(i, j) = std::min(distanceField(i, j),
                                                   shortestDistance(distanceField(i - 1, j), distanceField(i, j - 1)));
                }
            }
        }
        for (size_t i = distanceField.dimension[0] - 2; i >= 0; i--) {
            for (size_t j = 1; j < distanceField.dimension[1]; j++) {
                if (cellTypes(i, j) == cellType::FLUID_CELL) {
                    distanceField(i, j) = std::min(distanceField(i, j),
                                                   shortestDistance(distanceField(i + 1, j), distanceField(i, j - 1)));
                }
            }
        }

        for (size_t i = 1; i < distanceField.dimension[0]; i++) {
            for (size_t j = distanceField.dimension[1] - 2; j >= 0; j--) {
                if (cellTypes(i, j) == cellType::FLUID_CELL) {
                    distanceField(i, j) = std::min(distanceField(i, j),
                                                   shortestDistance(distanceField(i - 1, j), distanceField(i, j + 1)));
                }
            }
        }

        for (size_t i = distanceField.dimension[0] - 2; i >= 0; i--) {
            for (size_t j = distanceField.dimension[1] - 2; j >= 0; j--) {
                if (cellTypes(i, j) == cellType::FLUID_CELL) {
                    distanceField(i, j) = std::min(distanceField(i, j),
                                                   shortestDistance(distanceField(i + 1, j), distanceField(i, j + 1)));
                }
            }
        }
    }
}

void MACGrid::interpolateVelocityWithFastSweep() {
    sweepHorizontal(2);
    sweepVertical(2);
}

void MACGrid::sweepHorizontal(size_t iterationCount) {
    for(size_t k = 0; k < iterationCount; k++) {
        for (size_t i = 1; i < staggeredHorizontal.dimension[0]; i++) {
            for (size_t j = 1; j < staggeredHorizontal.dimension[1]; j++) {
                if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i - 1, j) == cellType::EMPTY_CELL) {
                    float di = distanceField(i, j) - distanceField(i - 1, j);
                    if (di < 0) {
                        continue;
                    }
                    float dj = distanceField(i, j) - distanceField(i, j - 1);
                    if (dj < 0) {
                        continue;
                    }
                    float coeff = (di + dj) < std::numeric_limits<float>::epsilon() ? 0.5 : di / (di + dj);
                    staggeredHorizontal(i, j) =
                            coeff * staggeredHorizontal(i - 1, j) + (1 - coeff) * staggeredHorizontal(i, j - 1);
                }
            }
        }
        for (size_t i = staggeredHorizontal.dimension[0] - 2; i >= 0; i--) {
            for (size_t j = 1; j < staggeredHorizontal.dimension[1]; j++) {
                if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i - 1, j) == cellType::EMPTY_CELL) {
                    float di = distanceField(i, j) - distanceField(i + 1, j);
                    if (di < 0) {
                        continue;
                    }
                    float dj = distanceField(i, j) - distanceField(i, j - 1);
                    if (dj < 0) {
                        continue;
                    }
                    float coeff = (di + dj) < std::numeric_limits<float>::epsilon() ? 0.5 : di / (di + dj);
                    staggeredHorizontal(i, j) =
                            coeff * staggeredHorizontal(i + 1, j) + (1 - coeff) * staggeredHorizontal(i, j - 1);
                }
            }
        }

        for (size_t i = 1; i < staggeredHorizontal.dimension[0]; i++) {
            for (size_t j = staggeredHorizontal.dimension[1] - 2; j >= 0; j--) {
                if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i - 1, j) == cellType::EMPTY_CELL) {
                    float di = distanceField(i, j) - distanceField(i - 1, j);
                    if (di < 0) {
                        continue;
                    }
                    float dj = distanceField(i, j) - distanceField(i, j + 1);
                    if (dj < 0) {
                        continue;
                    }
                    float coeff = (di + dj) < std::numeric_limits<float>::epsilon() ? 0.5 : di / (di + dj);
                    staggeredHorizontal(i, j) =
                            coeff * staggeredHorizontal(i - 1, j) + (1 - coeff) * staggeredHorizontal(i, j + 1);
                }
            }
        }

        for (size_t i = staggeredHorizontal.dimension[0] - 2; i >= 0; i--) {
            for (size_t j = staggeredHorizontal.dimension[1] - 2; j >= 0; j--) {
                if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i + 1, j) == cellType::EMPTY_CELL) {
                    float di = distanceField(i, j) - distanceField(i + 1, j);
                    if (di < 0) {
                        continue;
                    }
                    float dj = distanceField(i, j) - distanceField(i, j + 1);
                    if (dj < 0) {
                        continue;
                    }
                    float coeff = (di + dj) < std::numeric_limits<float>::epsilon() ? 0.5 : di / (di + dj);
                    staggeredHorizontal(i, j) =
                            coeff * staggeredHorizontal(i + 1, j) + (1 - coeff) * staggeredHorizontal(i, j + 1);
                }
            }
        }
    }
    for (size_t i = 0; i < staggeredHorizontal.dimension[1]; i++) {
        staggeredHorizontal(i, 0u) = staggeredHorizontal(i, 1u);
        staggeredHorizontal(i, staggeredHorizontal.dimension[0] - 1) = staggeredHorizontal(i, staggeredHorizontal[0] - 2);
    }
    for(size_t i = 0; i < staggeredHorizontal.dimension[0]; i++){
        staggeredHorizontal(0u, i) = staggeredHorizontal(1u, i);
        staggeredHorizontal(staggeredHorizontal.dimension[1] - 1, i) = staggeredHorizontal(staggeredHorizontal.dimension[1] - 2, i);
    }
}

void MACGrid::sweepVertical(size_t iterationCount) {

    for (size_t i = 1; i < staggeredVertical.dimension[0]; i++) {
        for (size_t j = 1; j < staggeredVertical.dimension[1]; j++) {
            if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i - 1, j) == cellType::EMPTY_CELL) {
                float di = distanceField(i, j) - distanceField(i, j-1);
                if (di < 0) {
                    continue;
                }
                float dj = distanceField(i, j) - distanceField(i-1, j);
                if (dj < 0) {
                    continue;
                }
                float coeff = (di + dj) < std::numeric_limits<float>::epsilon() ? 0.5 : di / (di + dj);
                staggeredVertical(i, j) =
                        coeff * staggeredVertical(i - 1, j) + (1 - coeff) * staggeredVertical(i, j - 1);
            }
        }
    }
    for (size_t i = staggeredVertical.dimension[0] - 2; i >= 0; i--) {
        for (size_t j = 1; j < staggeredVertical.dimension[1]; j++) {
            if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i - 1, j) == cellType::EMPTY_CELL) {
                float di = distanceField(i, j) - distanceField(i, j - 1);
                if (di < 0) {
                    continue;
                }
                float dj = distanceField(i, j) - distanceField(i + 1, j);
                if (dj < 0) {
                    continue;
                }
                float coeff = (di + dj) < std::numeric_limits<float>::epsilon() ? 0.5 : di / (di + dj);
                staggeredVertical(i, j) =
                        coeff * staggeredVertical(i + 1, j) + (1 - coeff) * staggeredVertical(i, j - 1);
            }
        }
    }

    for (size_t i = 1; i < staggeredVertical.dimension[0]; i++) {
        for (size_t j = staggeredVertical.dimension[1] - 2; j >= 0; j--) {
            if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i - 1, j) == cellType::EMPTY_CELL) {
                float di = distanceField(i, j) - distanceField(i, j + 1);
                if (di < 0) {
                    continue;
                }
                float dj = distanceField(i, j) - distanceField(i - 1, j);
                if (dj < 0) {
                    continue;
                }
                float coeff = (di + dj) < std::numeric_limits<float>::epsilon() ? 0.5 : di / (di + dj);
                staggeredVertical(i, j) =
                        coeff * staggeredVertical(i - 1, j) + (1 - coeff) * staggeredVertical(i, j + 1);
            }
        }
    }

    for (size_t i = staggeredVertical.dimension[0] - 2; i >= 0; i--) {
        for (size_t j = staggeredVertical.dimension[1] - 2; j >= 0; j--) {
            if (cellTypes(i, j) == cellType::EMPTY_CELL && cellTypes(i + 1, j) == cellType::EMPTY_CELL) {
                float di = distanceField(i, j) - distanceField(i, j + 1);
                if (di < 0) {
                    continue;
                }
                float dj = distanceField(i, j) - distanceField(i + 1, j);
                if (dj < 0) {
                    continue;
                }
                float coeff = (di + dj) < std::numeric_limits<float>::epsilon() ? 0.5 : di / (di + dj);
                staggeredVertical(i, j) =
                        coeff * staggeredVertical(i + 1, j) + (1 - coeff) * staggeredVertical(i, j + 1);
            }
        }
    }
    for (size_t i = 0; i < staggeredVertical.dimension[1]; i++) {
        staggeredVertical(i, 0u) = staggeredVertical(i, 1u);
        staggeredVertical(i, staggeredVertical.dimension[0] - 1) = staggeredVertical(i, staggeredVertical[0] - 2);
    }
    for(size_t i = 0; i < staggeredVertical.dimension[0]; i++){
        staggeredVertical(0u, i) = staggeredVertical(1u, i);
        staggeredVertical(staggeredVertical.dimension[1] - 1, i) = staggeredVertical(staggeredVertical.dimension[1] - 2, i);
    }
}

void MACGrid::updateBoundaries() {
    for(size_t i = 0; i < cellTypes.dimension[0]; i++){
        cellTypes(i, 0u) = cellType::SOLID_CELL;
        cellTypes(i, cellTypes.dimension[1]-1) = cellType::SOLID_CELL;
    }
    for(size_t i = 0; i < cellTypes.dimension[1]; i++){
        cellTypes(0u, i) = cellType::SOLID_CELL;
        cellTypes(cellTypes.dimension[0]-1, i) = cellType::SOLID_CELL;
    }

    for(size_t i = 0; i < staggeredHorizontal.dimension[0]; i++){
        staggeredHorizontal(i, 0u) = 0;
        staggeredHorizontal(i, 1u) = 0;
        staggeredHorizontal(i, staggeredHorizontal.dimension[1] - 1) = 0;
        staggeredHorizontal(i, staggeredHorizontal.dimension[1] - 2) = 0;
    }
    for(size_t i = 0; i < staggeredHorizontal.dimension[0]; i++){
        staggeredHorizontal(i, 0u) = 0;
        staggeredHorizontal(i, 1u) = 0;
        staggeredHorizontal(i, staggeredHorizontal.dimension[1] - 1) = 0;
        staggeredHorizontal(i, staggeredHorizontal.dimension[1] - 2) = 0;
    }
    for(size_t i = 0; i < staggeredVertical.dimension[1]; i++){
        staggeredVertical(0u, i) = 0;
        staggeredVertical(1u, i) = 0;
        staggeredVertical(staggeredVertical.dimension[0] - 1, i) = 0;
        staggeredVertical(staggeredVertical.dimension[0] - 2, i) = 0;
    }
}
