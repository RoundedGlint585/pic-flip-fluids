//
// Created by roundedglint585 on 3/6/21.
//

#include "MACGrid.h"

#include <cmath>

MACGrid::MACGrid(size_t xCellNumber, size_t yCellNumber, float cellSize): xCellNumber(xCellNumber), yCellNumber(yCellNumber), cellSize(cellSize){
    pressure = vcl::grid_2D<float>(xCellNumber, yCellNumber);
    pressure.fill(0);
    density = vcl::grid_2D<float>(xCellNumber, yCellNumber);
    density.fill(0);

    staggeredHorizontal = vcl::grid_2D<float>(xCellNumber, yCellNumber + 1);
    staggeredVertically = vcl::grid_2D<float>(xCellNumber + 1, yCellNumber);

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

barycentricCoordinate MACGrid::barycentricOnAxis(float coordinate) const {
    float cellIndex = coordinate/cellSize;
    return {static_cast<size_t>(cellIndex), cellIndex - std::floor(cellIndex) };
}

