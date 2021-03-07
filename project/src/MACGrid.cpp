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
    staggeredVertical = vcl::grid_2D<float>(xCellNumber + 1, yCellNumber);

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
    float cellsCoord = coordinate/cellSize;
    return {static_cast<size_t>(cellsCoord), cellsCoord - std::floor(cellsCoord) };
}

barycentricCoordinate MACGrid::barycentricOnX(float x) const {
    float cellsCoord = x / cellSize - 0.5 * cellSize;
    int cellIndex = static_cast<int>(cellsCoord);
    if(cellIndex < 0){
        return {0, 0.f};
    }else if(cellIndex > xCellNumber -2){
        return {xCellNumber - 2, 1.f};
    }else{
        return {static_cast<size_t>(cellsCoord), cellsCoord - std::floor(cellsCoord) };
    }
}

barycentricCoordinate MACGrid::barycentricOnY(float y) const {
    float cellsCoord = y / cellSize - 0.5 * cellSize;
    int cellIndex = static_cast<int>(cellsCoord);
    if(cellIndex < 0){
        return {0, 0.f};
    }else if(cellIndex > yCellNumber -2){
        return {yCellNumber - 2, 1.f};
    }else{
        return {static_cast<size_t>(cellsCoord), cellsCoord - std::floor(cellsCoord) };
    }
}

vcl::grid_2D<float> &MACGrid::getStaggeredHorizontal() {
    return staggeredHorizontal;
}

vcl::grid_2D<float> &MACGrid::getStaggeredVertical() {
    return staggeredVertical;
}
