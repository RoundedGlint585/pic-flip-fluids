//
// Created by roundedglint585 on 3/6/21.
//

#ifndef PROJECT_MACGRID_H
#define PROJECT_MACGRID_H


#include <cstddef>
#include <functional>
#include "vcl/vcl.hpp"

struct barycentricCoordinate{
    size_t index;
    float dist;
};
enum cellType{
    EMPTY_CELL,
    FLUID_CELL,// cell with at least one particle
    SOLID_CELL,
};
class MACGrid {
public:
    MACGrid(size_t xCellNumber=30, size_t yCellNumber=30, float cellSize=2);
    size_t getXCellNumber() const;
    size_t getYCellNumber() const;
    float getCellSize() const;
    vcl::grid_2D<float>& getU();
    vcl::grid_2D<float>& getV();
    vcl::grid_2D<float>& getDensity();
    vcl::grid_2D<cellType> &getCellTypes();

    barycentricCoordinate barycentricOnXUnsafe(float x) const;
    barycentricCoordinate barycentricOnYUnsafe(float x) const;
    barycentricCoordinate barycentricOnX(float x) const; // for Mac grid
    barycentricCoordinate barycentricOnY(float y) const; // TODO: sort of similar code, take a look on unification
    void updateDistanceField();
    void interpolateVelocityWithFastSweep();
    void updateBoundaries();

    vcl::grid_2D<float> getDivergence() const;
    void divFreeField();
    // this shit is tricky, i don't want to fix every part of fast sweeping, so this thing should do for in for and apply function
    // based on i, j, di, dj, basically just sent lambda
    void performSweep(int fromX, int toX, int fromY, int toY, const std::function<void(int, int, int, int)>& function);
private:

    void sweepHorizontal(size_t iterationCount);
    void sweepVertical(size_t iterationCount);
    vcl::grid_2D<float> pressure, density; //store in the center
    vcl::grid_2D<float> u, v; //velocity stored on the edges
    vcl::grid_2D<cellType> cellTypes;
    vcl::grid_2D<float> distanceField; // phi field from paper created using fast sweep
    size_t xCellNumber, yCellNumber;
    float cellSize;
};


#endif //PROJECT_MACGRID_H
