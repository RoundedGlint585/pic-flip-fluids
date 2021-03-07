//
// Created by roundedglint585 on 3/6/21.
//

#ifndef PROJECT_MACGRID_H
#define PROJECT_MACGRID_H


#include <cstddef>
#include "vcl/vcl.hpp"

struct barycentricCoordinate{
    size_t index;
    float dist;
};
enum cellType{
    EMPTY_CELL,
    FLUID_CELL, // cell with at least one particle
};
class MACGrid {
public:
    MACGrid(size_t xCellNumber, size_t yCellNumber, float cellSize);
    size_t getXCellNumber() const;
    size_t getYCellNumber() const;
    float getCellSize() const;
    vcl::grid_2D<float>& getStaggeredHorizontal();
    vcl::grid_2D<float>& getStaggeredVertical();
    vcl::grid_2D<cellType> &getCellTypes();

    barycentricCoordinate barycentricOnAxis(float coordinate) const;
    barycentricCoordinate barycentricOffsetedX(float x) const; // for Mac grid
    barycentricCoordinate barycentricOffsetedY(float y) const; // TODO: sort of similar code, take a look on unification
    void updateDistanceField();
private:
    vcl::grid_2D<float> pressure, density; //store in the center
    vcl::grid_2D<float> staggeredHorizontal, staggeredVertical; //store on the edges
    vcl::grid_2D<cellType> cellTypes;
    vcl::grid_2D<float> distanceField; // phi field from paper created using fast sweep
    size_t xCellNumber, yCellNumber;
    float cellSize;
};


#endif //PROJECT_MACGRID_H
