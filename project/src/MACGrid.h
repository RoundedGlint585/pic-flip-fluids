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

class MACGrid {
public:
    MACGrid(size_t xCellNumber, size_t yCellNumber, float cellSize);
    size_t getXCellNumber() const;
    size_t getYCellNumber() const;
    float getCellSize() const;

    vcl::grid_2D<float>& getStaggeredHorizontal();
    vcl::grid_2D<float>& getStaggeredVertical();

    barycentricCoordinate barycentricOnAxis(float coordinate) const;
    barycentricCoordinate barycentricOnX(float x) const; // for Mac grid
    barycentricCoordinate barycentricOnY(float y) const; // TODO: sort of similar code, take a look on unification

private:
    vcl::grid_2D<float> pressure, density; //store in the center
    vcl::grid_2D<float> staggeredHorizontal, staggeredVertical; //store on the edges
    size_t xCellNumber, yCellNumber;
    float cellSize;
};


#endif //PROJECT_MACGRID_H
