//
// Created by roundedglint585 on 3/6/21.
//

#ifndef PROJECT_MACGRID_H
#define PROJECT_MACGRID_H
#include <functional>
#include "vcl/vcl.hpp"

enum cellType{
    EMPTY_CELL,
    FLUID_CELL,
    SOLID_CELL,
};
struct barycentricCoords{
    size_t index;
    float offset;
};

struct MACGrid{
    size_t xCellCount{}, yCellCount{};
    float cellSize{};
    vcl::grid_2D<float> u, du;
    vcl::grid_2D<float> v, dv;
    vcl::grid_2D<cellType> cellTypes;
    vcl::grid_2D<float> distanceField;
    vcl::grid_2D<float> div;

    MACGrid() = default;

    MACGrid(int xCellCount, int yCellCount, float cellSize);

    barycentricCoords barycentricOnX(float x) const;
    barycentricCoords barycentricOnY(float y) const;

    void update(float dt);
    void updateBoundaries();
private:
    void saveFlipVelocities();
    void updateExternalForces(float dt);
    void updateDistanceField();
    void interpolateVelocities();
    void divFreeField();
    void updateVelocities();

    void sweepU();
    void sweepV();
    void calculateDiv();
    // this shit is tricky, i don't want to fix every part of fast sweeping, so this thing should do for in for and apply function
    // based on i, j, di, dj, basically just sent lambda
    void performSweep(int fromX, int toX, int fromY, int toY, const std::function<void(int, int, int, int)> &function);
};


#endif //PROJECT_MACGRID_H
