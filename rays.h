#ifndef SHADERS_H
#define SHADERS_H

#include <vector>

using namespace std;

struct Ray {
    //struct contains complex data related to ray intersections
    bool intersect = false; //not implemented
    std::vector<int> color = {20,20,20};
    std::vector<double> intersect_point = {0.0,0.0,0.0};
    std::vector<double> surface_normal = {0.0,0.0,0.0};
    std::vector<double> reflection_vec = {0.0,0.0,0.0};
    std::vector<int> reflection_color;
    double depth = 0.0;
    double frensel = 0.0;
};

std::vector<int> principled_bdsf(Ray ray);

#endif