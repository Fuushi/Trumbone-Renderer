#include "rays.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;


std::vector<int> principled_bdsf(Ray ray) {
    //
    if (!ray.intersect) {
        //no intersect, return sky
        return {50,20,20};
    }

    // no next reflection (maxed out), return material color
    if (ray.reflection_color.empty()) {
        return {200,200,200};
    }

    //here we can assume a reflection
    //blend between material color (depth) and reflect by frensel
    //TODO factor in frensel (just averaged for now)

    //get local color
    if (ray.intersect) {
        ray.color={200,200,200};
    } else {
        ray.color = {20,20,20};
    }

    //new color vector
    std::vector<int> color = {0,0,0};
    color[0] = (ray.color[0] + ray.reflection_color[0]) / 2;
    color[1] = (ray.color[1] + ray.reflection_color[1]) / 2;
    color[2] = (ray.color[2] + ray.reflection_color[2]) / 2;

    //return new color
    return color;
}
