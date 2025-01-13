#include "rays.h"
#include <cmath>
#include <vector>


std::vector<int> principled_bdsf(Ray ray) {
    //
    if (!ray.intersect) {return {20,20,20};}


    return {255,255,255};
}
