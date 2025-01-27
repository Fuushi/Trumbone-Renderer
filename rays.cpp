#include "rays.h"
#include "functions.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;


std::vector<int> mix_color(std::vector<int> c1, std::vector<int> c2, double fac) {
    //indev (will raise a warning until implemented)
};

std::vector<int> ambient_occlusion(Ray ray, World world) {
    //fast, but inaccurate lighting effect
    std::vector<int> sum = {0,0,0};
    for (int i = 0; i < world.lights.size(); i++) {
        double angle = get_vectors_angle(ray.surface_normal, world.lights[i].vec);
        //^ tangent = 0, parallel = 90, away = 180

        sum[0] = sum[0] + static_cast<int>(angle);
        sum[1] = sum[1] + static_cast<int>(angle);
        sum[2] = sum[2] + static_cast<int>(angle);

    }

    sum[0] = sum[0] / world.lights.size();
    sum[1] = sum[1] / world.lights.size();
    sum[2] = sum[2] / world.lights.size();

    return sum;
};

std::vector<int> principled_bdsf(Ray ray, Lux lux, World world) {
    //

    //calculate lux sum
    double sum = 0;
    for (int i = 0; i < lux.lighting_contributions.size(); i++) {
        //dropoff function

        //e^-x
        //double influence = pow(2.71828182, 0-lux.lighting_contributions[i].depth);
        //because type is sun, influence is always 1, split logit with if
        double influence = 1;

        //multiply brightness by influence to get lux at point
        //inrelation to n light
        if (!lux.lighting_contributions[i].obstructed) {
            sum = sum + (lux.lighting_contributions[i].brightness*influence);
        }
        

    };

    int int_sum = static_cast<int>(sum);
    return {int_sum,int_sum,int_sum};


    // disabled for debug
    const std::vector<int> material_color = ambient_occlusion(ray, world);
    const std::vector<int> sky_color = {20,25,20};
    if (!ray.intersect) {
        //no intersect, return sky
        return sky_color;
    }

    // no next reflection (maxed out), return material color
    if (ray.reflection_color.empty()) {
        return material_color;
    }

    //here we can assume a reflection
    //blend between material color (depth) and reflect by frensel
    //TODO factor in frensel (just averaged for now)

    //get local color
    if (ray.intersect) {
        //material color
        ray.color=material_color;
    } else {
        //sky color
        ray.color = sky_color;
    }

    //new color vector
    std::vector<int> color = {0,0,0};
    color[0] = (ray.color[0] + ray.reflection_color[0]) / 2;
    color[1] = (ray.color[1] + ray.reflection_color[1]) / 2;
    color[2] = (ray.color[2] + ray.reflection_color[2]) / 2;

    //return new color
    return color;
}
