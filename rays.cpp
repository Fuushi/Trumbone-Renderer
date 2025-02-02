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

    //Lux sums by color channel
    std::vector<double> sums = {0,0,0};
    for (int i = 0; i < lux.lighting_contributions.size(); i++) {
        //dropoff function

        //e^-x or 1
        //multiply by brightness to get Lux
        //for point source: pow(2.71828182, 0-lux.lighting_contributions[i].depth);
        //for sun: 1
        double influence = 1;

        

        //static cast color to double and normalize
        //multiply channel by Lux to get color Lux
        std::vector<double> color = set_magnitude(convert_vec_int_to_double(lux.lighting_contributions[i].color), 1.0);
        if (lux.lighting_contributions[i].obstructed) {
            continue;
        }
        
        // sum += color * brightness(lumine) * influence
        sums[0] = sums[0] + color[0] * ((lux.lighting_contributions[i].brightness*influence));
        sums[1] = sums[1] + color[1] * ((lux.lighting_contributions[i].brightness*influence));
        sums[2] = sums[2] + color[2] * ((lux.lighting_contributions[i].brightness*influence));
        
        
    };
    //max function for alphas
    for (int i=0; i<sums.size(); i++) {sums[i]=min(sums[i],255.0);};

    //define material colors
    std::vector<int> material_color = {128,128,128};

    //local color after illuminance
    std::vector<int> local_color = convert_vec_double_to_int(
        vector_multiply(
            convert_vec_int_to_double(material_color),
            sums
        )
    );

    const std::vector<int> sky_color = {80,90,80}; //20,25,20
    if (!ray.intersect) {
        //no intersect, return sky
        return sky_color;
    }

    // no next reflection (maxed out), return local color
    if (ray.reflection_color.empty()) {
        return local_color; //return local color
    }

    //here we can assume a reflection
    //blend between material color (depth) and reflect by frensel
    //TODO factor in frensel (just averaged for now)

    //get local color
    if (ray.intersect) {
        //material color
        ray.color=local_color;
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
