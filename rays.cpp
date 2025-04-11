#include "rays.h"
#include <cmath>
#include <vector>
#include <iostream>

#include "functions.h"
#include "mathHelper.h"

using namespace std;

int linear_interpolation(int a, int b, double fac) {
    int num = static_cast<int>((1-fac) * static_cast<double>(a) + fac * static_cast<double>(b));
    return num;
};

double limit_range(double num, double min, double max) {
    num = std::max(num, min);
    num = std::min(num, max);
    return num;
}

std::vector<int> mix_color(std::vector<int> c1, std::vector<int> c2, double fac) {
    //indev (will raise a warning until implemented)
    return std::vector<int> {
        linear_interpolation(c1[0], c2[0], fac),
        linear_interpolation(c1[1], c2[1], fac),
        linear_interpolation(c1[2], c2[2], fac),
    };
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

std::vector<int> principled_bdsf(Ray ray, Lux lux, World world, ShaderInputs shader_inputs) {
    //

    //Lux sums by color channel
    Vec3D sums;
    for (int i = 0; i < lux.lighting_contributions.size(); i++) {
        //dropoff function

        //e^-x or 1
        //multiply by brightness to get Lux
        //for point source: pow(2.71828182, 0-lux.lighting_contributions[i].depth);
        //for sun: 1
        double influence = 1;


        //static cast color to double and normalize
        //multiply channel by Lux to get color Lux
        Vec3D color(set_magnitude(convert_vec_int_to_double(lux.lighting_contributions[i].color), 1.0));
        if (lux.lighting_contributions[i].obstructed) {
            continue;
        }
        
        // sum += color * brightness(lumine) * influence
        sums.x = sums.x + color.x * ((lux.lighting_contributions[i].brightness*influence));
        sums.y = sums.y + color.y * ((lux.lighting_contributions[i].brightness*influence));
        sums.z = sums.z + color.z * ((lux.lighting_contributions[i].brightness*influence));
        
        
    };
    //max function for alphas
    sums.x = min(sums.x, 255.0);
    sums.y = min(sums.y, 255.0);
    sums.z = min(sums.z, 255.0);

    //define material colors from inputs
    iVec3D material_color = shader_inputs.material_color;

    //local color after illuminance
    iVec3D local_color_legacy = convert_vec_double_to_int( //TODO rewrite legacy code
        vector_multiply(
            convert_vec_int_to_double(material_color.to_vec()),
            sums.to_double()
        )
    );

    // multiplies material color by illuminance for each channel in double precision
    // before converting back to int (i need a function for this)
    iVec3D local_color = convert_vec_double_to_int((material_color.to_double()*sums).to_double());

    const std::vector<int> sky_color = world.sky_color; //20,25,20
    if (!ray.intersect) {
        //no intersect, return sky
        return sky_color;
    }

    // no next reflection (maxed out), return local color
    if (ray.reflection_color.empty()) {
        return local_color.to_vec(); //return local color
    }

    //here we can assume a reflection
    //blend between material color (depth) and reflect by frensel

    //get local color
    if (ray.intersect) {
        //material color
        ray.color=local_color.to_vec();
    } else {
        //sky color
        ray.color = sky_color;
    }

    //new color vector
    std::vector<int> color = {0,0,0};

    //mix colors by frensel fac
    color = mix_color(ray.reflection_color, ray.color, limit_range(ray.frensel/180+shader_inputs.gloss_diffuse_mix_fac, 0, 1));

    //return new color
    return color;
}
