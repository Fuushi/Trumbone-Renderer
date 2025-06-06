#include "rays.h"
#include <cmath>
#include <vector>
#include <iostream>

#include "functions.h"
#include "mathHelper.h"
#include "shaders.h"

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

//shader wrapper
std::vector<int> shader_wrapper(Ray ray, Lux lux, World world, ShaderInputs shader_inputs) {
    //run principled bdsf shader
    if (ray.intersect) {
        //run principled bdsf shader
        return geometry_shader(ray, lux, world, shader_inputs);
    } else {
        //run sky shader
        return sky_shader(ray);
    };
};

std::vector<int> sky_shader(Ray ray) {
    // Wrapper for the sky shader

    //call gradient sky shader from shaders.h
    std::vector<int> color = gradient_sky_shader(ray, {255,229,138}, {31,33,77}).to_vec(); //horizon color, zenith color
    
    //print_v(color);

    //return color
    return color; 

}; //indev

std::vector<int> geometry_shader(Ray ray, Lux lux, World world, ShaderInputs shader_inputs) {
    // Wrapper for the geometry shader (rename from principled_bdsf when moving logic to shaders.cpp)

    //drop inturrupts here...
    //return frensel_shader(ray).to_vec(); //frensel shader inturrupt

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
        Vec3D color(set_magnitude(convert_vec_int_to_double(lux.lighting_contributions[i].color.to_vec()), 1.0));
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

    const iVec3D sky_color = sky_shader(ray); //20,25,20
    if (!ray.intersect) {
        //no intersect, return sky
        return sky_color.to_vec(); //return sky color
    }

    // no next reflection (maxed out), return local color
    if (ray.reflection_color.empty) {
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
        ray.color = sky_color.to_vec();
    }

    //new color vector
    std::vector<int> color = {0,0,0};

    //mix colors by frensel fac
    color = mix_color(ray.reflection_color.to_vec(), ray.color.to_vec(), limit_range(ray.frensel/180+shader_inputs.gloss_diffuse_mix_fac, 0, 1));

    //return new color
    return color;
}
