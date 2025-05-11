
//includes
#include "shaders.h"

//namespace (debug)
#include <iostream>

//define shader functions

iVec3D gradient_sky_shader(const Ray& ray, iVec3D horizon_color, iVec3D zenith_color) {

    //calculate gradient factor
    double fac = abs(clamp(ray.euler.z, -1.0, 1.0));

    //print_v(ray.euler.to_double());

    //calculate gradient color
    iVec3D gradient_color = vector_interpolate(horizon_color, zenith_color, fac, 0);

    //return color
    return gradient_color;
};

iVec3D frensel_shader(const Ray& ray) {
    //inherit frensel factor
    double frensel_factor = ray.frensel;

    //calculate frensel inline and compare for debug
    

    //pack to iVec3D and return
    return {
        static_cast<int>(frensel_factor),
        static_cast<int>(frensel_factor),
        static_cast<int>(frensel_factor)
    };
};