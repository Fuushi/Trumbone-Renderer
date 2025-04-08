#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>  // Required for std::vector
#include "baseClasses.h"

using namespace std;

//print vec (debug functions)
void print_v(std::vector<double> vec);
void print_v(std::vector<float> vec);
void print_v(std::vector<int> vec);

//simple functions
bool pin_interpolator_range_checker(double t, int flag);

//Complex vector expressions
//returns the angle between 2 vectors in deg
double get_vectors_angle(const std::vector<double>& ray1, const std::vector<double>& ray2);

//reflects a vector along a normal
std::vector<double> reflect_vector_normal(std::vector<double> vec, std::vector<double> normal);

//applied rotation matrix based on degrees from centre
Vec3D calculate_ray_heading(
    const Vec3D& input_vector,
    const double& fov,
    const std::vector<double>& uv
);

#endif // FUNCTIONS_H
