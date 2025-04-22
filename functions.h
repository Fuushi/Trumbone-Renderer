#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>  // Required for std::vector
#include "baseClasses.h"

using namespace std;

//print vec (debug functions)
void print_v(std::vector<double> vec);
void print_v(std::vector<float> vec);
void print_v(std::vector<int> vec);

//vec3D overloads
void print_v(const Vec3D& vec);
void print_v(const iVec3D& vec);
void print_v(const Vec2D& vec);

//Orthogonol vector definition and methods
struct Orthogonol {
    Vec3D Forward;
    Vec3D Right;
    Vec3D Up;
};

Orthogonol get_orthogonol(const Vec3D& vec);

Vec3D rotate_vector_z(const Vec3D& vec, const double& rotation_rad);

double getUpDownAngleRadians(const Orthogonol& ortho);

Vec3D rodriques_formaula(const Vec3D& v, const Vec3D& rotation_axis, double rad);

//simple functions
bool pin_interpolator_range_checker(double t, int flag);

//Complex vector expressions
//returns the angle between 2 vectors in deg
double get_vectors_angle(const Vec3D& ray1, const Vec3D& ray2);

//reflects a vector along a normal
Vec3D reflect_vector_normal(Vec3D vec, Vec3D normal);

//applied rotation matrix based on degrees from centre
Vec3D calculate_ray_heading(
    const Vec3D& input_vector,
    const double& fov,
    const std::vector<double>& uv
);

#endif // FUNCTIONS_H
