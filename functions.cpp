#include "functions.h"
#include "mathHelper.h"
#include "baseClasses.h"
#include <cmath>  // Required for sqrt
#include <iostream>

#include <iostream>
#include <vector>

//print vec & overloads
void print_v(const std::vector<double> vec) {std::cout << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;}
void print_v(const std::vector<float> vec) {print_v(std::vector<double>(vec.begin(), vec.end()));}
void print_v(const std::vector<int> vec) {print_v(std::vector<double>(vec.begin(), vec.end()));}

//simple functions
bool pin_interpolator_range_checker(double t, int flag) {
    if (t < 0) {
        return false;
    }
    if (t > 1) {
        return false;
    }

    return true;
};

//Orthogonol vector definition and methods
struct Orthogonol {
    Vec3D Forward;
    Vec3D Right;
    Vec3D Up;
};
//gets orthogonal vector from euler with assumptions
Orthogonol get_orthogonol(const Vec3D& vec) {
    //get the orthoganol vector from a base vector making assumptions

    // rotate 90 degrees to the right
    Vec3D right(vec.y, -vec.x, 0);
    right++; //normalize

    //up is the cross product of the forward and right vectors
    Vec3D up(vec.cross(right));
    up++; //normalize

    //convert to legacy vectors before packing to orthogonol
    Orthogonol orthogonol;
    orthogonol.Forward=vec.to_double();
    orthogonol.Right=right.to_double();
    orthogonol.Up=up.to_double();
    return orthogonol;
};

//Complex vector operations
//gets the angle between 2 vectors in degrees
double get_vectors_angle(const std::vector<double>& ray1, const std::vector<double>& ray2) {

    // Compute the dot product of the vectors
    double dot_product = vector_dot_product(ray1, ray2);

    // Compute the cosine of the angle using the dot product formula
    double cos_theta = dot_product / (get_magnitude(ray1) * get_magnitude(ray2));

    // Clamp cos_theta to the range [-1, 1] to handle numerical precision issues
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    // Compute the angle in radians
    double angle_radians = std::acos(cos_theta);

    // Convert the angle to degrees
    double angle_degrees = rad_to_deg(angle_radians);

    return angle_degrees;
}

Vec3D rotate_vector_z(const Vec3D& vec, const double& rotation_rad) {
    
    // Define z-axis rotation matrix
    std::vector<std::vector<double>> r_z = {
        {cos(rotation_rad), -sin(rotation_rad), 0.0},
        {sin(rotation_rad), cos(rotation_rad), 0.0},
        {0.0, 0.0, 1.0}
    };

    // Apply z-axis rotation and return
    Vec3D rotated_vec = matrix_vector_multiplication(vec, r_z);
    rotated_vec++; //normalize
    return rotated_vec;
};

double get_z_rotation_rad(Vec3D vec) { //takes copy so i dont have to declare a cache explicitely
    //declare fwd vector (assumption)
    Vec3D Fwd = {1,0,0};

    //eliminate the z and normalize
    vec.z=0;
    vec++; //++ normalizes the vector

    //calculate angle
    double angle_radians = std::atan2(vec.y, vec.x);  // atan2 gives angle from X+ axis

    return angle_radians;
};

// Function to calculate the up/down angle in radians
double getUpDownAngleRadians(const Orthogonol& ortho) {
    // Extract Y and Z components of the Up vector
    double Ay = ortho.Up.x;
    double Az = ortho.Up.z;

    // Calculate angle in radians using atan2 for correct quadrant
    double angleRadians = std::atan2(Ay, -Az);

    return angleRadians;
}

// Rodrigues' rotation formula: Rotates `v` around axis `u` by `theta` radians
std::vector<double> rodriques_formaula(const Vec3D& v, const Vec3D& rotation_axis, double rad) {
    //calculate offsets
    double cos_theta = std::cos(rad);
    double sin_theta = std::sin(rad);
    double dot = v.dot(rotation_axis);
    
    //rotation matrix
    return {
        v.x * cos_theta + (rotation_axis.y * v.z - rotation_axis.z * v.y) * sin_theta + rotation_axis.x * dot * (1 - cos_theta),
        v.y * cos_theta + (rotation_axis.z * v.x - rotation_axis.x * v.z) * sin_theta + rotation_axis.y * dot * (1 - cos_theta),
        v.z * cos_theta + (rotation_axis.x * v.y - rotation_axis.y * v.x) * sin_theta + rotation_axis.z * dot * (1 - cos_theta)
    };
}

std::vector<double> reflect_vector_normal(std::vector<double> vec, std::vector<double> normal) {
    //normalize the vector
    std::vector<double> normalized_vec = normalize_vector(vec);

    //get the projection?
    std::vector<double> projection = scalar_multiply(normal, vector_dot_product(vec, normal));

    //calculate reflection
    std::vector<double> reflection = vector_subtract(vec, scalar_multiply(projection, 2.0));

    return reflection;
};

Vec3D calculate_ray_heading(
    const Vec3D& input_vector,
    const double& fov,
    const std::vector<double>& uv) 
    {
    // Calculate the heading of the ray in world space given
    //  -camera vector (world space)
    //  -uv (converted to camera normalized device coordinates)
    //  -fov (ignored for now)

    //convert uv2 to NDC (normalized device coordinates)
    // changes range from [0,-1] to [-1,1]
    std::vector<double> ndc = {
        (2*uv[1]) - 1,
        (2*(1-uv[0])) - 1
    };

    //calculate ray heading in camera space (cam space is looking to the x+)
    Vec3D cam_space(
        1, ndc[0], ndc[1]
    );
    cam_space++; //normalize

    //convert cam space to world space
    // establish rad deltas
    double rad_z_delta = get_z_rotation_rad(input_vector);
    
    // rotate along the Z
    Vec3D z_rotated_matrix = rotate_vector_z(cam_space, rad_z_delta);

    //get orthogonol vector
    Orthogonol orth = get_orthogonol(input_vector);
    
    //use ortho to get rotation angle
    double vertical_angle_rad = getUpDownAngleRadians(orth);
    
    //rotate along the orthoganal right of the camera
    std::vector<double> u_corrected_matrix = rodriques_formaula(z_rotated_matrix, orth.Right, vertical_angle_rad);

    //convert to Vec3D in place while working above
    Vec3D u_corrected_vec(u_corrected_matrix[0], u_corrected_matrix[1], u_corrected_matrix[2]);
    
    //normalize
    u_corrected_vec++;

    //return world space
    return u_corrected_vec;
};