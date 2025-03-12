#include "functions.h"
#include "mathHelper.h"
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
    std::vector<double> Forward;
    std::vector<double> Right;
    std::vector<double> Up;
};
//gets orthogonal vector from euler with assumptions
Orthogonol get_orthogonol(std::vector<double> vec) {
    //get the orthoganol vector from a base vector making assumptions

    // rotate 90 degrees to the right
    std::vector<double> right = normalize_vector({vec[1], -vec[0], 0});

    //up is the cross product
    std::vector<double> up = vector_cross_product(vec, right);
     
    Orthogonol orthogonol;
    orthogonol.Forward=vec;
    orthogonol.Right=right;
    orthogonol.Up=up;
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

std::vector<double> rotate_vector_z(const std::vector<double>& vec, const double& rotation_rad) {
    
    // Define z-axis rotation matrix
    std::vector<std::vector<double>> r_z = {
        {cos(rotation_rad), -sin(rotation_rad), 0.0},
        {sin(rotation_rad), cos(rotation_rad), 0.0},
        {0.0, 0.0, 1.0}
    };

    // Apply z-axis rotation and return
    return normalize_vector(matrix_vector_multiplication(vec, r_z));
};

double get_z_rotation_rad(std::vector<double> vec) { //takes copy so i dont have to declare a cache explicitely
    //declare fwd vector (assumption)
    std::vector<double> Fwd = {1,0,0};

    //eliminate the z and normalize
    vec[2]=0;
    vec = normalize_vector(vec);

    //calculate angle
    double angle_radians = std::atan2(vec[1], vec[0]);  // atan2 gives angle from X+ axis

    return angle_radians;
};

// Function to calculate the up/down angle in radians
double getUpDownAngleRadians(const Orthogonol& ortho) {
    // Extract Y and Z components of the Up vector
    double Ay = ortho.Up[1];
    double Az = ortho.Up[2];

    // Calculate angle in radians using atan2 for correct quadrant
    double angleRadians = std::atan2(Ay, -Az);

    return angleRadians;
}

// Rodrigues' rotation formula: Rotates `v` around axis `u` by `theta` radians
std::vector<double> rodriques_formaula(const std::vector<double>& v, const std::vector<double>& rotation_axis, double rad) {
    //calculate offsets
    double cos_theta = std::cos(rad);
    double sin_theta = std::sin(rad);
    double dot = vector_dot_product(v, rotation_axis);
    
    //rotation matrix
    return {
        v[0] * cos_theta + (rotation_axis[1] * v[2] - rotation_axis[2] * v[1]) * sin_theta + rotation_axis[0] * dot * (1 - cos_theta),
        v[1] * cos_theta + (rotation_axis[2] * v[0] - rotation_axis[0] * v[2]) * sin_theta + rotation_axis[1] * dot * (1 - cos_theta),
        v[2] * cos_theta + (rotation_axis[0] * v[1] - rotation_axis[1] * v[0]) * sin_theta + rotation_axis[2] * dot * (1 - cos_theta)
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

std::vector<double> calculate_ray_heading(
    const std::vector<double>& input_vector, 
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
    std::vector<double> cam_space = normalize_vector(
        {
            1, ndc[0], ndc[1]
        }
    );

    //convert cam space to world space
    //establish rad deltas
    double rad_z_delta = get_z_rotation_rad(input_vector);
    
    // rotate along the Z
    std::vector<double> z_rotated_matrix = rotate_vector_z(cam_space, rad_z_delta);

    //get orthogonol vector
    Orthogonol orth = get_orthogonol(input_vector);
    
    //use ortho to get rotation angle
    double vertical_angle_rad = getUpDownAngleRadians(orth);
    
    //rotate along the orthoganal right of the camera
    std::vector<double> u_corrected_matrix = rodriques_formaula(z_rotated_matrix, orth.Right, vertical_angle_rad);

    //return world space
    return normalize_vector(u_corrected_matrix);
};
