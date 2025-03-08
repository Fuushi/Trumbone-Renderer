#include "functions.h"
#include <cmath>  // Required for sqrt
#include <iostream>

#include <iostream>
#include <vector>

struct Orthogonol {
    std::vector<double> Forward;
    std::vector<double> Right;
    std::vector<double> Up;
};

void print_v(const std::vector<double> vec) {
    std::cout << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
}

void print_v(const std::vector<float> vec) {
    print_v(std::vector<double>(vec.begin(), vec.end()));  // Convert properly
}

void print_v(const std::vector<int> vec) {
    print_v(std::vector<double>(vec.begin(), vec.end()));  // Convert properly
}

double deg_to_rad(double deg) {
    return deg * (3.1415692 / 180);
};

double rad_to_deg(double rad) {
    return rad * (180 / 3.141592);
};

std::vector<int> convert_vec_double_to_int(std::vector<double> vec) {
    return {static_cast<int>(vec[0]), static_cast<int>(vec[1]), static_cast<int>(vec[2])};
};
std::vector<double> convert_vec_int_to_double(std::vector<int> vec) {
    return {static_cast<double>(vec[0]), static_cast<double>(vec[1]), static_cast<double>(vec[2])};
};

double get_magnitude(std::vector<double> vec) {
    return sqrt((vec[0] * vec[0]) + (vec[1] * vec[1]) + (vec[2] * vec[2]));
};

std::vector<double> set_magnitude(const std::vector<double>& vec, double target_magnitude) {
    // Compute the magnitude of the input vector
    double magnitude = get_magnitude(vec);

    // Create a new vector and normalize its elements
    std::vector<double> result = vec;
    result[0] = (result[0] / magnitude) * target_magnitude;
    result[1] = (result[1] / magnitude) * target_magnitude;
    result[2] = (result[2] / magnitude) * target_magnitude;

    return result;
};

double vector_difference(std::vector<double> vec1, std::vector<double> vec2) {
    //make pythagoras proud
    double sqare_sum = (vec1[0]-vec2[0])*(vec1[0]-vec2[0]) + ((vec1[1]-vec2[1])*(vec1[1]-vec2[1])) + (((vec1[2]-vec2[2])*(vec1[2]-vec2[2]))); 
    return sqrt(sqare_sum);
};

std::vector<double> normalize_vector(const std::vector<double>& vec) {
    // Calculate the magnitude of the vector
    double magnitude = 0.0;
    for (double val : vec) {
        magnitude += val * val;
    }
    magnitude = std::sqrt(magnitude);

    // Check for zero magnitude to avoid division by zero
    //we clown in this bitch (not handled to avoid compiler error)

    // Create the normalized vector
    std::vector<double> normalized_vec(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        normalized_vec[i] = vec[i] / magnitude;
    }

    return normalized_vec;
};

std::vector<double> vector_cross_product(std::vector<double> vec_a, std::vector<double> vec_b) {
    std::vector<double> cross_product_vector = {
        vec_a[1]*vec_b[2] - vec_a[2]*vec_b[1],
        vec_a[2]*vec_b[0] - vec_a[0]*vec_b[2],
        vec_a[0]*vec_b[1] - vec_a[1]*vec_b[0]
    };
    return cross_product_vector;
};


double vector_dot_product(std::vector<double> vec_a, std::vector<double> vec_b) {
    double product = 0;

    for (int i = 0; i < vec_a.size(); i++) {
        product = product + vec_a[i]*vec_b[i];
    };
    return product;
};


std::vector<double> matrix_subration(std::vector<double> vec_a, std::vector<double> vec_b) {
    std::vector<double> result(vec_a.size());

    for (int i = 0; i < vec_a.size(); i++) {
        result[i] = vec_a[i] - vec_b[i];
    };

    return result;
};

std::vector<double> scalar_multiply(const std::vector<double> vec, double scalar) {
    return {vec[0] * scalar, vec[1] * scalar, vec[2] * scalar};
};

std::vector<double> vector_add(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    return {vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]};
};

std::vector<double> vector_subtract(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    return {vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]};
};

std::vector<double> vector_multiply(std::vector<double> vec1, std::vector<double> vec2) {
    return {
        vec1[0]*vec2[0], vec1[1]*vec2[1], vec1[2]*vec2[2]
    };
}

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


double get_vectors_angle(std::vector<double> ray1, std::vector<double> ray2) {

    // Compute the dot product of the vectors
    double dot_product = ray1[0] * ray2[0] + ray1[1] * ray2[1] + ray1[2] * ray2[2];

    // Compute the magnitudes of the vectors
    double magnitude1 = std::sqrt(ray1[0] * ray1[0] + ray1[1] * ray1[1] + ray1[2] * ray1[2]);
    double magnitude2 = std::sqrt(ray2[0] * ray2[0] + ray2[1] * ray2[1] + ray2[2] * ray2[2]);

    // Compute the cosine of the angle using the dot product formula
    double cos_theta = dot_product / (magnitude1 * magnitude2);

    // Clamp cos_theta to the range [-1, 1] to handle numerical precision issues
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

    // Compute the angle in radians
    double angle_radians = std::acos(cos_theta);

    // Convert the angle to degrees
    double angle_degrees = angle_radians * (180.0 / 3.1415926535);

    return angle_degrees;
}

std::vector<double> matrix_vector_multiplication(
    std::vector<double> vector1,
    const std::vector<std::vector<double>>& matrix) 
{
    // Ensure dimensions match for multiplication
    // again, not handled to avoid a compiler error

    // Resultant vector of size equal to the number of rows in the matrix
    std::vector<double> result(matrix.size(), 0.0);

    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            result[i] += matrix[i][j] * vector1[j];
        }
    }

    return result;
};

std::vector<double> rotate_vector_z(std::vector<double> vec, double rotation_rad) {
    
    // Define z-axis rotation matrix
    std::vector<std::vector<double>> r_z = {
        {cos(rotation_rad), -sin(rotation_rad), 0.0},
        {sin(rotation_rad), cos(rotation_rad), 0.0},
        {0.0, 0.0, 1.0}
    };

    // Apply z-axis rotation first
    std::vector<double> z_corrected_matrix = normalize_vector(matrix_vector_multiplication(vec, r_z));

    return z_corrected_matrix;
};

double get_z_rotation_rad(std::vector<double> vec) {
    //declare fwd
    std::vector<double> Fwd = {1,0,0};

    //eliminate the z and normalize
    vec[2]=0;
    vec = normalize_vector(vec);

    //calculate angle
    double angle_radians = std::atan2(vec[1], vec[0]);  // atan2 gives angle from X+ axis

    return angle_radians;
};

//get vertical angle along Right
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
    
    double cos_theta = std::cos(rad);
    double sin_theta = std::sin(rad);
    double dot = v[0] * rotation_axis[0] + v[1] * rotation_axis[1] + v[2] * rotation_axis[2];

    return {
        v[0] * cos_theta + (rotation_axis[1] * v[2] - rotation_axis[2] * v[1]) * sin_theta + rotation_axis[0] * dot * (1 - cos_theta),
        v[1] * cos_theta + (rotation_axis[2] * v[0] - rotation_axis[0] * v[2]) * sin_theta + rotation_axis[1] * dot * (1 - cos_theta),
        v[2] * cos_theta + (rotation_axis[0] * v[1] - rotation_axis[1] * v[0]) * sin_theta + rotation_axis[2] * dot * (1 - cos_theta)
    };
}

std::vector<double> calculate_ray_heading(std::vector<double> input_vector, double fov, std::vector<double> uv) {
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
    //...

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

std::vector<double> reflect_vector_normal(std::vector<double> vec, std::vector<double> normal) {
    //normalize the vector
    std::vector<double> normalized_vec = normalize_vector(vec);

    //get the projection?
    std::vector<double> projection = scalar_multiply(normal, vector_dot_product(vec, normal));

    //calculate reflection
    std::vector<double> reflection = vector_subtract(vec, scalar_multiply(projection, 2.0));

    return reflection;
};

double linearInterpolate(double a, double b, double t) {
    return (1.0 - t) * a + t * b;
}

std::vector<double> vector_interpolate(std::vector<double> v1, std::vector<double> v2, double t, int flag) {
    return {
        linearInterpolate(v1[0],v2[0],t),
        linearInterpolate(v1[1],v2[1],t),
        linearInterpolate(v1[2],v2[2],t)
    };
};

bool pin_interpolator_range_checker(double t, int flag) {
    if (t < 0) {
        return false;
    }
    if (t > 1) {
        return false;
    }

    return true;
};
