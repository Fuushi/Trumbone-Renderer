#include "functions.h"
#include <cmath>  // Required for sqrt
#include <iostream>

#include <iostream>
#include <vector>

void print_v(const std::vector<double> vec) {
    std::cout << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;
}

void print_v(const std::vector<float> vec) {
    print_v(std::vector<double>(vec.begin(), vec.end()));  // Convert properly
}

void print_v(const std::vector<int> vec) {
    print_v(std::vector<double>(vec.begin(), vec.end()));  // Convert properly
}

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


std::vector<double> rotation_matrix_degrees(std::vector<double> input_vector, double fov, std::vector<double> uv) {
    //applies rotation matrix to a vector according to degrees, uv
    double theta_z_deg = (uv[0]*fov) - (fov / 2);
    double theta_x_deg = (uv[1]*fov) - (fov / 2);

    //convert to radians
    double theta_z_rad = theta_z_deg * (3.14159265358979323846264338327950288419716939937510582097494459230781640628 / 180);
    double theta_x_rad = theta_x_deg * (3.14159265358979323846264338327950288419716939937510582097494459230781640628 / 180);

    //define z-axis rotation matrix
    std::vector<std::vector<double>> r_z = {
        {cos(theta_z_rad), 0.0-sin(theta_z_rad), 0.0},
        {sin(theta_z_rad), cos(theta_z_rad), 0.0},
        {0.0, 0.0, 1.0}
    };

    //define x-axis rotation matrix
    std::vector<std::vector<double>> r_x = {
        {1.0, 0.0, 0.0},
        {0.0, cos(theta_x_rad), 0.0-sin(theta_x_rad)},
        {0.0, sin(theta_x_rad), cos(theta_x_rad)}
    };

    //apply matrix(s) via dot product
    //...
    std::vector<double> z_corrected_matrix = matrix_vector_multiplication(input_vector, r_z);

    std::vector<double> x_corrected_matrix = matrix_vector_multiplication(z_corrected_matrix, r_x);

    //normalize
    std::vector<double> normalized_vector = normalize_vector(x_corrected_matrix);
    
    return normalized_vector;
};


std::vector<double> reflect_vector_normal(std::vector<double> vec, std::vector<double> normal) {
    //normalize the vector
    std::vector<double> normalized_vec = normalize_vector(vec);

    //get the projection?
    std::vector<double> projection = scalar_multiply(normal, vector_dot_product(vec, normal));

    //calculate reflection
    std::vector<double> reflection = vector_subtract(vec, scalar_multiply(projection, 2.0));

    return reflection;
}


