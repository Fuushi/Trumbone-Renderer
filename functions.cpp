#include "functions.h"
#include <cmath>  // Required for sqrt

std::vector<int> convert_vec_double_to_int(std::vector<double> vec) {
    return {static_cast<int>(vec[0]), static_cast<int>(vec[1]), static_cast<int>(vec[2])};
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
