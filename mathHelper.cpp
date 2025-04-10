#include "mathHelper.h"
#include "baseClasses.h"
#include <vector>
#include <cmath>

std::vector<int> convert_vec_double_to_int(std::vector<double> vec) {
    return {static_cast<int>(vec[0]), static_cast<int>(vec[1]), static_cast<int>(vec[2])};
};
std::vector<double> convert_vec_int_to_double(std::vector<int> vec) {
    return {static_cast<double>(vec[0]), static_cast<double>(vec[1]), static_cast<double>(vec[2])};
};

double deg_to_rad(const double& deg) {
    return deg * (3.1415692 / 180);
};

double rad_to_deg(const double& rad) {
    return rad * (180 / 3.141592);
};

double get_magnitude(const std::vector<double>& vec) {
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

std::vector<double> normalize_vector(const std::vector<double>& vec) {
    return set_magnitude(vec, 1.0);
};

std::vector<double> scalar_multiply(const std::vector<double>& vec, const double& scalar) {
    return {vec[0] * scalar, vec[1] * scalar, vec[2] * scalar};
};

std::vector<double> vector_add(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    return {vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]};
};

std::vector<double> vector_subtract(const std::vector<double>& vec_a, const std::vector<double>& vec_b) {
    return vector_add(vec_a, scalar_multiply(vec_b, -1.0));
};

std::vector<double> vector_multiply(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    return {vec1[0]*vec2[0], vec1[1]*vec2[1], vec1[2]*vec2[2]};
}

double vector_difference(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    double sqare_sum = (vec1[0]-vec2[0])*(vec1[0]-vec2[0]) + ((vec1[1]-vec2[1])*(vec1[1]-vec2[1])) + (((vec1[2]-vec2[2])*(vec1[2]-vec2[2]))); 
    return sqrt(sqare_sum);
};

std::vector<double> vector_cross_product(const std::vector<double>& vec_a, const std::vector<double>& vec_b) {
    std::vector<double> cross_product_vector = {
        vec_a[1]*vec_b[2] - vec_a[2]*vec_b[1],
        vec_a[2]*vec_b[0] - vec_a[0]*vec_b[2],
        vec_a[0]*vec_b[1] - vec_a[1]*vec_b[0]
    };
    return cross_product_vector;
};

double vector_dot_product(const std::vector<double>& vec_a, const std::vector<double>& vec_b) {
    double product = 0;

    for (int i = 0; i < vec_a.size(); i++) {
        product = product + vec_a[i]*vec_b[i];
    };
    return product;
};

Vec3D matrix_vector_multiplication(
    const Vec3D& vector1,
    const std::vector<std::vector<double>>& matrix) //uses legacy vector for arbirary matrix size 
{
    //Multiplies an nx3 Matrix with a Vec3D
    // The result will be a Vec3D

    // Initialize the result vector
    Vec3D result(0,0,0);

    // Iterate through matrix rows (n)
    for (int i = 0; i < matrix.size(); ++i) {
        // matrix MUST be nx3

        result.x += matrix[i][0] * vector1.x;
        result.y += matrix[i][1] * vector1.y;
        result.z += matrix[i][2] * vector1.z;

    }

    return result; //return the resultant vector
};

double linearInterpolate(double a, double b, double t) {
    return (1.0 - t) * a + t * b;
}

std::vector<double> vector_interpolate(const std::vector<double>& v1, const std::vector<double>& v2, double t, int flag) {
    return {
        linearInterpolate(v1[0],v2[0],t),
        linearInterpolate(v1[1],v2[1],t),
        linearInterpolate(v1[2],v2[2],t)
    };
};