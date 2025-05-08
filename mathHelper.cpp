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
    //TODO depreciate this function in favor of the Vec3D version
    double sqare_sum = (vec1[0]-vec2[0])*(vec1[0]-vec2[0]) + ((vec1[1]-vec2[1])*(vec1[1]-vec2[1])) + (((vec1[2]-vec2[2])*(vec1[2]-vec2[2]))); 
    return sqrt(sqare_sum);
};

double vector_difference(const Vec3D& vec1, const Vec3D& vec2) {
    double sqare_sum = (vec1.x-vec2.x)*(vec1.x-vec2.x) + ((vec1.y-vec2.y)*(vec1.y-vec2.y)) + (((vec1.z-vec2.z)*(vec1.z-vec2.z))); 
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
    const std::vector<std::vector<double>>& matrix) //uses legacy vector for arbitrary matrix size 
{
    // Multiplies an nx3 Matrix with a Vec3D
    // The result will be a Vec3D

    // Initialize the result vector
    Vec3D result(0, 0, 0);

    // Compute the dot product of each row of the matrix with the vector
    result.x = matrix[0][0] * vector1.x + matrix[0][1] * vector1.y + matrix[0][2] * vector1.z;
    result.y = matrix[1][0] * vector1.x + matrix[1][1] * vector1.y + matrix[1][2] * vector1.z;
    result.z = matrix[2][0] * vector1.x + matrix[2][1] * vector1.y + matrix[2][2] * vector1.z;

    return result; // Return the resultant vector
};

double clamp(double value, double min, double max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

double linearInterpolate(double a, double b, double t) {
    return (1.0 - t) * a + t * b;
}
int linearInterpolate(int a, int b, double t) {
    return static_cast<int>((1.0 - t) * a + t * b);
};

Vec3D vector_interpolate(const Vec3D& v1, const Vec3D& v2, double t, int flag) {
    return {
        linearInterpolate(v1.x,v2.x,t),
        linearInterpolate(v1.y,v2.y,t),
        linearInterpolate(v1.z,v2.z,t)
    };
};
//integer overload
iVec3D vector_interpolate(const iVec3D& v1, const iVec3D& v2, double t, int flag) {
    return {
        static_cast<int>(linearInterpolate(v1.x,v2.x,t)),
        static_cast<int>(linearInterpolate(v1.y,v2.y,t)),
        static_cast<int>(linearInterpolate(v1.z,v2.z,t))
    };
};