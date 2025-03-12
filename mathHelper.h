#ifndef MATHHELPER_H
#define MATHHELPER_H

#include <vector>

//static casts a double vector to int vector & other
std::vector<int> convert_vec_double_to_int(std::vector<double> vec);
std::vector<double> convert_vec_int_to_double(std::vector<int> vec);

//converts degrees to radians & visca versa
double deg_to_rad(const double& deg);
double rad_to_deg(const double& rad);

//gets magnitude from a double vector
double get_magnitude(const std::vector<double>& vec);

//sets the magnitude of a double vector (set to 1.0 to normalize)
std::vector<double> set_magnitude(const std::vector<double>& vec, double target_magnitude);

//normalizes a vector by value
std::vector<double> normalize_vector(const std::vector<double>& vec);

//multiplies a double vec by a double scalar
std::vector<double> scalar_multiply(const std::vector<double>& vec, const double& scalar);

//adds 2 double vecs
std::vector<double> vector_add(const std::vector<double>& vec1, const std::vector<double>& vec2);

//subtracts double vec a from double vec b (todo add int overload)
std::vector<double> vector_subtract(const std::vector<double>& vec_a, const std::vector<double>& vec_b);

//multiplies 2 vecs
std::vector<double> vector_multiply(const std::vector<double>& vec1, const std::vector<double>& vec2);

//gets the distance between 2 double vectors -> double
double vector_difference(const std::vector<double>& vec1, const std::vector<double>& vec2);

//returns the cross product of 2 double vecs
std::vector<double> vector_cross_product(const std::vector<double>& vec_a, const std::vector<double>& vec_b);

//returns the dot product of 2 double vecs
double vector_dot_product(const std::vector<double>& vec_a, const std::vector<double>& vec_b);

//tbh i dont remember what this does but i know it's important
std::vector<double> matrix_vector_multiplication(
    const std::vector<double>& vector1,
    const std::vector<std::vector<double>>& matrix
); 

//animator functions
double linearInterpolate(double a, double b, double t); //self explanitory

std::vector<double> vector_interpolate(const std::vector<double>& v1, const std::vector<double>& v2, double t, int flag);


#endif