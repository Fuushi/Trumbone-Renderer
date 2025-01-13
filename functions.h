#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>  // Required for std::vector

using namespace std;

//static casts a double vector to int vector
std::vector<int> convert_vec_double_to_int(std::vector<double> vec);

//gets magnitude from a double vector
double get_magnitude(std::vector<double> vec);

//sets the magnitude of a double vector (set to 1.0 to normalize)
std::vector<double> set_magnitude(const std::vector<double>& vec, double target_magnitude);

//gets the difference between 2 double vectors -> double
double vector_difference(std::vector<double> vec1, std::vector<double> vec2);

//normalizes a vector by value
std::vector<double> normalize_vector(const std::vector<double>& vec);

//vector operations
//returns the cross product of 2 double vecs
std::vector<double> vector_cross_product(std::vector<double> vec_a, std::vector<double> vec_b);

//returns the dot product of 2 double vecs
double vector_dot_product(std::vector<double> vec_a, std::vector<double> vec_b);

//subtracts double vec a from double vec b (todo add int overload)
std::vector<double> matrix_subration(std::vector<double> vec_a, std::vector<double> vec_b);

//multiplies a double vec by a double scalar
std::vector<double> scalar_multiply(const std::vector<double> vec, double scalar);

//adds 2 double vecs
std::vector<double> vector_add(const std::vector<double>& vec1, const std::vector<double>& vec2);

//tbh i dont remember what this does but i know it's important
std::vector<double> matrix_vector_multiplication(
    std::vector<double> vector1,
    const std::vector<std::vector<double>>& matrix
); 
#endif // FUNCTIONS_H
