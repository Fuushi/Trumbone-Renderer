#ifndef BASECLASSES_H
#define BASECLASSES_H

#include <math.h> //ensure no double include issues

struct Vec3D {
    double x, y, z;
    
    //basic operators

    //addition
    Vec3D operator+(const Vec3D& other) const {
        return Vec3D(x + other.x, y + other.y, z + other.z);
    }

    //subtraction
    Vec3D operator-(const Vec3D& other) const {
        return Vec3D(x - other.x, y - other.y, z - other.z);
    }

    //scalar multiplication (to divide use the inverse)
    Vec3D operator%(double scalar) const {
        return Vec3D(x * scalar, y * scalar, z * scalar);
    }    

    //element wise multiply
    Vec3D operator*(const Vec3D& other) const {
        return Vec3D(x * other.x, y * other.y, z * other.z);
    }

    //element wise divide
    Vec3D operator/(const Vec3D& other) const {
        return Vec3D(
            other.x != 0 ? x / other.x : 0,
            other.y != 0 ? y / other.y : 0,
            other.z != 0 ? z / other.z : 0
        );
    }

    //normalize self (++)
    Vec3D& operator++() {
        double magnitude = sqrt(x * x + y * y + z * z);
        if (magnitude > 0) {
            x /= magnitude;
            y /= magnitude;
            z /= magnitude;
        }
        return *this;
    }

    //assignment operator
    Vec3D& operator=(const Vec3D& other) {
        if (this != &other) { // self-assignment check
            x = other.x;
            y = other.y;
            z = other.z;
        }
        return *this;
    }


    //get methods

    //magnitude
    double magnitude() const {
        return sqrt(x * x + y * y + z * z);
    }

    //compute methods

    //dot product
    double dot(const Vec3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    //cross product
    Vec3D cross(const Vec3D& other) const {
        return Vec3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
    
    //
    
    //normalize vector
    Vec3D normalize() const {
        double magnitude = sqrt(x * x + y * y + z * z);
        if (magnitude > 0) {
            return Vec3D(x / magnitude, y / magnitude, z / magnitude);
        }
        return *this; // Return the original vector if magnitude is zero (edge case)
    }

    //set magnitude (like normalize but to a specific value)
    Vec3D set_magnitude(double target_magnitude) const {
        double current_magnitude = magnitude();
        if (current_magnitude > 0) {
            return Vec3D(x * target_magnitude / current_magnitude, y * target_magnitude / current_magnitude, z * target_magnitude / current_magnitude);
        }
        return *this; // Return the original vector if magnitude is zero (edge case)
    }
    
    //....

    // constructor
    Vec3D() : x(0), y(0), z(0) {}
    Vec3D(double x, double y, double z) : x(x), y(y), z(z) {}
};


#endif