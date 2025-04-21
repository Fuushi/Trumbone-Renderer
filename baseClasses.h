#ifndef BASECLASSES_H
#define BASECLASSES_H

#include <math.h> //ensure no double include issues
#include <vector> //for legacy conversion (remove later)

struct Vec3D {
    double x, y, z; 
    bool empty = false; //used to check if vector is empty (0,0,0)
    
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
    Vec3D& operator++(int) {
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
    
    //...

    //convert to double
    std::vector<double> to_double() const {
        return {x, y, z};
    }

    // constructor
    Vec3D() : x(0), y(0), z(0) {}
    Vec3D(double x, double y, double z) : x(x), y(y), z(z) {}
    Vec3D(const std::vector<double>& vec) {
        if (vec.size() >= 3) {
            x = vec[0];
            y = vec[1];
            z = vec[2];
        } else {
            x = y = z = 0; // Default to zero if vector size is insufficient
        }
    }
};

struct Vec2D {
    double x, y;
    bool empty = false; //used to check if vector is empty (0,0)

    //basic operators

    //addition
    Vec2D operator+(const Vec2D& other) const {
        return Vec2D(x + other.x, y + other.y);
    }

    //subtraction
    Vec2D operator-(const Vec2D& other) const {
        return Vec2D(x - other.x, y - other.y);
    }

    //scalar multiplication (to divide use the inverse)
    Vec2D operator%(double scalar) const {
        return Vec2D(x * scalar, y * scalar);
    }

    //element wise multiply
    Vec2D operator*(const Vec2D& other) const {
        return Vec2D(x * other.x, y * other.y);
    }

    //element wise divide
    Vec2D operator/(const Vec2D& other) const {
        return Vec2D(
            other.x != 0 ? x / other.x : 0,
            other.y != 0 ? y / other.y : 0
        );
    }

    //normalize self (++)
    Vec2D& operator++(int) {
        double magnitude = sqrt(x * x + y * y);
        if (magnitude > 0) {
            x /= magnitude;
            y /= magnitude;
        }
        return *this;
    }

    //assignment operator
    Vec2D& operator=(const Vec2D& other) {
        if (this != &other) { // self-assignment check
            x = other.x;
            y = other.y;
        }
        return *this;
    }

    //get operators

    //magnitude
    double magnitude() const {
        return sqrt(x * x + y * y);
    }

    //compute operators
    //...

    //convert to double
    std::vector<double> to_double() const {
        return {x, y};
    }

    //constructor(s)
    Vec2D() : x(0), y(0) {} // Default constructor
    Vec2D(double x, double y) : x(x), y(y) {} // Constructor with parameters
    Vec2D(const std::vector<double>& vec) { //convert from legacy vector
        //ensure vector size, else default to zero
        if (vec.size() >= 2) {
            x = vec[0];
            y = vec[1];
        } else {
            x = y = 0; // Default to zero if vector size is insufficient
        }
    }
};

struct iVec3D {
    int x, y, z;
    bool empty = false;

    //basic operators

    // addition
    iVec3D operator+(const iVec3D& other) const {
        return iVec3D(x + other.x, y + other.y, z + other.z);
    }

    // subtraction
    iVec3D operator-(const iVec3D& other) const {
        return iVec3D(x - other.x, y - other.y, z - other.z);
    }

    // scalar multiplication (to divide use the inverse)
    iVec3D operator%(int scalar) const {
        return iVec3D(x * scalar, y * scalar, z * scalar);
    }

    // element wise multiply
    iVec3D operator*(const iVec3D& other) const {
        return iVec3D(x * other.x, y * other.y, z * other.z);
    }

    // element wise divide
    iVec3D operator/(const iVec3D& other) const {
        return iVec3D(
            other.x != 0 ? x / other.x : 0,
            other.y != 0 ? y / other.y : 0,
            other.z != 0 ? z / other.z : 0
        );
    }

    // normalize self (++)
    // resultant magnitude will always be 255
    iVec3D& operator++(int) {
        double magnitude = sqrt(x * x + y * y + z * z);
        if (magnitude > 0) {
            double scale = (255.0) / magnitude;
            x = static_cast<int>(x * scale);
            y = static_cast<int>(y * scale);
            z = static_cast<int>(z * scale);
        }
        return *this;
    }

    // assignment operator
    iVec3D& operator=(const iVec3D& other) {
        if (this != &other) { // self-assignment check
            x = other.x;
            y = other.y;
            z = other.z;
        }
        return *this;
    }

    // legacy cast (convert to std::vector<int>)
    std::vector<int> to_vec() const {
        return {x, y, z};
    }

    // convert to double
    Vec3D to_double() const {
        return Vec3D(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z));
    }

    // constructor(s)
    iVec3D() : x(0), y(0), z(0) {} // Default constructor
    iVec3D(int x, int y, int z) : x(x), y(y), z(z) {} // Constructor with parameters
    iVec3D(const std::vector<int>& vec) {
        if (vec.size() >= 3) {
            x = vec[0];
            y = vec[1];
            z = vec[2];
        } else {
            x = y = z = 0; // Default to zero if vector size is insufficient
        }
    }
};

struct iVec2D {
    int x, y;
    bool empty = false;

    //basic operators

    // addition
    iVec2D operator+(const iVec2D& other) const {
        return iVec2D(x + other.x, y + other.y);
    }

    // subtraction
    iVec2D operator-(const iVec2D& other) const {
        return iVec2D(x - other.x, y - other.y);
    }

    // scalar multiplication (to divide use the inverse)
    iVec2D operator%(int scalar) const {
        return iVec2D(x * scalar, y * scalar);
    }

    // element wise multiply
    iVec2D operator*(const iVec2D& other) const {
        return iVec2D(x * other.x, y * other.y);
    }

    // element wise divide
    iVec2D operator/(const iVec2D& other) const {
        return iVec2D(
            other.x != 0 ? x / other.x : 0,
            other.y != 0 ? y / other.y : 0
        );
    }

    // normalize self (++)
    // resultant magnitude will always be 255
    iVec2D& operator++(int) {
        double magnitude = sqrt(x * x + y * y);
        if (magnitude > 0) {
            double scale = (255.0) / magnitude;
            x = static_cast<int>(x * scale);
            y = static_cast<int>(y * scale);
        }
        return *this;
    }

    // assignment operator
    iVec2D& operator=(const iVec2D& other) {
        if (this != &other) { // self-assignment check
            x = other.x;
            y = other.y;
        }
        return *this;
    }

    // legacy cast (convert to std::vector<int>)
    std::vector<int> to_vec() const {
        return {x, y};
    }

    // convert to double
    std::vector<double> to_double() const {
        return {static_cast<double>(x), static_cast<double>(y)};
    }

    // constructor(s)
    iVec2D() : x(0), y(0) {} // Default constructor
    iVec2D(int x, int y) : x(x), y(y) {} // Constructor with parameters
    // no legacy casting for ivecs
};


#endif