/*#################################################################################
 * @file objects.h
 * This header file contains definitions of globally referenced objects and classes
 * used in the raytracer.
 * It includes the Ray struct, which contains data related to ray intersections,
#################################################################################*/


#ifndef OBJECTS_H
#define OBJECTS_H

#include "baseClasses.h"

struct Ray {
    //struct contains complex data related to ray intersections
    Vec3D origin;
    Vec3D euler;
    bool intersect = false; //not implemented?
    int objectID = -1; //object ID of the intersected object (-1 for none)
    double depth = 0.0;
    double frensel = 0.0;
    iVec3D color = {20,25,20}; //20,25,20
    Vec3D intersect_point = {0.0,0.0,0.0};
    Vec3D surface_normal = {0.0,0.0,0.0};
    iVec3D reflection_color;
    Vec3D reflection_vec = {0.0,0.0,0.0};

};

#endif