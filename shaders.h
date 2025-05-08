#ifndef SHADERS_H
#define SHADERS_H

#include "baseClasses.h"
#include "objects.h"
#include "functions.h"
#include "mathHelper.h"

//define no namespace

iVec3D gradient_sky_shader(const Ray& ray, iVec3D horizon_color, iVec3D zenith_color);


#endif