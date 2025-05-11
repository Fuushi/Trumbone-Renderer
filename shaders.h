#ifndef SHADERS_H
#define SHADERS_H

#include "baseClasses.h"
#include "objects.h"
#include "functions.h"
#include "mathHelper.h"

//define no namespace

// sky shaders

//gradient sky shader preset (node)
iVec3D gradient_sky_shader(const Ray& ray, iVec3D horizon_color, iVec3D zenith_color);

// geometry shaders

//frensel shader (should be able to use precomputed values)
iVec3D frensel_shader(const Ray& ray);

//diffuse geometry shader 

//gloss geometry shader

//fast specular shader (ignores gi)

//precise specular shader (gi)

//principal shader

#endif