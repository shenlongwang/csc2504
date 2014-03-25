/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <cstdlib>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
    Vector3D v_l = _pos - ray.intersection.point;
    v_l.normalize();
    // double perturb = 0.25;
    // Vector3D v_rand = Vector3D(perturb*double(rand())/RAND_MAX, perturb*double(rand())/RAND_MAX, perturb*double(rand())/RAND_MAX);
    // Vector3D v_n = ray.intersection.normal+v_rand;
    Vector3D v_n = ray.intersection.normal;
    v_n.normalize();
    Vector3D v_v = -ray.dir;
    v_v.normalize();
    Vector3D v_h = 2*v_l.dot(v_n)*v_n-v_l;
    v_h.normalize();
    double n = ray.intersection.mat->specular_exp;
    Colour col_d = fmax(0, v_n.dot(v_l))*ray.intersection.mat->diffuse*_col_diffuse;
    Colour col_a = ray.intersection.mat->ambient*_col_ambient;
    Colour col_s = pow(fmax(0, v_v.dot(v_h)), n)*ray.intersection.mat->specular*_col_specular;
    ray.col = col_d + col_a + col_s;   
    // ray.col = ray.intersection.mat->diffuse;
    // ray.col = col_d + col_a;
    ray.col.clamp();
}

