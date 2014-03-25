/*********************************************************** Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <cstdlib>
#include "light_source.h"

// #define GLOSS

void PointLight::shade( Ray3D& ray, const bool shadow) {
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
    if (shadow)
        ray.col = 0.3*ray.intersection.mat->diffuse*_col_ambient;
    else
    {
        Colour col_d = fmax(0, v_n.dot(v_l))*ray.intersection.mat->diffuse*_col_diffuse;
        Colour col_a = ray.intersection.mat->ambient*_col_ambient;
#ifdef GLOSS
        Colour col_s;
        int gloss_size = 4;
        double gloss_weight = 1.0 / double(gloss_size);
        double gloss_mag = 0.4;
        double r_j = 0.2;
        double r_i = 0.4;
        Vector3D gloss_u = v_h.cross(v_n);
        Vector3D gloss_v = v_h.cross(gloss_u);
        gloss_u.normalize();
        gloss_v.normalize();
        for (int i = 0; i < gloss_size; i++)
        {
            r_i = ((double) rand() / (RAND_MAX) - 0.5)*gloss_mag;
            r_j = ((double) rand() / (RAND_MAX) - 0.5)*gloss_mag;
            v_h = v_h + r_i * gloss_u + r_j * gloss_v;
            v_h.normalize();
            col_s = col_s + gloss_weight * pow(fmax(0, v_v.dot(v_h)), n)*ray.intersection.mat->specular*_col_specular;
        }
#else
        Colour col_s = pow(fmax(0, v_v.dot(v_h)), n)*ray.intersection.mat->specular*_col_specular;
#endif
        
        ray.col = col_d + col_a + col_s; 
    }
    ray.col.clamp();
}

