/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

//      Vector3D v_d = worldToModel*ray.dir;
//      Point3D p_e = worldToModel*ray.origin;
//      Point3D p_a = Point3D(-0.5, -0.5, 0.0);
//      Point3D p_b = Point3D(0.5, -0.5, 0.0);
//      Point3D p_c = Point3D(-0.5, 0.5, 0.0);
//      Vector3D ab = p_a - p_b;
//      Vector3D ac = p_a - p_c;
//      Vector3D ae = p_a - p_e;
//      float a = ab[0];
//      float b = ab[1];
//      float c = ab[2];
//      float d = ac[0];
//      float e = ac[1];
//      float f = ac[2];
//      float g = v_d[0];
//      float h = v_d[1];
//      float i = v_d[2];
//      float j = ae[0];
//      float k = ae[1];
//      float l = ae[2];
//      float M = a*(e*i-h*f)+b*(g*f-d*i)+c*(d*h-e*g);
//      float t = -(f*(a*k-j*b)+e*(j*c-a*l)+d*(b*l-k*c))/M;
//      if (!ray.intersection.none && ray.intersection.t_value < t)
//      {
//          return false;
//      }
//      if (t < 0)
//      {
//          return false;
//      }
//      float gamma = (i*(a*k-j*b)+h*(j*c-a*l)+g*(b*l-k*c))/M;
//      if ((gamma < 0)||(gamma>1))
//      {
//          return false;
//      }
//      float beta = (j*(e*i-h*f)+k*(g*f-d*i)+l*(d*h-e*g))/M;
//      if ((beta < 0)||(beta>1))
//      {
//          return false;
//      }
//      
//      ray.intersection.point = modelToWorld * (p_e + t*v_d);
//      ray.intersection.normal = transNorm(modelToWorld, Vector3D(0.0, 0.0, 1.0));  
//      ray.intersection.none = false;
//      ray.intersection.t_value = t;
//  	return true;
      Vector3D v_d = worldToModel*ray.dir;
      Point3D p_o = worldToModel*ray.origin;
      float t = - p_o[2]/v_d[2];
      Point3D p_i = p_o + t * v_d;
      
      if (!ray.intersection.none && ray.intersection.t_value < t)
      {
          return false;
      }
      if (t < 0)
      {
          return false;
      }
      if ((p_i[0] > 0.5)||(p_i[1]>0.5))
      {
          return false;
      }
      if ((p_i[0] < -0.5)||(p_i[1]<-0.5))
      {
          return false;
      }
      Vector3D v_n(0.0, 0.0, 1.0);
      v_n.normalize();
      ray.intersection.point = modelToWorld * p_i;
      // ray.intersection.normal = transNorm(modelToWorld, v_n); 
      ray.intersection.normal = transNorm(worldToModel, v_n); 
      ray.intersection.normal.normalize(); 
      ray.intersection.none = false;
      ray.intersection.t_value = t;
      return true;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
    
    Vector3D d = worldToModel*ray.dir;
    Point3D e = worldToModel*ray.origin;
    Point3D c;
    float A = d.dot(d);
    float B = d.dot(e-c);
    float C = (e-c).dot(e-c)-1.0;
    if ((B*B - A*C) < 0)
    {
        return false;
    }
    else
    {
        float t1 = (-B - sqrt(B*B - A*C)) / A;
        float t2 = (-B + sqrt(B*B - A*C)) / A;
        float t = fmin(t1, t2);
        if (t > 0)
        {
            if (!ray.intersection.none && ray.intersection.t_value < t)
            {
                return false;
            }
            Vector3D n = e + t*d - c;
            n.normalize();
            ray.intersection.point = modelToWorld * (e + t*d);
            ray.intersection.normal = modelToWorld * n;
            ray.intersection.none = false;
            ray.intersection.t_value = t;
            return true;
        }
        else
            return false;
    }
}

