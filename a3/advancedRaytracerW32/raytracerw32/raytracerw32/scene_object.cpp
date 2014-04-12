/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"
#define PI 3.1415926535
#define EPSILON 1e-6 
// #define TEXTURE
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
            Point3D p_n = e + t*d;
            Vector3D n = p_n - c;
            n.normalize();
            ray.intersection.uv[0] = 0.5 + atan2(p_n[2], p_n[0])/(2*PI); 
            ray.intersection.uv[1] = 0.5 - asin(p_n[1])/PI;
            // std::cout<<ray.intersection.mat->normal_ind<<std::endl;
#ifdef TEXTURE
            if (ray.intersection.mat->normal_ind)
            {
                Colour n_bump = ray.intersection.mat->normalmap->getColor(ray.intersection.uv);
                Vector3D n_z(0.0, 0.0, 0.1);
                Vector3D n_x = n_z.cross(n);
                n_x.normalize();
                n_z = n_x.cross(n);
                n_z.normalize();
                Vector3D n_new = (n_bump[0]-0.5)*3*n_x + (n_bump[1]-0.5)*3*n_z + (n_bump[2]-0.5)*n;
                //Vector3D n_new = (n_bump[0]-0.5)*n_x + (n_bump[1]-0.5)*n_z + (n_bump[2]-0.5)*n;
                n_new.normalize();
                n = n_new;
            }
#endif
            ray.intersection.point = modelToWorld * p_n;
            ray.intersection.normal = modelToWorld * n;
            ray.intersection.none = false;
            ray.intersection.t_value = t;
            return true;
        }
        else
            return false;
    }
}
bool EnvirSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
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
        float t = fmax(t1, t2);
        if (t > 0)
        {
            if (!ray.intersection.none)
            {
                return false;
            }
            Point3D p_n = e + t*d;
            Vector3D n = c - p_n;
            n.normalize();
            ray.intersection.uv[0] = 0.5 + atan2(p_n[2], p_n[0]+0.000001)/(2*PI); 
            ray.intersection.uv[1] = 0.5 - asin(p_n[1])/PI;
            ray.intersection.point = modelToWorld * p_n;
            ray.intersection.normal = modelToWorld * n;
            ray.intersection.t_value = 10000;
            return true;
        }
        else
            return false;
    }
}

bool UnitRing::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
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
      Vector3D v_d = worldToModel*ray.dir;
      Point3D p_o = worldToModel*ray.origin;
      float t = - p_o[2]/v_d[2];
      Point3D p_i = p_o + t * v_d;
      double radius = sqrt(pow(p_i[0],2)+pow(p_i[1],2));
      
      if (!ray.intersection.none && ray.intersection.t_value < t)
      {
          return false;
      }
      if (t < 0)
      {
          return false;
      }
      if (radius < 0.7)
      {
          return false;
      }
      if (radius > 1.0)
      {
          return false;
      }
       
      Vector3D v_n(0.0, 0.0, 1.0);
      v_n.normalize();
      ray.intersection.point = modelToWorld * p_i;
      ray.intersection.normal = transNorm(worldToModel, v_n); 
      ray.intersection.normal.normalize(); 
      ray.intersection.uv[0] = (radius - 0.7)/0.3; 
      ray.intersection.uv[1] = 0.5;
      ray.intersection.none = false;
      ray.intersection.t_value = t;
      return true;
    
}


bool GeneralObject::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for GeneralObject, which is
	// bounded by unique cube with vertices (1, 1, 0), 
	// (-1, 1, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

      Vector3D dir = worldToModel*ray.dir;
      Vector3D normal;
      Point3D origin = worldToModel*ray.origin;
      Point3D intersect;
      
      // face x = 1 
      float t = 1e5;
      float t_temp = 1e5;
      float t_min = 0;
      float epsilon = 1e-15;
      bool hit_box = false;
      
      if (!ray.intersection.none)
        t = ray.intersection.t_value;
      
      t_temp =  (1 - origin[0])/(dir[0]+epsilon);
      if ((t_temp > t_min) && (t_temp < t))
      {
          intersect = origin + t_temp * dir;
          if ((intersect[1] >= -1 && intersect[1]<=1) && (intersect[2] >= -1 && intersect[2]<=1))
          {
                t = t_temp;
                hit_box = true;
                normal = Vector3D(1, 0, 0);
          }
      }
      t_temp =  (1 - origin[1])/(dir[1]+epsilon);
      if ((t_temp > t_min) && (t_temp < t))
      {
          intersect = origin + t_temp * dir;
          if ((intersect[0] >= -1 && intersect[0]<=1) && (intersect[2] >= -1 && intersect[2]<=1))
          {
                t = t_temp;
                hit_box = true;
                normal = Vector3D(0, 1, 0);
          }
      }
      t_temp =  (1 - origin[2])/(dir[2]+epsilon);
      if ((t_temp > t_min) && (t_temp < t))
      {
          intersect = origin + t_temp * dir;
          if ((intersect[0] >= -1 && intersect[0]<=1) && (intersect[1] >= -1 && intersect[1]<=1))
          {
                t = t_temp;
                hit_box = true;
                normal = Vector3D(0, 0, 1);
          }
      }
      t_temp =  (-1 - origin[0])/(dir[0]+epsilon);
      if ((t_temp > t_min) && (t_temp < t))
      {
          intersect = origin + t_temp * dir;
          if ((intersect[1] >= -1 && intersect[1]<=1) && (intersect[2] >= -1 && intersect[2]<=1))
          {
                t = t_temp;
                hit_box = true;
                normal = Vector3D(-1, 0, 0);
          }
      }
      t_temp =  (-1 - origin[1])/(dir[1]+epsilon);
      if ((t_temp > t_min) && (t_temp < t))
      {
          intersect = origin + t_temp * dir;
          if ((intersect[0] >= -1 && intersect[0]<=1) && (intersect[2] >= -1 && intersect[2]<=1))
          {
                t = t_temp;
                hit_box = true;
                normal = Vector3D(0, -1, 0);
          }
      }
      t_temp =  (-1 - origin[2])/(dir[2]+epsilon);
      if ((t_temp > t_min) && (t_temp < t))
      {
          intersect = origin + t_temp * dir;
          if ((intersect[0] >= -1 && intersect[0]<=1) && (intersect[1] >= -1 && intersect[1]<=1))
          {
                t = t_temp;
                hit_box = true;
                normal = Vector3D(0, 0, -1);
          }
      }
      // std::cout<<t<<std::endl;
      if (hit_box)
      {
         
          Point3D p_0; 
          Point3D p_1; 
          Point3D p_2;
          Point3D p_x;
          Vector3D n_face;
          Vector3D v_10;
          Vector3D v_21;
          Vector3D v_02;
          double t_triangle = 0;
          double t_min = 10000;
          Vector3D n_min;
          vector<Vector3f> verts = getMesh()->getVerts();
          vector<Face> faces = getMesh()->getFaces();
          vector<Vector3f> normals = getMesh()->getNormals();
          for (unsigned int i = 0; i < faces.size(); i++)
          {
                p_0[0] = verts[faces[i].verts.coords[0]].coords[0];
                p_0[1] = verts[faces[i].verts.coords[0]].coords[1];
                p_0[2] = verts[faces[i].verts.coords[0]].coords[2];
                p_1[0] = verts[faces[i].verts.coords[1]].coords[0];
                p_1[1] = verts[faces[i].verts.coords[1]].coords[1];
                p_1[2] = verts[faces[i].verts.coords[1]].coords[2];
                p_2[0] = verts[faces[i].verts.coords[2]].coords[0];
                p_2[1] = verts[faces[i].verts.coords[2]].coords[1];
                p_2[2] = verts[faces[i].verts.coords[2]].coords[2];
                n_face[0] = normals[faces[i].normals.coords[0]].coords[0]; 
                n_face[1] = normals[faces[i].normals.coords[0]].coords[1]; 
                n_face[2] = normals[faces[i].normals.coords[0]].coords[2];
                // std::cout<<n_face<<std::endl;
                // std::cout<<p_0<<p_1<<p_2<<std::endl;
                n_face.normalize();
                t_triangle = - (n_face.dot(origin - p_0))/(dir.dot(n_face)+EPSILON);
                if ((t_triangle < 0) || (t_triangle > t_min))
                    continue;
                p_x = origin + t_triangle * dir;
                v_10 = p_1 - p_0;
                v_02 = p_0 - p_2;
                v_21 = p_2 - p_1;
                if ((n_face.dot(v_10.cross(p_x - p_0))>-EPSILON) && (n_face.dot(v_21.cross(p_x - p_1))>-EPSILON) && (n_face.dot(v_02.cross(p_x - p_2))>-EPSILON))
                {
                      std::cout<<t_triangle<<std::endl;
                      t_min = t_triangle;
                      n_min = n_face;
                }
          }
		  if (!ray.intersection.none && ray.intersection.t_value < t_min)
		  {
			  return false;
		  }
          if (t_min < 10000)
          {

              p_x = origin + t_min * dir;
              ray.intersection.point = modelToWorld * p_x;
              ray.intersection.normal = transNorm(worldToModel, n_min); 
              ray.intersection.normal.normalize(); 
              ray.intersection.none = false;
              ray.intersection.t_value = t_min;
              return true;
          }
      }
      return false; 
}

bool UnitCube::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
	const Matrix4x4& modelToWorld) {
	// TODO: implement intersection code for GeneralObject, which is
	// bounded by unique cube with vertices (1, 1, 0), 
	// (-1, 1, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	Vector3D dir = worldToModel*ray.dir;
	Vector3D normal;
	Point3D origin = worldToModel*ray.origin;
	Point3D intersect;

	// face x = 1 
	float t = 1e5;
	float t_temp = 1e5;
	float t_min = 0;
	float epsilon = 1e-15;
	bool hit_box = false;

	if (!ray.intersection.none)
		t = ray.intersection.t_value;

	t_temp = (1 - origin[0]) / (dir[0] + epsilon);
	if ((t_temp > t_min) && (t_temp < t))
	{
		intersect = origin + t_temp * dir;
		if ((intersect[1] >= -1 && intersect[1] <= 1) && (intersect[2] >= -1 && intersect[2] <= 1))
		{
			t = t_temp;
			hit_box = true;
			normal = Vector3D(1, 0, 0);
		}
	}
	t_temp = (1 - origin[1]) / (dir[1] + epsilon);
	if ((t_temp > t_min) && (t_temp < t))
	{
		intersect = origin + t_temp * dir;
		if ((intersect[0] >= -1 && intersect[0] <= 1) && (intersect[2] >= -1 && intersect[2] <= 1))
		{
			t = t_temp;
			hit_box = true;
			normal = Vector3D(0, 1, 0);
		}
	}
	t_temp = (1 - origin[2]) / (dir[2] + epsilon);
	if ((t_temp > t_min) && (t_temp < t))
	{
		intersect = origin + t_temp * dir;
		if ((intersect[0] >= -1 && intersect[0] <= 1) && (intersect[1] >= -1 && intersect[1] <= 1))
		{
			t = t_temp;
			hit_box = true;
			normal = Vector3D(0, 0, 1);
		}
	}
	t_temp = (-1 - origin[0]) / (dir[0] + epsilon);
	if ((t_temp > t_min) && (t_temp < t))
	{
		intersect = origin + t_temp * dir;
		if ((intersect[1] >= -1 && intersect[1] <= 1) && (intersect[2] >= -1 && intersect[2] <= 1))
		{
			t = t_temp;
			hit_box = true;
			normal = Vector3D(-1, 0, 0);
		}
	}
	t_temp = (-1 - origin[1]) / (dir[1] + epsilon);
	if ((t_temp > t_min) && (t_temp < t))
	{
		intersect = origin + t_temp * dir;
		if ((intersect[0] >= -1 && intersect[0] <= 1) && (intersect[2] >= -1 && intersect[2] <= 1))
		{
			t = t_temp;
			hit_box = true;
			normal = Vector3D(0, -1, 0);
		}
	}
	t_temp = (-1 - origin[2]) / (dir[2] + epsilon);
	if ((t_temp > t_min) && (t_temp < t))
	{
		intersect = origin + t_temp * dir;
		if ((intersect[0] >= -1 && intersect[0] <= 1) && (intersect[1] >= -1 && intersect[1] <= 1))
		{
			t = t_temp;
			hit_box = true;
			normal = Vector3D(0, 0, -1);
		}
	}
	if (hit_box)
	{
		// std::cout<<t<<std::endl;
		Point3D p_x;
		p_x = origin + t * dir;
		ray.intersection.point = modelToWorld * p_x;
		ray.intersection.normal = transNorm(worldToModel, normal);
		ray.intersection.normal.normalize();
		ray.intersection.none = false;
		ray.intersection.t_value = t;
		return true;
	}

	return false;

}



bool UnitCone::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
	const Matrix4x4& modelToWorld) {
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
	d.normalize();
	Point3D c;
	float A = -pow(d[0], 2) - pow(d[1], 2) + pow(d[2], 2);
	float B = (-2*e[0]*d[0]-2*e[1]*d[1]+2*(e[2]-1)*d[2]);
	float C = -pow(e[0], 2) - pow(e[1], 2) + pow(e[2] - 1, 2);
	if ((B*B - 4*A*C) < 0)
	{
		return false;
	}
	else
	{
		float t1 = (-B - sqrt(B*B - 4*A*C)) / (2*A);
		float t2 = (-B + sqrt(B*B - 4*A*C)) / (2*A);
		float t = fmin(t1, t2);
		if (t > 0)
		{
			
			if (!ray.intersection.none && ray.intersection.t_value < t)
			{
				return false;
			}
			Point3D p_n = e + t*d;
			if ((p_n[2] > 1) || (p_n[2] < 0))
			{
				return false;
			}
			Vector3D n;
			n[0] = p_n[0];
			n[1] = p_n[1];
			n[2] = sqrt(pow(p_n[0],2)+pow(p_n[1],2));
			n.normalize();
			if (d.dot(n)>0)
			{
				n = (-1)*n;
			}
			ray.intersection.point = modelToWorld * p_n;
			ray.intersection.normal = modelToWorld * n;
			ray.intersection.none = false;
			ray.intersection.t_value = t;
			return true;
		}
		else
			return false;
	}
}


bool UnitWindow::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
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
      if ((p_i[0] > 0.1)&&(p_i[0]<0.15) && (p_i[1] > 0.1)&&(p_i[1]<0.15))
      {
          return false;
      }
      if ((p_i[0] > 0.1)&&(p_i[0]<0.15) && (p_i[1] < -0.1)&&(p_i[1]>-0.15))
      {
          return false;
      }
      if ((p_i[0] < -0.1)&&(p_i[0]>-0.15) && (p_i[1] > 0.1)&&(p_i[1]<0.15))
      {
          return false;
      }
      if ((p_i[0] < -0.1)&&(p_i[0]>-0.15) && (p_i[1] < -0.1)&&(p_i[1]>-0.15))
      {
          return false;
      }

      Vector3D v_n(0.0, 0.0, 1.0);
      v_n.normalize();
      ray.intersection.point = modelToWorld * p_i;
      ray.intersection.normal = transNorm(worldToModel, v_n); 
      ray.intersection.normal.normalize(); 
      ray.intersection.none = false;
      ray.intersection.t_value = t;
      return true;
}
