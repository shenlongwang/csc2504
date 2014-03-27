/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <omp.h>
#include <ctime>
#include <iostream>
#include <cstdlib>

#define SHADOW true 
#define PHONG false 
#define DEPTH 4


// #define SCENE1
#define SCENE2

#define DEPTH_OF_FIELD
// #define GLOSS
// #define SOFT_SHADOW
// #define ANTI_ALIASING

#define REFLECTION
#define REFRACTION 
Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	_modelToWorld = _modelToWorld*node->trans;
	_worldToModel = node->invtrans*_worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	_worldToModel = node->trans*_worldToModel;
	_modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
    Colour rayAccumulater(0.0, 0.0, 0.0);
    int light_ind = 1;
	for (light_ind = 1;light_ind < 100;light_ind++) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// Implement shadows here if needed.
        Ray3D rayLight;
        Colour raySoftShadowColour;
#ifdef SOFT_SHADOW
        double light_size = 0.5;
        int soft_size = 20;
#else
        double light_size = 0.0;
        int soft_size = 1;
#endif
        double softWeight = 1.0 / double(soft_size);
        Vector3D rayDir1 = Vector3D(0.0, 1.0, 0.0);
        Vector3D rayDir2 = Vector3D(1.0, 0.0, 0.0);
        double r_i = 0.0;
        double r_j = 0.0;
        for (int i = 0; i < soft_size; i++)
        {
            r_i = (((double) rand() / (RAND_MAX))-0.5)*light_size;
            r_j = (((double) rand() / (RAND_MAX))-0.5)*light_size;
            rayLight.dir = curLight->light->get_position() + r_i * rayDir1 + r_j * rayDir2 - ray.intersection.point;
            // rayLight.dir.normalize();
            // std::cout<<rayLight.dir<<std::endl;
            rayLight.origin = ray.intersection.point+0.001*rayLight.dir;
            traverseScene(_root, rayLight); 
            if (!rayLight.intersection.none && rayLight.intersection.t_value > 0.0)
                curLight->light->shade(ray, SHADOW);
        	else
                curLight->light->shade(ray, PHONG);
            raySoftShadowColour = raySoftShadowColour + softWeight * ray.col;
        }
        rayAccumulater = rayAccumulater + raySoftShadowColour;
		curLight = curLight->next;
	}
    rayAccumulater.clamp();
    ray.col = rayAccumulater;
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::initAntiAliasPixelBuffer() {
    int antiWidth = _scrWidth * anti_size;
    int antiHeight = _scrHeight * anti_size;
    int numbytes = antiWidth * antiHeight * sizeof(unsigned char);
	anti_rbuffer = new unsigned char[numbytes];
	anti_gbuffer = new unsigned char[numbytes];
	anti_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < antiHeight; i++) {
		for (int j = 0; j < antiWidth; j++) {
			anti_rbuffer[i*antiWidth+j] = 0;
			anti_gbuffer[i*antiWidth+j] = 0;
			anti_bbuffer[i*antiWidth+j] = 0;
		}
	}
}

void Raytracer::antiDownsampling() {
    int antiWidth = _scrWidth * anti_size;
    int antiHeight = _scrHeight * anti_size;
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
    for (int i = 0; i < _scrHeight; i++) {
    	for (int j = 0; j < _scrWidth; j++) {
            double temp_r = 0.0, temp_g = 0.0, temp_b = 0.0;
            for (int c = 0; c < anti_size; c++)
            {
                for (int r = 0; r < anti_size; r++)
                {
                    temp_r = temp_r + (anti_rbuffer[(i*anti_size+c)*_scrWidth + j*anti_size + r]);
                    temp_g = temp_g + (anti_gbuffer[(i*anti_size+c)*_scrWidth + j*anti_size + r]);
                    temp_b = temp_b + (anti_bbuffer[(i*anti_size+c)*_scrWidth + j*anti_size + r]);
            	}
            }
            temp_r = temp_r / (anti_size*anti_size);
            temp_g = temp_g / (anti_size*anti_size);
            temp_b = temp_b / (anti_size*anti_size);

    		_rbuffer[i*_scrWidth+j] = fmin(255.0, fmax(temp_r, 0.0));
    		_gbuffer[i*_scrWidth+j] = fmin(255.0, fmax(temp_g, 0.0));
    		_bbuffer[i*_scrWidth+j] = fmin(255.0, fmax(temp_b, 0.0));
        
        }
    }
    delete anti_rbuffer;
    delete anti_gbuffer;
    delete anti_bbuffer;
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray, int depth ) {
	Colour col(0.0, 0.0, 0.0); 
	Colour col_reflect(0.0, 0.0, 0.0); 
	Colour col_refract(0.0, 0.0, 0.0); 
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
//*  	if (!ray.intersection.none) {
//*  #ifdef REFLECTION
//*          if (ray.intersection.mat->specular_exp > 0 && ray.intersection.mat->reflect > 0)
//*          {
//*              Ray3D reflect_ray;
//*              double reflect_rate = ray.intersection.mat->reflect;
//*              Vector3D reflect_v = ray.intersection.point - ray.origin;
//*              reflect_v.normalize();
//*              ray.intersection.normal.normalize();
//*              reflect_ray.dir = reflect_v - 2*((reflect_v.dot(ray.intersection.normal))*ray.intersection.normal);
//*              // reflect_ray.dir.normalize();
//*              reflect_ray.origin = ray.intersection.point+0.001*reflect_ray.dir;
//*  	        traverseScene(_root, reflect_ray); 
//*  	        if (!reflect_ray.intersection.none) {
//*                  computeShading(reflect_ray);
//*  		        computeShading(ray); 
//*  		        ray.col = ray.col+reflect_rate*ray.intersection.mat->specular*reflect_ray.col;
//*                  ray.col.clamp();
//*                  col = ray.col;
//*                  return col;
//*              }
//*              else
//*              {
//*                  computeShading(ray); 
//*                  col = ray.col;
//*                  col.clamp();
//*                  return col;
//*              }
//*  	    }
//*          else 
//*          {
//*              computeShading(ray); 
//*              col = ray.col;
//*              col.clamp();
//*              return col;
//*          }
//*  #else
//*  		computeShading(ray); 
//*          col = ray.col;
//*  #endif
//*  	}
//*  
//*  	// You'll want to call shadeRay recursively (with a different ray, 
//*  	// of course) here to implement reflection/refraction effects.  
//*  
//*  	return col; 
   	if ((!ray.intersection.none) && (depth != 0)) {
    #ifdef REFLECTION
            double reflect_rate = ray.intersection.mat->reflect;
            if (ray.intersection.mat->specular_exp > 0 && reflect_rate > 0)
            {
                Ray3D reflect_ray;
                Vector3D reflect_v = ray.intersection.point - ray.origin;
                reflect_v.normalize();
                ray.intersection.normal.normalize();
                reflect_ray.dir = reflect_v - 2*((reflect_v.dot(ray.intersection.normal))*ray.intersection.normal);
                reflect_ray.dir.normalize();
                reflect_ray.origin = ray.intersection.point+0.001*reflect_ray.dir;
                col_reflect = shadeRay(reflect_ray, depth - 1);
            }
    #endif
    #ifdef REFRACTION
            if ( ray.intersection.mat->refract > 0)
            {
                Vector3D refract_d = ray.intersection.point - ray.origin;
                refract_d.normalize();
                Vector3D refract_n = ray.intersection.normal;
                refract_n.normalize();
                double n_div_nt = 1/ray.intersection.mat->refract_ind;
                double d_dot_n = refract_d.dot(refract_n);
                double c = -d_dot_n;
                double cost2 = 1 - pow(n_div_nt,2)*(1-pow(d_dot_n,2));
                if (cost2 > 0.0) {
                    Ray3D refract_ray;
                    Vector3D refract_t = n_div_nt * refract_d + (n_div_nt * c - sqrt(cost2)) * refract_n;
                    refract_ray.dir = refract_t;
                    refract_ray.origin = ray.intersection.point+0.0001*refract_ray.dir;
                    col_refract = shadeRay(refract_ray, depth - 1);
                    col_refract.clamp();
                }
                
// Comment at 08:01 pm Mar 26;                
//                  if (d_dot_n < 0) {
//                      refract_d.normalize();
//                      refract_n.normalize();
//                      Vector3D refract_t1 = n_div_nt*(refract_d - d_dot_n*refract_n);
//                      Vector3D refract_t2 = sqrt(1 - pow(n_div_nt,2)*(1-pow(d_dot_n,2)));
//                      Vector3D refract_t = refract_t1 + refract_t2;
//                      refract_ray.dir = refract_t;
//                      refract_ray.origin = ray.intersection.point+0.0001*refract_ray.dir;
//                  }
//                  else {
//                      double n_div_nt = ray.intersection.mat->refractive_ind;
//                      refract_d.normalize();
//                      refract_n = -refract_n;
//                      refract_n.normalize();
//                      Vector3D refract_t1 = n_div_nt*(refract_d - d_dot_n*refract_n);
//                      Vector3D refract_t2 = sqrt(1 - pow(n_div_nt,2)*(1-pow(d_dot_n,2)));
//                      Vector3D refract_t = refract_t1 + refract_t2;
//                      refract_ray.dir = refract_t;
//                      refract_ray.origin = ray.intersection.point+0.0001*refract_ray.dir;
//                      if (!refract_ray.intersection.none)
//                          c = refract_t.dot(-refract_n);
//                      else {
//                          ray.col = refract_col * col_reflect;
//                          ray.col.clamp();
//                          return ray.col;
//                      }
//                  }
//                  double R_0 = (refract_n - 1)^2 / (refract_n + 1)^2;
//                  double R = R_0 + (1 - R_0) * (1 - c)^5;
//                  computeShading(refract_ray, depth - 1);
//      		    col_refract = R*refract_ray.col;
//                  col_refract.clamp();
//                  }
            }
    
    #endif
    		computeShading(ray);
    #ifdef REFLECTION
            ray.col = ray.col + reflect_rate*ray.intersection.mat->specular*col_reflect;
    #endif
    #ifdef REFRACTION
            ray.col = ray.col + col_refract;
    #endif
            ray.col.clamp();
            col = ray.col;
            return col;
    	}
    
        if ((!ray.intersection.none) && (depth == 0)) {
    		computeShading(ray);
            col = ray.col;
            col.clamp();
            return col;
        }
    	// You'll want to call shadeRay recursively (with a different ray, 
    	// of course) here to implement reflection/refraction effects.  
    
    	return col; 
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);
	
    Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);
    int tempHeight = _scrHeight;
    int tempWidth = _scrWidth;
	double factor = (double(tempHeight)/2)/tan(fov*M_PI/360.0);
	// Construct a ray for each pixel.
#ifdef ANTI_ALIASING
    anti_size = 4;
#else
    anti_size = 1;
#endif
    double r_i = 0.0;
    double r_j = 0.0;
    double anti_step = 1.0/double(anti_size);
    double anti_weight = 1.0/double(anti_size*anti_size);
#ifdef DEPTH_OF_FIELD
    int dof_size = 40;
    double d_j = 0.0;
    double d_i = 0.0;
    double dof_mag = 1*1e-1;
    double dof_weight = 1.0/double(dof_size);
#endif
    for (int i = 0; i < tempHeight; i++) {
		for (int j = 0; j < tempWidth; j++) {
			// Sets up ray origin and direction in view space, 
			// image plane is at z = -1.
			Point3D origin(0, 0, 0);
			Point3D imagePlane;
            Colour col;
            // Antialiasing with jittering sampling. 
            // Reference: Sectoin 13.4.1 Fundamentals of Computer Graphics
            for (int anti_i = 0; anti_i < anti_size; anti_i++)
                for (int anti_j = 0; anti_j < anti_size; anti_j++)
                {
#ifdef ANTI_ALIASING
                    r_i = ((double) rand() / (RAND_MAX));
                    r_j = ((double) rand() / (RAND_MAX));
#endif
            		imagePlane[0] = (-double(tempWidth)/2 + j + (r_j+anti_j)*anti_step)/factor;
            		imagePlane[1] = (-double(tempHeight)/2 + i + (r_i+anti_i)*anti_step)/factor;
            		imagePlane[2] = -1;

            		// TODO: Convert ray to world space and call 
            		// shadeRay(ray) to generate pixel colour. 	
            		
            		Ray3D ray;
#ifdef DEPTH_OF_FIELD
                    Colour col_dof;
                    for (int dof_i = 0; dof_i < dof_size; dof_i++)
                    {     
                        d_i = ((double) rand() / (RAND_MAX) - 0.5)*dof_mag;
                        d_j = ((double) rand() / (RAND_MAX) - 0.5)*dof_mag;
	                    viewToWorld = initInvViewMatrix(eye + d_i*up + d_j*w, view, up);
                        ray.origin = viewToWorld * origin;
                        ray.dir = viewToWorld * imagePlane - ray.origin;
                        ray.dir.normalize();
                        col_dof = col_dof + dof_weight*shadeRay(ray, DEPTH);
                    }     
            		col = col + anti_weight*col_dof; 
#else
                    ray.origin = viewToWorld * origin;
                    ray.dir = viewToWorld * imagePlane - ray.origin;
                    ray.dir.normalize();
            		col = col + anti_weight*shadeRay(ray, 3); 
#endif
                }
            _rbuffer[i*width+j] = int(col[0]*255);
		    _gbuffer[i*width+j] = int(col[1]*255);
		    _bbuffer[i*width+j] = int(col[2]*255);
            std::cout<< "\r Rendering pixel("<<i<<", "<<j<<")...";
		}
        if (i%10 == 0)
            std::cout<<"Rendering "<<i<<"th Row..."<<std::endl;
	}
	flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 320; 
	int height = 240; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}
    time_t timer1, timer2;
    time(&timer1);
	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2, 1.0, 0.0, 0.0 );
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.5, 0.0, 0.0  );
	Material red( Colour(0.0, 0.0, 0.0), Colour(0.9, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.5, 0.0, 0.0  );
	Material yellow( Colour(0.0, 0.0, 0.0), Colour(0.9, 0.9, 0.0), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.5, 0.0, 0.0  );
	Material green( Colour(0.0, 0.0, 0.0), Colour(0.3, 0.9, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.5, 0.0, 0.0  );
	Material cyan( Colour(0.0, 0.0, 0.0), Colour(0.9, 0.0, 0.9), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.5, 0.0, 0.0  );
	Material blue( Colour(0.0, 0.0, 0.0), Colour(0.3, 0.3, 0.9), 
			Colour(0.316228, 0.316228, 0.616228), 
			52.8, 1.0, 0.0, 0.0  );
	Material black( Colour(0.0, 0.0, 0.0), Colour(0.4, 0.4, 0.4), 
			Colour(0.316228, 0.316228, 0.316228), 
			52.8, 1.0, 0.0, 0.0  );
	Material mirror( Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0), 
			Colour(0.8, 0.8, 0.8), 
			12.8, 1.0, 0.0, 0.0  );
	Material glass( Colour(0.0, 0.0, 0.0), Colour(0.4, 0.4, 0.4), 
			Colour(0.8, 0.6, 0.8), 
			52.8, 0.2, 1.0, 1.4 );

#ifdef SCENE1
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 5), 
				Colour(0.6, 0.6, 0.6) ) );
	raytracer.addLightSource( new PointLight(Point3D(-4, -4, 4), 
				Colour(0.4, 0.4, 0.4) ) );
	raytracer.addLightSource( new PointLight(Point3D(4, 4, 4), 
				Colour(0.4, 0.4, 0.4) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* plane_right = raytracer.addObject( new UnitSquare(), &red );
	SceneDagNode* plane_up = raytracer.addObject( new UnitSquare(), &yellow );
	SceneDagNode* plane_left = raytracer.addObject( new UnitSquare(), &mirror);
	SceneDagNode* plane_down = raytracer.addObject( new UnitSquare(), &green );
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* sphere_2 = raytracer.addObject( new UnitSphere(), &blue );
	SceneDagNode* sphere_3 = raytracer.addObject( new UnitSphere(), &black );
	SceneDagNode* sphere_4 = raytracer.addObject( new UnitSphere(), &glass );
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 12.0, 12.0, 12.0 };
	raytracer.translate(sphere, Vector3D(0, 0, -5));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(sphere_2, Vector3D(-4, -0, 2));	
	raytracer.translate(sphere_3, Vector3D(3, 3, 4));	
	raytracer.translate(sphere_4, Vector3D(0, 0, 9));	
	
    raytracer.translate(plane, Vector3D(0, 0, -7));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	
	raytracer.rotate(plane_up, 'z', 45); 
    raytracer.translate(plane_up, Vector3D(0, 6, -2));	
	raytracer.rotate(plane_up, 'x', 90); 
	raytracer.scale(plane_up, Point3D(0, 0, 0), factor2);
	
    raytracer.rotate(plane_left, 'z', 45); 
	raytracer.translate(plane_left, Vector3D(-6, 0, -2));	
	raytracer.rotate(plane_left, 'y', 90); 
	raytracer.scale(plane_left, Point3D(0, 0, 0), factor2);
	
	raytracer.rotate(plane_right, 'z', 45); 
	raytracer.translate(plane_right, Vector3D(6, 0, -2));	
	raytracer.rotate(plane_right, 'y', -90); 
	raytracer.scale(plane_right, Point3D(0, 0, 0), factor2);
	
	raytracer.rotate(plane_down, 'z', 45); 
    raytracer.translate(plane_down, Vector3D(0, -6, -2));	
	raytracer.rotate(plane_down, 'x', -90); 
	raytracer.scale(plane_down, Point3D(0, 0, 0), factor2);
	// Render the scene, feel free to make the image smaller for
	// testing purposes.
	// Camera parameters.
	Vector3D up(0, 1, 0);
	double fov = 60;
	Point3D eye(0, 0, 14);
	Vector3D view(0, 0, -9);
    raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	Point3D eye2(3, 3, 4);
	Vector3D view2(-3, -3, -3);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
#endif	
#ifdef SCENE2
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(-6, 8, 0), 
				Colour(0.4, 0.4, 0.4) ) );
	raytracer.addLightSource( new PointLight(Point3D(4, 8, 2), 
				Colour(0.4, 0.4, 0.4) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &gold );
	// SceneDagNode* plane_2 = raytracer.addObject( new UnitSquare(), &glass );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &glass );
	SceneDagNode* sphere_2 = raytracer.addObject( new UnitSphere(), &blue );
	SceneDagNode* sphere_3 = raytracer.addObject( new UnitSphere(), &mirror );
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 3.0, 3.0, 3.0 };
	double factor2[3] = { 32.0, 32.0, 32.0 };
	raytracer.translate(sphere, Vector3D(-2, 4, -1));	
	raytracer.translate(sphere_2, Vector3D(0, 2, 1));	
	raytracer.translate(sphere_3, Vector3D(1, 6, -0.5));	

	raytracer.rotate(plane, 'x', -90); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	
	// raytracer.translate(plane_2, Vector3D(0, 5, 0));	
	// raytracer.rotate(plane_2, 'x', -90); 
	// raytracer.scale(plane_2, Point3D(0, 0, 0), factor1);
	
    Vector3D up(0, 0, 1);
	double fov = 20;
	Point3D eye(0, 10, 0);
	Vector3D view(0, -1, 0);
    raytracer.render(width, height, eye, view, up, fov, "view1_scene2.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	double fov2 = 60;
    Vector3D up2(0, 1, 3);
	Point3D eye2(0, 8, -2);
	Vector3D view2(0, -3, 1);
	raytracer.render(width, height, eye2, view2, up2, fov2, "view2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
#endif	
    
    return 0;
}

