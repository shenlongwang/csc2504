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
	
	return node;;
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
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// Implement shadows here if needed.
        Ray3D rayLight;
        rayLight.dir = curLight->light->get_position() - ray.intersection.point;
        rayLight.dir.normalize();
        rayLight.origin = ray.intersection.point+0.001*rayLight.dir;
	    traverseScene(_root, rayLight); 
        if (!rayLight.intersection.none && rayLight.intersection.t_value > 0.0)
            curLight->light->shade(ray, SHADOW);
		else 
            curLight->light->shade(ray, PHONG);
		curLight = curLight->next;
	}
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

Colour Raytracer::shadeRay( Ray3D& ray ) {
	Colour col(0.0, 0.0, 0.0); 
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray); 
		col = ray.col; 
        // std::cout<<col<<std::endl;
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
    anti_size = 3;
    double anti_step = 1.0/double(anti_size);
    double anti_weight = 1.0/double((anti_size*2-1)*(anti_size*2-1));
    std::cout<<anti_weight<<std::endl;
    for (int i = 0; i < tempHeight; i++) {
		for (int j = 0; j < tempWidth; j++) {
			// Sets up ray origin and direction in view space, 
			// image plane is at z = -1.
			Point3D origin(0, 0, 0);
			Point3D imagePlane;
            Colour col;
            for (int anti_i = -anti_size+1; anti_i < anti_size; anti_i++)
// #pragma omp parallel for schedule(dynamic, 1)
                for (int anti_j = -anti_size+1; anti_j < anti_size; anti_j++)
                {
            		imagePlane[0] = (-double(tempWidth)/2 + 0.5 + j+anti_j*anti_step)/factor;
            		imagePlane[1] = (-double(tempHeight)/2 + 0.5 + i+anti_i*anti_step)/factor;
            		imagePlane[2] = -1;

            		// TODO: Convert ray to world space and call 
            		// shadeRay(ray) to generate pixel colour. 	
            		
            		Ray3D ray;
                    ray.origin = viewToWorld * origin;
                    ray.dir = viewToWorld * imagePlane - ray.origin;
                    ray.dir.normalize();
            		col = col + anti_weight*shadeRay(ray); 
                }
            _rbuffer[i*width+j] = int(col[0]*255);
		    _gbuffer[i*width+j] = int(col[1]*255);
		    _bbuffer[i*width+j] = int(col[2]*255);
		}
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
	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2 );
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8 );

	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 5), 
				Colour(0.9, 0.9, 0.9) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	raytracer.translate(sphere, Vector3D(0, 0, -5));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane, Vector3D(0, 0, -7));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
	return 0;
}

