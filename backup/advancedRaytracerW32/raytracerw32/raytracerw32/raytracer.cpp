/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "util.h"
#include "bmp_io.h"
#include "meshObject.h"
#include <cmath>
#include <omp.h>
#include <ctime>
#include <iostream>
#include <cstdlib>

//---------------------------------------------------------------------------------
// define parameters for shading, does not mean disable PHONG shading
#define SHADOW true 
#define PHONG false 
#define ENVIRLIGHT 2 
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// define the depth of reflection and refraction. The deeper, the more time comsuming
#define DEPTH 4
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// define the scene to render, please only enable one each time.

// #define SCENE1
// #define SCENE2
// #define SCENE3
// #define SCENE4
// #define SCENE5
// #define SCENE6
// #define SCENE7
// #define SCENE8
#define SCENE9

//---------------------------------------------------------------------------------

// define the resolution.
#define SMALL
// #define SMALL_CUBE
#define LARGE_CUBE

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// define whether environment mapping is used. If no, disable that.
// #define ENVIR
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// define whether depth of field effects is on.
// #define DEPTH_OF_FIELD
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// define whether GLOSS reflection is on.
#define GLOSSY
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// define whether soft shadow is on
#define SOFT_SHADOW
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// define whether anti-aliasing is on 
#define ANTI_ALIASING
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// define whether reflection is on 
#define REFLECTION
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// define whether refraction is on 
 #define REFRACTION 
//---------------------------------------------------------------------------------

using std::vector;

//double* multi_vector(double s, const double* v)
//{
//  double t[3] = {s*v[0], s*v[1], s*v[2]};
//  return t;
//}

void dumpvectorfaces(vector<Face> vec)
{
    for (unsigned int i=0; i < vec.size(); i+=3)
    {
        printf("Vertex Points: (%i,%i,%i)\n",vec[i].verts.coords[0],vec[i].verts.coords[1],vec[i].verts.coords[2]);
        printf("Texture Coords: (%i,%i,%i)\n",vec[i].texCoords.coords[0],vec[i].texCoords.coords[1],vec[i].texCoords.coords[2]);
        printf("Normal Vector: (%i,%i,%i)\n",vec[i].normals.coords[0],vec[i].normals.coords[1],vec[i].normals.coords[2]);
    }
}

void dumpvectorf(vector<Vector3f> vec)
{
    for (unsigned int i=0; i < vec.size(); ++i)
    {
        printf("Vertex Points: (%f,%f,%f)\n",vec[i].coords[0],vec[i].coords[1],vec[i].coords[2]);
    }
}

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

Texture::Texture() {

}

Texture::~Texture() {
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

bool Texture::readImage(char *file_name){
    
	// Pixel buffer.
    long int imHeight = (long int) _texHeight; 
    bool error_flag = bmp_read(file_name, &_texWidth, &imHeight, &_rbuffer , &_gbuffer, &_bbuffer);
    _texHeight = (unsigned long int) imHeight;
    if (error_flag)
        std::cout<<"Texture Image Reading Error!"<<std::endl;
    else
        std::cout<<"Texture Image Reading Success! Image size:"<<_texHeight<<_texWidth<<std::endl;
    
    return error_flag;
}

Colour Texture::getColor(double uv[2]){
   
    double v = 1-fmin(1.0, fmax(uv[0], 0.0));
    double u = 1-fmin(1.0, fmax(uv[1], 0.0));
    int i = (int)(u*(_texHeight-1));
    int j = (int)(v*(_texWidth-1));
    // std::cout<<i<<", "<<j<<std::endl; 
    Colour texColor(0.0, 0.0, 0.0);
    texColor[0] = (1/255.0)*((double)_rbuffer[i*_texWidth+j]);
   	texColor[1] = (1/255.0)*((double)_gbuffer[i*_texWidth+j]);
   	texColor[2] = (1/255.0)*((double)_bbuffer[i*_texWidth+j]);
    // std::cout<< texColor<<std::endl;
    return texColor;
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

void Raytracer::scale2( SceneDagNode* node, Point3D origin, double factor ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor;
	scale[0][3] = origin[0] - factor * origin[0];
	scale[1][1] = factor;
	scale[1][3] = origin[1] - factor * origin[1];
	scale[2][2] = factor;
	scale[2][3] = origin[2] - factor * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor;
	scale[0][3] = origin[0] - 1/factor * origin[0];
	scale[1][1] = 1/factor;
	scale[1][3] = origin[1] - 1/factor * origin[1];
	scale[2][2] = 1/factor;
	scale[2][3] = origin[2] - 1/factor * origin[2];
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


void Raytracer::generatePhotonMap(PhotonMap* photon_map) {
	LightListNode* curLight = _lightSource;
	// photon_map->initialPhotonMap(300);
	// photon_map->displayPhotoMap();
	int light_ind = 1;
	double light_size = 0.4;
	double light_power = 0;
	Point3D light_pos;
	Colour light_col;
	Vector3D photon_dir;
	Vector3D light_u(1.0, 0.0, 0.0);
	Vector3D light_v(0.0, 1.0, 0.0);
	Vector3D light_dir(0.0, 0.0, -1.0);
	double r_i, r_j, r_k;
	int depth = 0;
	int photon_num_unit = 120000;
	int photon_num = 0;
	int photon_counter = 0;
	Ray3D photon_ray;
	for (light_ind = 1; light_ind < 100; light_ind++) 
	{
		
		if (curLight == NULL) 
			break;
		light_pos = curLight->light->get_position();
		light_col = curLight->light->get_col();
		photon_num = (int) photon_num_unit*(sqrt(pow(light_col[0], 2) + pow(light_col[1], 2) + pow(light_col[2], 2))/3);
		while (photon_counter < photon_num)
		{
			photon_ray.intersection.none = true;
			photon_ray.intersection.t_value = 50000;
			r_i = (((double)rand() / (RAND_MAX)) - 0.5) * light_size;
			r_j = (((double)rand() / (RAND_MAX)) - 0.5) * light_size;
			photon_ray.origin = light_pos + r_i * light_u + r_j * light_v;
			do
			{
				r_i = (((double)rand() / (RAND_MAX)) - 0.5)*2.0;
				r_j = (((double)rand() / (RAND_MAX)) - 0.5)*2.0;
				r_k = (((double)rand() / (RAND_MAX)) - 0.5)*2.0;
			} while (pow(r_i, 2) + pow(r_j, 2) + pow(r_k, 2) > 1.0);
			photon_ray.dir = r_k * light_dir + r_i * light_u + r_j * light_v;
			photon_ray.dir.normalize();
			Photon photon_temp(photon_ray.origin, light_col, photon_ray.dir, 0);
			tracePhoton(&photon_temp, 4);
			if (photon_temp.type == ABSORPTION)
			{
				// Specular and diffuse
				if ((photon_temp.specular_history != 0) || (photon_temp.diffuse_history != 0))
				// Only specular
				// if ((photon_temp.specular_history != 0))
				// Only diffuse
				// if ((photon_temp.diffuse_history != 0))
				{
					photon_temp.power.clamp();
					photon_map->addExistPhoton(photon_temp);
					photon_counter++;
					if (!photon_counter % 1000)
						cout << "Already generated " << photon_counter << " photons" << std::endl;
				}
			}
			
			// traverseScene(_root, photon_ray);
			//if (!photon_ray.intersection.none && photon_ray.intersection.t_value > 0.0)
			//{
			//	photon_map->addNewPhoton(photon_ray.intersection.point, photon_ray.intersection.mat->diffuse, photon_ray.dir, 0);
			//	photon_counter++;
			//}
		}
		curLight = curLight->next;
	}
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
        int soft_size = 30;
#else
        double light_size = 0.0;
        int soft_size = 1;
#endif
        double softWeight = 1.0 / double(soft_size);
        Vector3D rayDir1 = Vector3D(0.0, 1.0, 0.0);
        Vector3D rayDir2 = Vector3D(1.0, 0.0, 0.0);
        double r_i = 0.0;
        double r_j = 0.0;
		double t_light = 0.0;
        for (int i = 0; i < soft_size; i++)
        {
            r_i = (((double) rand() / (RAND_MAX))-0.5)*light_size;
            r_j = (((double) rand() / (RAND_MAX))-0.5)*light_size;
            if (!ray.intersection.none && ray.intersection.t_value > 0.0)
            {
                rayLight.dir = curLight->light->get_position() + r_i * rayDir1 + r_j * rayDir2 - ray.intersection.point;
				t_light = rayLight.dir.dot(rayLight.dir);
                rayLight.origin = ray.intersection.point+0.0001*ray.intersection.normal;
                rayLight.dir.normalize(); 
				// to avoid intersection after the light
				t_light = sqrt(t_light / (rayLight.dir.dot(rayLight.dir)));
                traverseScene(_root, rayLight);
				if (!rayLight.intersection.none && rayLight.intersection.t_value > 0.0 && rayLight.intersection.t_value < t_light)
                {
                    curLight->light->shade(ray, SHADOW);
        	    }
                else
                    curLight->light->shade(ray, PHONG);
            }
#ifdef ENVIR
        	else
                curLight->light->shade(ray, SHADOW);
#endif    
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

    		_rbuffer[i*_scrWidth + j] = (unsigned char) fmin(255.0, fmax(temp_r, 0.0));
			_gbuffer[i*_scrWidth + j] = (unsigned char) fmin(255.0, fmax(temp_g, 0.0));
			_bbuffer[i*_scrWidth + j] = (unsigned char) fmin(255.0, fmax(temp_b, 0.0));
        
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
    
   	if ((!ray.intersection.none) && (depth != 0)) {
    // if (ray.intersection.mat->texture_ind)
        // std::cout<<ray.intersection.mat->texture->_texWidth<<std::endl;
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
		#ifdef GLOSSY
				if (ray.intersection.mat->glossy_ind)
				{
					Vector3D u;
					Vector3D v;
					Ray3D reflect_temp;
					Colour colour_temp(0.0, 0.0, 0.0);
					reflect_temp.origin = reflect_ray.origin;
					int glossy_num = 80;
					double glossy_rate = 0.5;
					double r_i = 0;
					double r_j = 0;
					double cos_theta = 0;
					double glossy_weight = 0;
					for (int i = 0; i < glossy_num; i++)
					{
						r_i = (((double)rand() / (RAND_MAX)) - 0.5)*glossy_rate;
						r_j = (((double)rand() / (RAND_MAX)) - 0.5)*glossy_rate;
						u = reflect_ray.dir.cross(ray.intersection.normal);
						v = reflect_ray.dir.cross(u);
						u.normalize();
						v.normalize();
						reflect_temp.dir = reflect_ray.dir + r_i * u + r_j * v;
						reflect_temp.dir.normalize();
						if (reflect_temp.dir.dot(ray.intersection.normal) < 0.0)
						{
							reflect_temp.dir = reflect_ray.dir - r_i * u - r_j * v;
							reflect_temp.dir.normalize();
						}
						// for glossy reflections, reduce the depth to reduce computational cost
						colour_temp = shadeRay(reflect_temp, fmax(depth - 2, 0));
						colour_temp.clamp();
						cos_theta = pow(reflect_temp.dir.dot(reflect_ray.dir),2);
						col_reflect = col_reflect + cos_theta * colour_temp;
						glossy_weight = glossy_weight + cos_theta;
					}
					if (glossy_weight > 0)
					{
						glossy_weight = 1 / (glossy_weight + DBL_EPSILON);
						col_reflect = glossy_weight * col_reflect;
					}
					else
						col_reflect = shadeRay(reflect_ray, depth - 1);
					col_refract.clamp();

				}
				else
				{
					col_reflect = shadeRay(reflect_ray, depth - 1);
					col_reflect.clamp();
				}

		#else
				col_reflect = shadeRay(reflect_ray, depth - 1);
				col_reflect.clamp();
		#endif

            }
    
	#endif             

    #ifdef REFRACTION
			if (ray.intersection.mat->refract > 0)
			{

				Vector3D refract_d = ray.intersection.point - ray.origin;
				refract_d.normalize();
				Vector3D refract_n = ray.intersection.normal;
				refract_n.normalize();
				double n_div_nt = 1 / ray.intersection.mat->refract_ind;
				double d_dot_n = refract_d.dot(refract_n);
				if (d_dot_n > 0) 
				{
					n_div_nt = ray.intersection.mat->refract_ind;
				}
				double c = -d_dot_n;
				double cost2 = 1 - pow(n_div_nt, 2)*(1 - pow(d_dot_n, 2));
				if (cost2 > 0.0) {
					Ray3D refract_ray;
					Vector3D refract_t = n_div_nt * refract_d + (n_div_nt * c - sqrt(cost2)) * refract_n;
					refract_ray.dir = refract_t;
					refract_ray.origin = ray.intersection.point + 0.0001*refract_ray.dir;
					col_refract = shadeRay(refract_ray, depth - 1);
					col_refract.clamp();
				}

            }
    
    #endif
    		computeShading(ray);
    
    #ifdef REFLECTION
            ray.col = ray.col + reflect_rate*ray.intersection.mat->specular*col_reflect;
    #endif
    
    #ifdef REFRACTION
            ray.col = ray.col + ray.intersection.mat->refract*col_refract;
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

#ifdef ENVIR        
        if (ray.intersection.none) {
            computeShading(ray);
            col = ray.col;
            col.clamp();
            return col;
        }
#endif
    	// You'll want to call shadeRay recursively (with a different ray, 
    	// of course) here to implement reflection/refraction effects.  
    
    	return col; 
}	

void Raytracer::tracePhoton(Photon* photon, int depth)
{
	Ray3D photon_ray;
	photon_ray.intersection.none = true;
	photon_ray.intersection.t_value = 40000;
	photon_ray.origin = photon->pos;
	photon_ray.dir = photon->dir;
	traverseScene(_root, photon_ray);
	double p_diffuse, p_specular, p_absorption, p_sum, p_rand;
	double numerator, denominator;
	Colour col_diffuse, col_specular;
	Vector3D reflect_v, refract_v;
	if (!photon_ray.intersection.none && photon_ray.intersection.t_value > 0.0)
	{
		double toss = (double)rand() / (RAND_MAX) * 15.0;
		denominator = fmax(fmax(photon->power[0], photon->power[1]), photon->power[2]) + DBL_EPSILON;
		// the probability of diffuse reflection
		col_diffuse = photon_ray.intersection.mat->diffuse * photon->power;
		numerator = fmax(fmax(col_diffuse[0], col_diffuse[1]), col_diffuse[2]);
		p_diffuse = numerator / denominator;
		// the probability of specular reflection
		col_specular = (photon_ray.intersection.mat->specular * photon->power);
		numerator = fmax(fmax(col_specular[0], col_specular[1]), col_specular[2]);
		p_specular = numerator / denominator;
		if ((p_specular < 0.2) && (photon_ray.intersection.t_value > toss))
		{
			photon->type = OFF;
			return;
		}
		// Russian roulette
		p_sum = p_specular + p_diffuse;
		if (p_sum > 1)
		{
			p_specular = p_specular / p_sum;
			p_diffuse = p_diffuse / p_sum;
			p_sum = 1;
		}
		p_absorption = 1 - p_sum;
		p_rand = (double)rand() / (RAND_MAX);
		if (p_rand <= p_diffuse*0.5)
		{
			photon_ray.intersection.normal.normalize();
			double r_i = (double)rand() / (RAND_MAX);
			double r_j = (double)rand() / (RAND_MAX) - 0.5;
			double r_k = (double)rand() / (RAND_MAX) - 0.5;
			reflect_v = photon_ray.intersection.point - photon_ray.origin;
			reflect_v.normalize();
			reflect_v = reflect_v - 2 * ((reflect_v.dot(photon_ray.intersection.normal))*photon_ray.intersection.normal);
			reflect_v.normalize();
			Vector3D u = reflect_v.cross(photon_ray.dir);
			Vector3D v = reflect_v.cross(v);
			u.normalize();
			v.normalize();
			reflect_v = r_i * reflect_v + r_j * u + r_k * v;
			if (reflect_v.dot(photon_ray.intersection.normal))
				reflect_v = r_i * reflect_v - r_j * u - r_k * v;
			photon->dir = reflect_v;
			photon->dir.normalize();
			photon->pos = photon_ray.intersection.point + 0.00001*photon->dir;
			//photon->power = (1 / (p_diffuse + DBL_EPSILON) / denominator )*col_diffuse;

			if (depth <= 0)
			{
				photon->type = ABSORPTION;
				photon->power.clamp();
				// photon->power = photon_ray.intersection.mat->diffuse;
				// photon->type = ABSORPTION;
			}
			else
			{
				photon->power = (1 / denominator)*col_diffuse;
				photon->power.clamp();
				photon->diffuse_history++;
				photon->type = DIFFUSE;
				tracePhoton(photon, depth - 1);
			}
		}
		else if (p_rand <= p_diffuse*0.5 + p_specular)
		{
			if (depth <= 0)
			{
				photon->type = OFF;
			}
			else
			{
				if (photon_ray.intersection.mat->refract <= photon_ray.intersection.mat->reflect)
				{	// reflection
					reflect_v = photon_ray.intersection.point - photon_ray.origin;
					reflect_v.normalize();
					photon_ray.intersection.normal.normalize();
					photon->dir = reflect_v - 2 * ((reflect_v.dot(photon_ray.intersection.normal))*photon_ray.intersection.normal);
					photon->dir.normalize();
					photon->pos = photon_ray.intersection.point + 0.00001*photon->dir;
					if (photon_ray.intersection.mat->reflect > 0)
						photon->power = (1 / (p_specular + DBL_EPSILON) / denominator)*col_specular;
					photon->power.clamp();
					photon->specular_history++;
					photon->type = SPECULAR;
				}
				else
				{ // refraction
					refract_v = photon_ray.intersection.point - photon_ray.origin;
					refract_v.normalize();
					photon_ray.intersection.normal.normalize();
					double n_div_nt = 1 / photon_ray.intersection.mat->refract_ind;
					double d_dot_n = refract_v.dot(photon_ray.intersection.normal);
					if (d_dot_n > 0)
					{
						n_div_nt = photon_ray.intersection.mat->refract_ind;
					}
					
					double c = -d_dot_n;
					double cost2 = 1 - pow(n_div_nt, 2)*(1 - pow(d_dot_n, 2));
					if (cost2 > 0.0)
						refract_v = n_div_nt * refract_v + (n_div_nt * c - sqrt(cost2)) * photon_ray.intersection.normal;
					else
					{
						photon->type = OFF;
						return;
					}

					photon->dir = refract_v;
					photon->dir.normalize();
					photon->power = (1 / (p_specular + DBL_EPSILON) / denominator)*col_specular;
					photon->pos = photon_ray.intersection.point + 0.00001*photon->dir;
					photon->power.clamp();
					photon->specular_history++;
					photon->type = SPECULAR;
				}
				tracePhoton(photon, depth - 1);
			}


		}
		else if(depth < 4) // To make sure direct photon cannot be absorbed.
		{
			photon->pos = photon_ray.intersection.point;
			photon->dir = photon_ray.intersection.normal;
			photon->dir.normalize();
			//photon->power = col_diffuse;
			photon->power.clamp();
			photon->type = ABSORPTION;
		}

	}
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
	anti_size = 3;
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
            std::cout<< "Rendering pixel("<<i<<", "<<j<<")..."<<std::endl;
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
            		col = col + anti_weight*shadeRay(ray, DEPTH); 
#endif
                }


            _rbuffer[i*width+j] = int(col[0]*255);
		    _gbuffer[i*width+j] = int(col[1]*255);
		    _bbuffer[i*width+j] = int(col[2]*255);
		}
        if (i%1 == 0)
            std::cout<<"Rendering "<<i<<"th Row..."<<std::endl;
	}
	flushPixelBuffer(fileName);
}

void Raytracer::renderPhotonMap(int width, int height, Point3D eye, Vector3D view,
	Vector3D up, double fov, char* fileName, PhotonMap* photon_map) {
	Matrix4x4 viewToWorld;
	Matrix4x4 worldToView;
	_scrWidth = width;
	_scrHeight = height;
	initPixelBuffer();

	viewToWorld = initInvViewMatrix(eye, view, up);
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	double factor = (double(height) / 2) / tan(fov*M_PI / 360.0);

	// Construct a ray for each pixel.
	Point3D origin(0, 0, 0);
	Photon photon_temp;
	Vector3D photon_eye(0, 0, 0);
	Ray3D eye_ray;
	Point3D normalized_photon(0, 0, -1);
	origin = viewToWorld * origin;
	int pixel_i, pixel_j;
	for (int i = 0; i < photon_map->getSize(); i++)
	{
		photon_temp = photon_map->getPhoton(i);
		photon_eye = photon_temp.pos - origin;
		eye_ray.dir = -photon_eye;
		eye_ray.origin = photon_temp.pos + 0.00001 * eye_ray.dir;
		eye_ray.intersection.none = true;
		traverseScene(_root, eye_ray);
		if (eye_ray.intersection.none)
		{
			photon_eye[1] = (photon_eye[1] / (-photon_eye[0] + DBL_EPSILON));
			photon_eye[2] = (photon_eye[2] / (-photon_eye[0] + DBL_EPSILON));
			photon_eye[0] = (photon_eye[0] / (-photon_eye[0] + DBL_EPSILON));
			normalized_photon = origin + photon_eye;
			pixel_i = (int) (normalized_photon[2] * factor + double(width) / 2);
			pixel_j = (int) (normalized_photon[1] * factor + double(height) / 2);
			photon_temp.power.clamp();
			if ((pixel_i < height) && (pixel_j < width) && (pixel_i > 0) && (pixel_j > 0))
			{
				_rbuffer[pixel_i*width + pixel_j] = int(photon_temp.power[0] * 255);
				_gbuffer[pixel_i*width + pixel_j] = int(photon_temp.power[1] * 255);
				_bbuffer[pixel_i*width + pixel_j] = int(photon_temp.power[2] * 255);
			}
		}
	}

	flushPixelBuffer(fileName);

}

void Raytracer::renderWithPhoton(int width, int height, Point3D eye, Vector3D view,
	Vector3D up, double fov, char* fileName, PhotonMap* photon_map) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);
	
	int knn_num = 100;
	double knn_radius = 1;
	Colour photon_col;

	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);
	int tempHeight = _scrHeight;
	int tempWidth = _scrWidth;
	double factor = (double(tempHeight) / 2) / tan(fov*M_PI / 360.0);
	// Construct a ray for each pixel.
#ifdef ANTI_ALIASING
	anti_size = 3;
#else
	anti_size = 1;
#endif
	double r_i = 0.0;
	double r_j = 0.0;
	double anti_step = 1.0 / double(anti_size);
	double anti_weight = 1.0 / double(anti_size*anti_size);
	vector<Photon> photon_knn;
	for (int i = 0; i < tempHeight; i++) {
		for (int j = 0; j < tempWidth; j++) {
			// Sets up ray origin and direction in view space, 
			// image plane is at z = -1.
			Point3D origin(0, 0, 0);
			Point3D imagePlane;
			Colour col(0.0, 0.0, 0.0);
			std::cout << "Rendering pixel(" << i << ", " << j << ")..." << std::endl;
			// Antialiasing with jittering sampling. 
			// Reference: Sectoin 13.4.1 Fundamentals of Computer Graphics


			for (int anti_i = 0; anti_i < anti_size; anti_i++)
			for (int anti_j = 0; anti_j < anti_size; anti_j++)
			{
#ifdef ANTI_ALIASING
				r_i = ((double)rand() / (RAND_MAX));
				r_j = ((double)rand() / (RAND_MAX));
#endif
				imagePlane[0] = (-double(tempWidth) / 2 + j + (r_j + anti_j)*anti_step) / factor;
				imagePlane[1] = (-double(tempHeight) / 2 + i + (r_i + anti_i)*anti_step) / factor;
				imagePlane[2] = -1;

				// TODO: Convert ray to world space and call 
				// shadeRay(ray) to generate pixel colour. 	

				Ray3D ray;
				ray.origin = viewToWorld * origin;
				ray.dir = viewToWorld * imagePlane - ray.origin;
				ray.dir.normalize();
				col = col + anti_weight*shadeRay(ray, DEPTH);
				traverseScene(_root, ray);
				if (!ray.intersection.none && (ray.intersection.mat->reflect == 0) && (ray.intersection.mat->refract == 0))
				{
					photon_col = (1 / 1200.0) * photon_map->findKNN(ray.intersection.point, knn_num);
					// std::cout << photon_col << std::endl;
					photon_col.clamp();
					col = col + anti_weight*photon_col;
					col.clamp();
				}
			}


			_rbuffer[i*width + j] = int(col[0] * 255);
			_gbuffer[i*width + j] = int(col[1] * 255);
			_bbuffer[i*width + j] = int(col[2] * 255);
		}
		if (i % 1 == 0)
			std::cout << "Rendering " << i << "th Row..." << std::endl;
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
#ifdef SMALL
	int width = 480; 
	int height = 320; 
#else
	int width = 1024; 
	int height = 768; 
	// int width = 1920; 
	// int height = 1080; 
#endif
#ifdef SMALL_CUBE
	width = 160;
	height = 160;
#endif
#ifdef LARGE_CUBE
	width = 640;
	height = 640;
#endif


	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}
    time_t timer1, timer2;
    time(&timer1);
	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			12.8, 0.0, 0.0, 0.0 );
	Material gold_mirror(Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
		Colour(0.628281, 0.555802, 0.366065),
		12.8, 0.6, 0.0, 0.0);
	Material jade( Colour(0, 0, 0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.0, 0.0, 0.0  );
	Material jade_mirror(Colour(0, 0, 0), Colour(0.54, 0.89, 0.63),
		Colour(0.316228, 0.316228, 0.316228),
		12.8, 0.5, 0.0, 0.0);
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
			32.8, 0.3, 0.0, 0.0  );
	Material black( Colour(0.0, 0.0, 0.0), Colour(0.4, 0.4, 0.4), 
			Colour(0.316228, 0.316228, 0.316228), 
			52.8, 1.0, 0.0, 0.0  );
	Material mirror( Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0), 
			Colour(0.8, 0.8, 0.8), 
			32.8, 1.0, 0.0, 0.0  );
	Material glass( Colour(0.0, 0.0, 0.0), Colour(0.4, 0.4, 0.4), 
			Colour(0.8, 0.6, 0.8), 
			52.8, 0.1, 1.0, 1.4 );
	Material blue_glass( Colour(0.0, 0.0, 0.0), Colour(0.2, 0.2, 0.4), 
			Colour(0.2, 0.2, 0.2), 
			52.8, 0, 0.5, 1.0 );
	Material red_glass( Colour(0.0, 0.0, 0.0), Colour(0.5, 0.3, 0.0), 
			Colour(0.2, 0.2, 0.2), 
			52.8, 0, 0.5, 1.0 );

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
    raytracer.render(width, height, eye, view, up, fov, "scene1_view1.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	Point3D eye2(3, 3, 4);
	Vector3D view2(-3, -3, -3);
	raytracer.render(width, height, eye2, view2, up, fov, "scene1_view2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
#endif	
#ifdef SCENE2
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(-6, 20, 0), 
				Colour(0.4, 0.4, 0.4) ) );
	raytracer.addLightSource( new PointLight(Point3D(4, 28, 2), 
				Colour(0.4, 0.4, 0.4) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &gold );
	// SceneDagNode* plane_2 = raytracer.addObject( new UnitSquare(), &glass );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &glass );
	SceneDagNode* sphere_2 = raytracer.addObject( new UnitSphere(), &mirror );
	SceneDagNode* sphere_3 = raytracer.addObject( new UnitSphere(), &blue );
    
    plane->mat->texture_ind = false; 
    sphere->mat->texture_ind = false; 
    sphere_2->mat->texture_ind = false; 
    sphere_3->mat->texture_ind = false; 
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 3.0, 3.0, 3.0 };
	double factor2[3] = { 32.0, 32.0, 32.0 };
	raytracer.translate(sphere, Vector3D(-1, 10, -1));	
	raytracer.translate(sphere_2, Vector3D(0, 5, 0));	
	raytracer.translate(sphere_3, Vector3D(1, 15, -0.5));	

	raytracer.rotate(plane, 'x', -90); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	
	// raytracer.translate(plane_2, Vector3D(0, 5, 0));	
	// raytracer.rotate(plane_2, 'x', -90); 
	// raytracer.scale(plane_2, Point3D(0, 0, 0), factor1);
	
    Vector3D up(0, 0, 1);
	double fov = 60;
	Point3D eye(0, 17, 0);
	Vector3D view(0, -1, 0);
    raytracer.render(width, height, eye, view, up, fov, "scene2_view1.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	double fov2 = 60;
    Vector3D up2(0, 1, 3);
	Point3D eye2(0, 18, -6);
	Vector3D view2(0, -3, 1);
	raytracer.render(width, height, eye2, view2, up2, fov2, "scene2_view2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
#endif	
#ifdef SCENE3
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 0), 
				Colour(0.8, 0.8, 0.8) ) );
	raytracer.addLightSource( new PointLight(Point3D(0, 30, 30), 
				Colour(0.4, 0.4, 0.4) ) );
	raytracer.addLightSource( new PointLight(Point3D(0, 30, -30), 
				Colour(0.4, 0.4, 0.4) ) );
	raytracer.addLightSource( new PointLight(Point3D(0, 80, 0), 
				Colour(0.8, 0.8, 0.8) ) );

	Material mat_mercury( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_venus( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_earth( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.9), 
			Colour(0.316228, 0.316228, 0.616228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_mars( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_jupiter( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_saturn( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_uranus( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_neptune( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_galaxy( Colour(0.3, 0.3, 0.3), Colour(0.1, 0.1, 0.1), 
			Colour(0.016228, 0.016228, 0.016228), 
			1.0, 0.0, 0.0, 0.0  );
	Material mat_saturnring( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	Material mat_uranusring( Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0.2, 0.0, 0.0  );
	// Add a unit square into the scene with material mat.
	SceneDagNode* mercury = raytracer.addObject( new UnitSphere(), &mat_mercury );
	SceneDagNode* venus = raytracer.addObject( new UnitSphere(), &mat_venus );
	SceneDagNode* earth = raytracer.addObject( new UnitSphere(), &mat_earth );
	SceneDagNode* mars = raytracer.addObject( new UnitSphere(), &mat_mars );
	SceneDagNode* jupiter = raytracer.addObject( new UnitSphere(), &mat_jupiter);
	SceneDagNode* saturn = raytracer.addObject( new UnitSphere(), &mat_saturn );
	SceneDagNode* uranus = raytracer.addObject( new UnitSphere(), &mat_uranus );
	SceneDagNode* neptune = raytracer.addObject( new UnitSphere(), &mat_neptune );
	SceneDagNode* saturn_ring = raytracer.addObject( new UnitRing(), &mat_saturnring );
	SceneDagNode* uranus_ring = raytracer.addObject( new UnitRing(), &mat_uranusring );
	SceneDagNode* galaxy = raytracer.addObject( new EnvirSphere(), &mat_galaxy );
	
	// Apply some transformations to the unit square.
    // double radius[8] = {5, 12, 13, 7, 142, 120, 51, 49};
    double radius[8] = {40, 45, 53, 40, 162, 120, 81, 69};
    double distance[8] = {20, 200, 440, 700, 1300, 2000, 2500, 3000};
    // double rad_factor[3] = {0.01, 0.01, 0.01};
    Vector3D dis_factor(0.000, 0.02, 0.0);
    double rad_factor = 0.03;
    // 	double radius[8] = {1, 1, 1, 1, 1, 1, 1, 1};
	// double distance[8] = {1, 3, 5, 7, 9, 11, 13};
    //	Vector3D dis_factor(0.0, 1.0, 0.0);
	
    // Add a texture; 
    Texture mercurytex;
    mercurytex.readImage("./Textures/mercurymap.bmp");
    mercury->mat->texture = &mercurytex;  
    mercury->mat->texture_ind = true; 
    // Add a texture; 
    Texture venustex;
    venustex.readImage("./Textures/venusmap.bmp");
    venus->mat->texture = &venustex;  
    venus->mat->texture_ind = true; 
    // Add a texture; 
    Texture earthtex;
    earthtex.readImage("./Textures/earthmap1k.bmp");
    earth->mat->texture = &earthtex;  
    earth->mat->texture_ind = true; 
    // Add a texture; 
    Texture marstex;
    marstex.readImage("./Textures/mars_1k_color.bmp");
    mars->mat->texture = &marstex;  
    mars->mat->texture_ind = true; 
    // Add a texture; 
    Texture jupitertex;
    jupitertex.readImage("./Textures/jupitermap.bmp");
    jupiter->mat->texture = &jupitertex;  
    jupiter->mat->texture_ind = true; 
    // Add a texture; 
    Texture saturntex;
    saturntex.readImage("./Textures/saturnmap.bmp");
    saturn->mat->texture = &saturntex;  
    saturn->mat->texture_ind = true; 
    // Add a texture; 
    Texture uranustex;
    uranustex.readImage("./Textures/uranusmap.bmp");
    uranus->mat->texture = &uranustex;  
    uranus->mat->texture_ind = true; 
    // Add a texture; 
    Texture neptunetex;
    neptunetex.readImage("./Textures/neptunemap.bmp");
    neptune->mat->texture = &neptunetex;  
    neptune->mat->texture_ind = true; 
    // Add a texture; 
    Texture galaxytex;
    galaxytex.readImage("./Textures/galaxy.bmp");
    galaxy->mat->texture = &galaxytex;  
    galaxy->mat->texture_ind = true; 
    
    // Add a texture; 
    Texture saturnringtex;
    saturnringtex.readImage("./Textures/saturnringcolor.bmp");
    saturn_ring->mat->texture = &saturnringtex;  
    saturn_ring->mat->texture_ind = true; 
    
    // Add a texture; 
    Texture uranusringtex;
    uranusringtex.readImage("./Textures/uranusringcolour.bmp");
    uranus_ring->mat->texture = &uranusringtex;  
    uranus_ring->mat->texture_ind = true; 
    
   //   // Add a texture; 
   //   Texture saturnringtrans;
   //   saturnringtrans.readImage("./Textures/saturnringtrans.bmp");
   //   saturn_ring->mat->transparent = &saturnringtrans;  
   //   saturn_ring->mat->transparent_ind = true; 
   //   
   //   // Add a texture; 
   //   Texture uranusringtrans;
   //   uranusringtrans.readImage("./Textures/uranusringtrans.bmp");
   //   uranus_ring->mat->transparent = &uranusringtrans;  
   //   uranus_ring->mat->transparent_ind = true; 
    
    raytracer.translate(mercury, (distance[0]*dis_factor));	
    raytracer.translate(venus, (distance[1]*dis_factor));	
    raytracer.translate(earth, (distance[2]*dis_factor));	
    raytracer.translate(mars, (distance[3]*dis_factor));	
    raytracer.translate(jupiter, (distance[4]*dis_factor));	
    raytracer.translate(saturn, (distance[5]*dis_factor));	
    raytracer.translate(saturn_ring, (distance[5]*dis_factor));	
    raytracer.translate(uranus, (distance[6]*dis_factor));	
    raytracer.translate(uranus_ring, (distance[6]*dis_factor));	
    raytracer.translate(neptune, (distance[7]*dis_factor));	
    raytracer.translate(galaxy, (distance[4]*dis_factor));	
   
   raytracer.scale2(mercury, Point3D(0, 0, 0), rad_factor*radius[0]);
   raytracer.scale2(venus, Point3D(0, 0, 0), rad_factor*radius[1]);
   raytracer.scale2(earth, Point3D(0, 0, 0), rad_factor*radius[2]);
   raytracer.scale2(mars, Point3D(0, 0, 0), rad_factor*radius[3]);
   raytracer.scale2(jupiter, Point3D(0, 0, 0), rad_factor*radius[4]);
   raytracer.scale2(saturn, Point3D(0, 0, 0), rad_factor*radius[5]);
   raytracer.scale2(saturn_ring, Point3D(0, 0, 0), rad_factor*radius[5]*1.6);
   raytracer.scale2(uranus, Point3D(0, 0, 0), rad_factor*radius[6]);
   raytracer.scale2(uranus_ring, Point3D(0, 0, 0), rad_factor*radius[6]*1.6);
   raytracer.scale2(neptune, Point3D(0, 0, 0), rad_factor*radius[7]);
   raytracer.scale2(galaxy, Point3D(0, 0, 0), 50);
    
    // raytracer.scale(mercury, Point3D(0, 0, 0), multi_vector(radius[0],rad_factor));
    // raytracer.scale(venus, Point3D(0, 0, 0), multi_vector(radius[1],rad_factor));
    // raytracer.scale(earth, Point3D(0, 0, 0), multi_vector(radius[2],rad_factor));
    // raytracer.scale(mars, Point3D(0, 0, 0), multi_vector(radius[3],rad_factor));
    // raytracer.scale(jupiter, Point3D(0, 0, 0), multi_vector(radius[4],rad_factor));
    // raytracer.scale(saturn, Point3D(0, 0, 0), multi_vector(radius[5],rad_factor));
    // raytracer.scale(uranus, Point3D(0, 0, 0), multi_vector(radius[6],rad_factor));
    // raytracer.scale(neptune, Point3D(0, 0, 0), multi_vector(radius[7],rad_factor));
	
	// raytracer.translate(plane_2, Vector3D(0, 5, 0));	
	
    raytracer.rotate(earth, 'x', 75); 
    raytracer.rotate(saturn_ring, 'x', 75); 
	raytracer.rotate(saturn_ring, 'z', 30); 
	raytracer.rotate(uranus_ring, 'x', 90); 
	raytracer.rotate(uranus_ring, 'z', 15); 
	// raytracer.scale(plane_2, Point3D(0, 0, 0), factor1);
	
    Vector3D up(-1, 0, 0);
	double fov = 60;
	Point3D eye(0, 28, 44);
	Vector3D view(0, 0, -1);
    
    raytracer.render(width, height, eye, view, up, fov, "solar_1.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	double fov2 = 60;
    Vector3D up2(-1, 0, 0);
	Point3D eye2(0, 7, 12);
	Vector3D view2(0, 0, -1);
	raytracer.render(width, height, eye2, view2, up2, fov2, "solar_2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
#endif	
#ifdef SCENE4
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 0), 
				Colour(0.8, 0.8, 0.8) ) );
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 5), 
				Colour(0.4, 0.4, 0.4) ) );
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 0, -5), 
				Colour(0.4, 0.4, 0.4) ) );
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(5, 0, 0), 
				Colour(0.4, 0.4, 0.4) ) );
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(-5, 0, 0), 
				Colour(0.4, 0.4, 0.4) ) );
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 0), 
				Colour(0.6, 0.6, 0.6) ) );

	Material mat_earth( Colour(0.0, 0.0, 0.0), Colour(0.3, 0.3, 0.9), 
			Colour(0.316228, 0.316228, 0.316228), 
			52.8, 0.3, 0.0, 0.0  );
	Material mat_cloud( Colour(0.0, 0.0, 0.0), Colour(0.3, 0.3, 0.9), 
			Colour(0.316228, 0.316228, 0.316228), 
			52.8, 0.3, 0.0, 0.0  );
	Material mat_light( Colour(0.0, 0.0, 0.0), Colour(0.3, 0.3, 0.9), 
			Colour(0.616228, 0.616228, 0.616228), 
			52.8, 0.3, 0.0, 0.0  );
	// Add a unit square into the scene with material mat.
	SceneDagNode* earth = raytracer.addObject( new UnitSphere(), &mat_earth );
	SceneDagNode* earth_cloud = raytracer.addObject( new UnitSphere(), &mat_cloud );
	SceneDagNode* earth_light = raytracer.addObject( new UnitSphere(), &mat_light );
	
    Vector3D dis_factor(0, 3, 0);
    raytracer.translate(earth, dis_factor);	
    raytracer.translate(earth_light, (2.0*dis_factor));	
    // Add a texture; 
    Texture earthtex;
    earthtex.readImage("./Textures/earthmap1k.bmp");
    earth->mat->texture = &earthtex;  
    earth->mat->texture_ind = true; 
    // Add a texture; 
   Texture earthcloud;
   earthcloud.readImage("./Textures/earthcloudmap.bmp");
   earth_cloud->mat->texture = &earthcloud;  
   earth_cloud->mat->texture_ind = true; 
   // Add a texture; 
   Texture earthlight;
   earthlight.readImage("./Textures/earthlights1k.bmp");
   earth_light->mat->texture = &earthlight;  
   earth_light->mat->texture_ind = true; 
    

    Vector3D up(0, 0, -1);
	double fov = 30;
	Point3D eye(10, 3, 0);
	Vector3D view(-1, 0, 0);
    raytracer.render(width, height, eye, view, up, fov, "earth_1.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	double fov2 = 60;
    Vector3D up2(-1, 0, 0);
	Point3D eye2(0, 3, -5);
	Vector3D view2(0, 0, 1);
	raytracer.render(width, height, eye2, view2, up2, fov2, "earth_2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
#endif	
#ifdef SCENE5
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(-6, 10, 0), 
				Colour(0.3, 0.3, 0.3) ) );
	raytracer.addLightSource( new PointLight(Point3D(4, 18, 2), 
				Colour(0.3, 0.3, 0.3) ) );

	Material mat_galaxy( Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0), 
			Colour(0.016228, 0.016228, 0.016228), 
			1.0, 0.0, 0.0, 0.0  );
	Material mirror2(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(0.8, 0.8, 0.8),
		32.8, 1.0, 0.0, 0.0);
	mirror2.glossy_ind = true;
	// Defines a material for shading.
	Material gold_floor(Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
		Colour(0.628281, 0.555802, 0.366065),
		12.8, 0.6, 0.0, 0.0);
	Material black_floor(Colour(0.3, 0.3, 0.3), Colour(0.3, 0.3, 0.3),
		Colour(0.628281, 0.555802, 0.366065),
		12.8, 1.0, 0.0, 0.0);
	gold_floor.glossy_ind = true;
	black_floor.glossy_ind = true;
	// Add a unit square into the scene with material mat.
	SceneDagNode* plane = raytracer.addObject(new UnitSquare(), &black_floor);
	// SceneDagNode* plane_2 = raytracer.addObject( new UnitSquare(), &glass );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &mirror2 );
	SceneDagNode* sphere_2 = raytracer.addObject( new UnitSphere(), &mirror );
	SceneDagNode* sphere_3 = raytracer.addObject( new UnitSphere(), &glass );
	SceneDagNode* cube = raytracer.addObject( new UnitCube(), &mirror );
	SceneDagNode* galaxy = raytracer.addObject( new EnvirSphere(), &mat_galaxy );
    
    plane->mat->texture_ind = false; 
    sphere->mat->texture_ind = false; 
    sphere_2->mat->texture_ind = false; 
    sphere_3->mat->normal_ind = false; 
    plane->mat->normal_ind = false; 
    sphere->mat->normal_ind = false; 
    sphere_2->mat->normal_ind = false; 
    sphere_3->mat->normal_ind = false; 
	sphere->mat->glossy_ind = true;
    
    Texture galaxytex;
    galaxytex.readImage("./Textures/skysmall.bmp");
    galaxy->mat->texture = &galaxytex;  
    galaxy->mat->texture_ind = true; 
	
	// Apply some transformations to the unit square.
	double factor1[3] = { 3.0, 3.0, 3.0 };
	double factor2[3] = { 8.0, 8.0, 8.0 };
	double factor3[3] = { 0.6, 0.6, 0.6 };
	raytracer.translate(cube, Vector3D(-2.5, 0.5, -2.5));	
	raytracer.translate(sphere, Vector3D(1, 2.5, -1));	
	raytracer.translate(sphere_2, Vector3D(-1, 3, 1));	
	raytracer.translate(sphere_3, Vector3D(0, 4, 0));	
    raytracer.scale2(galaxy, Point3D(0, 0, 0), 150);

	raytracer.rotate(plane, 'x', -90); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	raytracer.scale(cube, Point3D(0, 0, 0), factor3);
	
   //  raytracer.translate(plane_2, Vector3D(0, 5, 0));	
   //  raytracer.rotate(plane_2, 'x', -90); 
   //  raytracer.scale(plane_2, Point3D(0, 0, 0), factor1);
	double fov3 = 60;
	Vector3D up3(0, 1, 0);
	Point3D eye3(0, 3, -7);
	Vector3D view3(0, 0, 1);
	raytracer.render(width, height, eye3, view3, up3, fov3, "scene5_view3.bmp");
	time(&timer2);
	double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;

	timer1 = timer2;
	double fov2 = 60;
	Vector3D up2(0, 1, 3);
	Point3D eye2(0, 9, -3);
	Vector3D view2(0, -3, 1);
	raytracer.render(width, height, eye2, view2, up2, fov2, "scene5_view2.bmp");
	time(&timer2);
	seconds = difftime(timer2, timer1);
	std::cout << "third image rendering finished. It takes " << seconds << " seconds." << std::endl;

	timer1 = timer2;
	Vector3D up(0, 0, 1);
	double fov = 60;
	Point3D eye(0, 9, 0);
	Vector3D view(0, -1, 0);
	raytracer.render(width, height, eye, view, up, fov, "scene5_view1.bmp");
	time(&timer2);
	// Render it from a different point of view.
	seconds = difftime(timer2, timer1);
	std::cout << "Second image rendering finished. It takes " << seconds << " seconds." << std::endl;

#endif	
#ifdef SCENE6
	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, -10, 0), 
				Colour(0.8, 0.8, 0.8) ) );
    raytracer.addLightSource( new PointLight(Point3D(0, -30, 30), 
      			Colour(0.6, 0.6, 0.6) ) );
    raytracer.addLightSource( new PointLight(Point3D(-5, 5, 0), 
      			Colour(0.4, 0.4, 0.4) ) );
    raytracer.addLightSource( new PointLight(Point3D(5, 5, 0), 
      			Colour(0.4, 0.4, 0.4) ) );
   //  raytracer.addLightSource( new PointLight(Point3D(0, 80, 0), 
   //    			Colour(0.8, 0.8, 0.8) ) );

	Material mat_earth( Colour(0.3, 0.3, 0.3), Colour(0.6, 0.6, 0.6), 
			Colour(0.316228, 0.316228, 0.616228), 
			12.8, 0.0, 0.0, 0.0  );
	Material mat_moon( Colour(0.3, 0.3, 0.3), Colour(0.6, 0.6, 0.6), 
			Colour(0.316228, 0.316228, 0.616228), 
			12.8, 0.0, 0.0, 0.0  );
	Material mat_galaxy( Colour(0.3, 0.3, 0.3), Colour(0.1, 0.1, 0.1), 
			Colour(0.016228, 0.016228, 0.016228), 
			1.0, 0.0, 0.0, 0.0  );
	// Add a unit square into the scene with material mat.
	SceneDagNode* earth = raytracer.addObject( new UnitSphere(), &mat_earth );
 	SceneDagNode* moon = raytracer.addObject( new UnitSphere(), &mat_moon );
	SceneDagNode* galaxy = raytracer.addObject( new EnvirSphere(), &mat_galaxy );
	
	// Apply some transformations to the unit square.
    // double radius[8] = {5, 12, 13, 7, 142, 120, 51, 49};
    double radius[8] = {40, 45, 53, 40, 162, 120, 81, 69};
    double distance[8] = {20, 200, 440, 700, 1300, 2000, 2500, 3000};
    // double rad_factor[3] = {0.01, 0.01, 0.01};
    Vector3D dis_factor(0.000, 0.02, 0.0);
    double rad_factor = 0.03;
    // 	double radius[8] = {1, 1, 1, 1, 1, 1, 1, 1};
	// double distance[8] = {1, 3, 5, 7, 9, 11, 13};
    //	Vector3D dis_factor(0.0, 1.0, 0.0);
	
    // Add a texture; 
    Texture earthtex;
    earthtex.readImage("./Textures/earthmap1k.bmp");
    earth->mat->texture = &earthtex;  
    earth->mat->texture_ind = true; 
    // Add a texture; 
    Texture earthnormal;
    earthnormal.readImage("./Textures/earthnormal.bmp");
    earth->mat->normalmap = &earthnormal;  
    earth->mat->normal_ind = true; 
    earth->mat->motion_ind = true; 
    earth->mat->motion_speed = 0.01; 
    earth->mat->motion_direction = Vector3D(1, 1, 1); 
    // Add a texture; 
    Texture moontex;
    moontex.readImage("./Textures/moonmap1k.bmp");
    moon->mat->texture = &moontex;  
    moon->mat->texture_ind = true; 
    // Add a texture; 
    Texture moonnormal;
    moonnormal.readImage("./Textures/moonnormal.bmp");
    moon->mat->normalmap = &moonnormal;  
    moon->mat->normal_ind = true; 
    // Add a texture; 
    Texture galaxytex;
    galaxytex.readImage("./Textures/starfield.bmp");
    galaxy->mat->texture = &galaxytex;  
    galaxy->mat->texture_ind = true; 
    
    raytracer.translate(moon, Vector3D(1, 3, 2));	
    raytracer.translate(galaxy, (distance[4]*dis_factor));	
   
    raytracer.scale2(earth, Point3D(0, 0, 0), rad_factor*radius[2]);
    raytracer.scale2(moon, Point3D(0, 0, 0), 0.8);
    raytracer.scale2(galaxy, Point3D(0, 0, 0), 50);
    
    raytracer.rotate(earth, 'z', 90); 
    raytracer.rotate(earth, 'x', 45); 
    // raytracer.rotate(earth, 'y', 90); 
	
    Vector3D up(-1, 0, 0);
	double fov = 60;
	Point3D eye(0, 0, 4);
	Vector3D view(0, 0, -1);
    
    raytracer.render(width, height, eye, view, up, fov, "scene6_view1.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	double fov2 = 60;
    Vector3D up2(0, 0, 1);
	Point3D eye2(-8, 2, 0);
	Vector3D view2(1, 0, 0);
	raytracer.render(width, height, eye2, view2, up2, fov2, "scene6_view2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	double fov3 = 60;
    Vector3D up3(0, 0, 1);
	Point3D eye3(3, 3, 2);
	Vector3D view3(-1, 0, 0);
	raytracer.render(width, height, eye3, view3, up3, fov3, "scene6_view3.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	double fov4 = 60;
    Vector3D up4(-1, 0, 0);
	Point3D eye4(0, 0, -3);
	Vector3D view4(0, 0, 1);
	raytracer.render(width, height, eye4, view4, up4, fov4, "scene6_view4.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
#endif

#ifdef SCENE7

    MeshObject teapot_mesh("./obj/bunny.obj");
    teapot_mesh.normalizeVectorCoord();
    printf("VERTEX\n");
    dumpvectorf(teapot_mesh.getVerts());
    printf("FACES\n");
    dumpvectorfaces(teapot_mesh.getFaces());
    printf("NORMALS\n");
    dumpvectorf(teapot_mesh.getNormals());
    printf("TEXCOORDS\n");
    dumpvectorf(teapot_mesh.getTexCoords());
    printf("NUMBER OF VERTEX:%d\n", teapot_mesh.getVerts().size());
    printf("NUMBER OF FACE:%d\n", teapot_mesh.getFaces().size()); 

	raytracer.addLightSource( new PointLight(Point3D(0, 0, 15), 
				Colour(0.6, 0.6, 0.6) ) );
	raytracer.addLightSource( new PointLight(Point3D(4, 4, 4), 
				Colour(0.4, 0.4, 0.4) ) );

	// Add a unit square into the scene with material mat.
	SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
	SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold );
	SceneDagNode* cube = raytracer.addObject( new GeneralObject(), &gold );
    	
	raytracer.rotate(cube, 'x', -45); 
	raytracer.rotate(cube, 'y', 45); 
	raytracer.rotate(cube, 'z', 45); 
   
    cube->obj->setMesh(&teapot_mesh);
    plane->mat->texture_ind = false; 
    cube->mat->texture_ind = false; 
    sphere->mat->texture_ind = false; 
    plane->mat->normal_ind = false; 
    cube->mat->normal_ind = false; 
    sphere->mat->normal_ind = false; 
	
    // Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor2[3] = { 12.0, 12.0, 12.0 };
	double factor3[3] = { 3.0, 3.0, 3.0 };
	
    raytracer.translate(sphere, Vector3D(0, 0, -5));	
	raytracer.rotate(sphere, 'x', -45); 
	raytracer.rotate(sphere, 'z', 45); 
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

    raytracer.translate(plane, Vector3D(0, 0, -7));	
	raytracer.rotate(plane, 'z', 45); 
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);
	raytracer.scale(cube, Point3D(0, 0, 0), factor3);
	
	Vector3D up(0, 1, 0);
	double fov = 60;
	Point3D eye(0, 0, 14);
	Vector3D view(0, 0, -9);
    raytracer.render(width, height, eye, view, up, fov, "scene7_view1.bmp");
	time(&timer2);
	// Render it from a different point of view.
    double seconds = difftime(timer2, timer1);
    std::cout<<"First image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;
    timer1 = timer2;
	Vector3D up2(0, 0, 1);
	Point3D eye2(-5, 0, 0);
	Vector3D view2(15, 0, 0);
	raytracer.render(width, height, eye2, view2, up2, fov, "scene7_view2.bmp");
	time(&timer2);
    seconds = difftime(timer2, timer1);
    std::cout<<"Second image rendering finished. It takes "<<seconds<<" seconds."<<std::endl;

#endif

#ifdef SCENE8
	// Defines a point light source.
	raytracer.addLightSource(new PointLight(Point3D(0, 0, 5.8),
		Colour(1.2, 1.2, 1.2)));
	Material cornell_gray(Colour(0.0, 0.0, 0.0), Colour(0.6, 0.6, 0.7),
		Colour(0.316228, 0.316228, 0.316228),
		12.8, 0.0, 0.0, 0.0);
	Material cornell_red(Colour(0.0, 0.0, 0.0), Colour(1.0, 0.0, 0.4),
		Colour(0.316228, 0.316228, 0.316228),
		12.8, 0.0, 0.0, 0.0);
	Material cornell_yellow(Colour(0.0, 0.0, 0.0), Colour(0.7, 0.9, 0.1),
		Colour(0.316228, 0.316228, 0.316228),
		12.8, 0.0, 0.0, 0.0);
	Material cornell_blue(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.4, 1.0),
		Colour(0.316228, 0.316228, 0.316228),
		12.8, 0.0, 0.0, 0.0);
	Material cornell_mirror(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(0.8, 0.8, 0.8),
		52.8, 1.0, 0.0, 0.0);
	Material cornell_glass(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(1.0, 1.0, 1.0),
		52.8, 0.1, 1.0, 1.25);
	Material cornell_glossy(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(0.8, 0.8, 0.8),
		52.8, 1.0, 0.0, 1.4);
	cornell_glossy.glossy_ind = true;
	// Add a unit square into the scene with material mat.
	SceneDagNode* plane_right = raytracer.addObject(new UnitSquare(), &cornell_red);
	SceneDagNode* plane_up = raytracer.addObject(new UnitSquare(), &cornell_gray);
	SceneDagNode* plane_left = raytracer.addObject(new UnitSquare(), &cornell_blue);
	SceneDagNode* plane_down = raytracer.addObject(new UnitSquare(), &cornell_gray);
	SceneDagNode* plane_front = raytracer.addObject(new UnitSquare(), &cornell_gray);

	SceneDagNode* sphere_mirror = raytracer.addObject(new UnitSphere(), &cornell_mirror);
	SceneDagNode* sphere_glossy = raytracer.addObject(new UnitSphere(), &cornell_glass);
	// SceneDagNode* cone = raytracer.addObject(new UnitCone(), &cornell_mirror);
	SceneDagNode* cube = raytracer.addObject(new UnitCube(), &cornell_mirror);

	// Apply some transformations to the unit square.
	double factor1[3] = { 2.0, 2.0, 2.0 };
	double factor2[3] = { 12.0, 12.0, 12.0 };
	double factor3[3] = { 1.2, 1.2, 3.0 };

	raytracer.translate(sphere_mirror, Vector3D(-3, -3, -4));
	// raytracer.translate(sphere_glass, Vector3D(4.5, 0, -5));
	raytracer.translate(sphere_glossy, Vector3D(1, 4, -4));
	// raytracer.translate(cone, Vector3D(2, -2, -6));
	raytracer.translate(cube, Vector3D(3, 0, -5));
	raytracer.rotate(cube, 'z', 30);

	raytracer.scale(sphere_mirror, Point3D(0, 0, 0), factor1);
	raytracer.scale(sphere_glossy, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane_down, Vector3D(0, 0, -6));
	raytracer.scale(plane_down, Point3D(0, 0, 0), factor2);
	
	raytracer.translate(plane_up, Vector3D(0, 0, 6));
	raytracer.rotate(plane_up, 'y', 180);
	raytracer.scale(plane_up, Point3D(0, 0, 0), factor2);
	
	raytracer.translate(plane_front, Vector3D(-6, 0, 0));
	raytracer.rotate(plane_front, 'y', 90);
	raytracer.scale(plane_front, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_left, Vector3D(0, 6, 0));
	raytracer.rotate(plane_left, 'x', 90);
	raytracer.scale(plane_left, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_right, Vector3D(0, -6, 0));
	raytracer.rotate(plane_right, 'x', -90);
	raytracer.scale(plane_right, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.
	// Camera parameters.
	Vector3D up(0, 0, 1);
	double fov = 60;
	Point3D eye(16, 0, 0);
	Vector3D view(-1, 0, 0);
	std::cout << "Generating photon map..." << std::endl;
	PhotonMap photon_map;
	raytracer.generatePhotonMap(&photon_map);
	//std::cout << "Rendering photon map..." << std::endl;
	//raytracer.renderWithPhoton(width, height, eye, view, up, fov, "scene8_global.bmp", &photon_map);
	raytracer.renderPhotonMap(width, height, eye, view, up, fov, "photon.bmp", &photon_map);
	std::cout << "Finished" << std::endl;
#endif

#ifdef SCENE9
	PhotonMap photon_map;
	int photon_num = 20000;
	// photon_map.initialPhotonMap(photon_num);
	// photon_map.displayPhotoMap();
	// Defines a point light source.
	raytracer.addLightSource(new PointLight(Point3D(0, 0, 5.8),
		Colour(1.2, 1.2, 1.2)));
	Material cornell_gray(Colour(0.0, 0.0, 0.0), Colour(0.5, 0.5, 0.5),
		Colour(0.0, 0.0, 0.0),
		12.8, 0.0, 0.0, 0.0);
	Material cornell_red(Colour(0.0, 0.0, 0.0), Colour(1.0, 0.0, 0.4),
		Colour(0.0, 0.0, 0.0),
		12.8, 0.0, 0.0, 0.0);
	Material cornell_yellow(Colour(0.0, 0.0, 0.0), Colour(0.7, 0.9, 0.0),
		Colour(0.116228, 0.116228, 0.116228),
		12.8, 0.0, 0.0, 0.0);
	Material cornell_blue(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.4, 1.0),
		Colour(0.0, 0.0, 0.0),
		12.8, 0.0, 0.0, 0.0);
	Material cornell_mirror(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(1.0, 1.0, 1.0),
		52.8, 1.0, 0.0, 0.0);
	Material cornell_glass(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(1.0, 1.0, 1.0),
		52.8, 0.1, 1.0, 1.28);
	Material cornell_glossy(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(0.8, 0.8, 0.8),
		52.8, 1.0, 0.0, 1.8);
	cornell_glossy.glossy_ind = true;
	// Add a unit square into the scene with material mat.
	SceneDagNode* plane_right = raytracer.addObject(new UnitSquare(), &cornell_red);
	SceneDagNode* plane_up = raytracer.addObject(new UnitSquare(), &cornell_gray);
	SceneDagNode* plane_left = raytracer.addObject(new UnitSquare(), &cornell_blue);
	SceneDagNode* plane_down = raytracer.addObject(new UnitSquare(), &cornell_gray);
	SceneDagNode* plane_front = raytracer.addObject(new UnitSquare(), &cornell_gray);

	SceneDagNode* sphere_mirror = raytracer.addObject(new UnitSphere(), &cornell_mirror);
	SceneDagNode* sphere_glossy = raytracer.addObject(new UnitSphere(), &cornell_glass);
	SceneDagNode* cube = raytracer.addObject(new UnitCube(), &cornell_mirror);

	// Apply some transformations to the unit square.
	double factor1[3] = { 2.0, 2.0, 2.0 };
	double factor2[3] = { 12.0, 12.0, 12.0 };
	double factor3[3] = { 1.2, 1.2, 3.0 };

	raytracer.translate(sphere_mirror, Vector3D(-3, -3, -4));
	// raytracer.translate(sphere_glass, Vector3D(4.5, 0, -5));
	raytracer.translate(sphere_glossy, Vector3D(1, 3.5, -4));
	// raytracer.translate(cone, Vector3D(2, -2, -6));
	raytracer.translate(cube, Vector3D(3, 0, -5));
	raytracer.rotate(cube, 'z', 30);

	//raytracer.translate(sphere_mirror, Vector3D(-3, -3, -4));
	//// raytracer.translate(sphere_glass, Vector3D(4.5, 0, -5));
	//raytracer.translate(sphere_glossy, Vector3D(-1, 2, -4));
	//raytracer.translate(cube, Vector3D(3, 4, -5));
	//raytracer.rotate(cube, 'z', 30);

	raytracer.scale(sphere_mirror, Point3D(0, 0, 0), factor1);
	raytracer.scale(sphere_glossy, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane_down, Vector3D(0, 0, -6));
	raytracer.scale(plane_down, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_up, Vector3D(0, 0, 6));
	raytracer.rotate(plane_up, 'y', 180);
	raytracer.scale(plane_up, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_front, Vector3D(-6, 0, 0));
	raytracer.rotate(plane_front, 'y', 90);
	raytracer.scale(plane_front, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_left, Vector3D(0, 6, 0));
	raytracer.rotate(plane_left, 'x', 90);
	raytracer.scale(plane_left, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_right, Vector3D(0, -6, 0));
	raytracer.rotate(plane_right, 'x', -90);
	raytracer.scale(plane_right, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.
	// Camera parameters.
	Vector3D up(0, 0, 1);
	double fov = 60;
	Point3D eye(16, 0, 0);
	Vector3D view(-1, 0, 0);
	std::cout << "Generating photon map..." << std::endl;
	raytracer.generatePhotonMap(&photon_map);
	std::cout << "Rendering photon map..." << std::endl;
	raytracer.renderWithPhoton(width, height, eye, view, up, fov, "scene8_global.bmp", &photon_map);
	raytracer.renderPhotonMap(width, height, eye, view, up, fov, "photon.bmp", &photon_map);
	std::cout << "Finished" << std::endl;
	// photon_map.displayPhotoMap();
#endif


    return 0;
}

