/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		   light source classes

***********************************************************/

#include "util.h"
#include "photon.h"

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
	virtual void shade( Ray3D&, bool ) = 0;
	virtual Point3D get_position() const = 0; 
	virtual Colour get_col() const = 0;
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col) {}
	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular) {}
	void shade( Ray3D& ray, bool shadow );
	Point3D get_position() const { return _pos; }
	Colour get_col() const { return _col_diffuse;  }
private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};


// A point light is defined by its position in world space and its
// colour.
class SquareLight : public LightSource {
public:
	SquareLight(Point3D pos, Colour col) : _pos(pos), _col_ambient(col),
		_col_diffuse(col), _col_specular(col) {}
	SquareLight(Point3D pos, Colour ambient, Colour diffuse, Colour specular)
		: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse),
		_col_specular(specular) {}
	void shade(Ray3D& ray, bool shadow);
	Point3D get_position() const { return _pos; }
	Colour get_col() const { return _col_diffuse; }
	//void emit_photon(double num);
private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse;
	Colour _col_specular;
};
