/*
 *      ABCluster
 *
 * Copyright (c) 2015 Jun Zhang
 * All rights reserved.
 *
 * The source code in this file is part of ABCluster, and is provided under 
 * a licence and may be used, copied, transmitted, or stored only in accord 
 * with that licence.
 *
 */

#include "Vec3D.h"
#include <limits>
#include <string>
#include <cmath>
#include <cstdio>

using namespace geom;
using namespace std;

// Threshold below which two double precision numbers are equal.
const double Vec3D::threshold = 1.00E-08;

Vec3D::Vec3D(void) : x(0.0), y(0.0), z(0.0) 
{ 
    /* Nothing to do */ 
}

Vec3D::Vec3D(double tx, double ty, double tz) : x(tx), y(ty), z(tz)
{ 
    /* Nothing to do */ 
}

Vec3D::Vec3D(const Vec3D& v0) : x(v0.x), y(v0.y), z(v0.z) 
{ 
    /* Nothing to do */ 
}

Vec3D::~Vec3D(void)
{ 
    /* Nothing to do */ 
}

//
Vec3D Vec3D::operator + (const Vec3D& v0) const    
{
    return Vec3D(x+v0.x, y+v0.y, z+v0.z); 
}

Vec3D& Vec3D::operator += (const Vec3D& v0)
{ 
    x += v0.x; y += v0.y; z += v0.z;
    return *this;
}

Vec3D Vec3D::operator - (const Vec3D& v0) const
{ 
    return Vec3D(x-v0.x, y-v0.y, z-v0.z); 
}

Vec3D& Vec3D::operator -= (const Vec3D& v0)
{ 
    x -= v0.x; y -= v0.y; z -= v0.z;
    return *this;
}

Vec3D Vec3D::operator * (double t) const
{ 
    return Vec3D(x*t, y*t, z*t); 
}

Vec3D& Vec3D::operator *= (double t)
{ 
    x *= t; y *= t; z *= t;
    return *this;
}

Vec3D Vec3D::operator / (double t) const
{ 
    return Vec3D(x/t, y/t, z/t); 
}

Vec3D& Vec3D::operator /= (double t)
{ 
    x /= t; y /= t; z /= t;
    return *this;
}

double Vec3D::operator & (const Vec3D& v0) const
{ 
    return x*v0.x+y*v0.y+z*v0.z; 
}

Vec3D Vec3D::operator * (const Vec3D& v0) const
{ 
    return Vec3D(y*v0.z-z*v0.y, z*v0.x-x*v0.z, x*v0.y-y*v0.x); 
}

Vec3D& Vec3D::operator *= (const Vec3D& v0)
{ 
    double tx, ty, tz;
    tx = y*v0.z-z*v0.y; ty = z*v0.x-x*v0.z; tz = x*v0.y-y*v0.x;
    x = tx; y = ty; z = tz;
    return *this;
}

Vec3D& Vec3D::operator = (const Vec3D& v0)
{
    x = v0.x; y = v0.y; z = v0.z;
    return *this;    
}

bool Vec3D::operator == (const Vec3D& v0) const
{
    return (fabs(x-v0.x) < threshold && fabs(y-v0.y) < threshold &&
		fabs(z-v0.z) < threshold);
}

bool Vec3D::operator != (const Vec3D& v0) const
{
    return !(*this == v0);
}

bool Vec3D::operator < (const Vec3D& v0) const
{
    if(x < v0.x) return true;
    if(fabs(x-v0.x) < threshold && y < v0.y) return true;
    if(fabs(x-v0.x) < threshold && fabs(y-v0.y) < threshold && z < v0.z) 
	return true;
    return false;
}

bool Vec3D::operator <= (const Vec3D& v0) const
{ 
    return (*this < v0) || (*this == v0); 
}

bool Vec3D::operator > (const Vec3D& v0) const
{ 
    return !(*this <= v0); 
}

bool Vec3D::operator >= (const Vec3D& v0) const
{
    return !(*this < v0); 
}

double Vec3D::norm2(void) const
{
    return x*x+y*y+z*z; 
}

double Vec3D::norm2(const Vec3D& v0) const
{
    return (*this-v0).norm2(); 
}

double Vec3D::norm(void) const
{ 
    return sqrt(this->norm2()); 
}

double Vec3D::norm(const Vec3D& v0) const    
{
    return sqrt(this->norm2(v0)); 
}

Vec3D Vec3D::normI(void) const
{
    double t = this->norm();
    return (*this)/t;
}

bool Vec3D::isFinite(void) const 
{
    double t = this->norm2();
    return (std::isfinite(x/t) && std::isfinite(y/t) && std::isfinite(z/t));
}

string Vec3D::tostr(void) const
{
    char buffer[64];
    sprintf(buffer, "%10.5f %10.5f %10.5f", x, y, z);
    return string(buffer);
}




