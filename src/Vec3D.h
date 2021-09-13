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

#ifndef    __VEC3D_H__
#define    __VEC3D_H__

#include <string>

namespace geom {

using namespace std;

class Vec3D {
public:
    // Constructor
    // 
    // Construct a zero vector (0,0,0).
    Vec3D(void); 
    // 
    // Construct a vector from its components.
    Vec3D(double tx, double ty, double tz);
    //
    // Construct a vector from another vector.
    Vec3D(const Vec3D& v0);
    //
    // Destructor
    ~Vec3D(void);
    
    // Operators
    //
    // Vector addition.
    Vec3D operator + (const Vec3D& v0) const;
    // 
    // Vector addition.
    Vec3D& operator += (const Vec3D& v0);
    //
    // Vector substraction.
    Vec3D operator - (const Vec3D& v0) const;
    //
    // Vector substraction.
    Vec3D& operator -= (const Vec3D& v0);
    //
    // Vector scaling: Multiplication.
    Vec3D operator * (double t) const;
    //
    // Vector scaling: Multiplication.
    Vec3D& operator *= (double t);
    //
    // Vector scaling: Dividion. 
    // Warning: if t == 0, then an infinite vector...
    Vec3D operator / (double t) const;
    //
    // Vector scaling: Dividion. 
    // Warning: if t == 0, then an infinite vector...
    Vec3D& operator /= (double t);
    //
    // Vector inner Product.
    double operator & (const Vec3D& v0) const;
    //
    // Vector outer Product.
    Vec3D operator * (const Vec3D& v0) const;
    //
    // Vector outer Product.
    Vec3D& operator *= (const Vec3D& v0);
    //
    // Assignment.
    Vec3D& operator = (const Vec3D& v0);
    //
    // Vector Order: == != < <= > >=
    // To determine if v and v0 are equal within threshold.
    //         The order of vector is defined as:
    //         1. if x < x0;
    //         2. if x = x0 and y < y0;
    //         3. if x = x0, y = y0 and z < z0;
    //         then v < v0.
    bool operator == (const Vec3D& v0) const;
    bool operator != (const Vec3D& v0) const;
    bool operator < (const Vec3D& v0) const;
    bool operator <= (const Vec3D& v0) const;
    bool operator > (const Vec3D& v0) const;
    bool operator >= (const Vec3D& v0) const;

    // Norms
    // 
    // Compute the square of the norm of the vector by x^2+y^2+z^2.
    double norm2(void) const;
    // 
    // Compute the square of the norm of the vector by x^2+y^2+z^2.    
    double norm2(const Vec3D& v0) const;
    // 
    // Compute the norm of the vector by sqrt(x^2+y^2+z^2).
    double norm(void) const;    
    // 
    // Compute the norm of the vector with another vector by sqrt(x^2+y^2+z^2).
    double norm(const Vec3D& v0) const;
    //
    // Get the norm '1' vector.
    // Warning: an infinite vector is obtained if norm is zero.
    Vec3D normI(void) const;
    //
    // To see if a vector is finite.
    bool isFinite(void) const;
    
    // Generate the plain text form of the vector.
    string tostr(void) const;

    // Threshold below which two double precision numbers are equal.
    static const double threshold; 

    // Vector componets
    double x;
    double y;
    double z;
};

}

#endif //  __VEC3D_H__
