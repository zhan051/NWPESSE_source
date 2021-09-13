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

#ifndef    __STRUCTUREGENERATOR_H__
#define    __STRUCTUREGENERATOR_H__

#include "Vec3D.h"
#include "RMCluster.h"
#include "lbfgs.h"
#include <tuple>
#include <utility>
#include <vector>
#include <string>

namespace abcluster {

using namespace geom;    
using namespace std;

class StructureGenerator {
public:
    StructureGenerator(void);
    ~StructureGenerator(void); 

    void Initialization(const vector<string>& structuretypes, vector<string>& descriptions);
    void OneCluster(bool docoarseopt, RMCluster& rmc);

    void AddClusterForTraining(double energy, const vector<ExtCoord>& eCoords);


private:
    enum OpType { Pos, InBox, OutBox, InSphere, OnSurface, Micelle, Micelle4P, Micelle6P, Layer }; // + circle
    struct Op {
        OpType type;
        vector<double> real_params;
        vector<int> int_params;
    };
    const double Pi2 = 6.28318530718;
    const double Ang2Rad = 0.01745329251;

    enum CoordType { Random, Fix };    

private:
    bool initflag;
    vector<Op> ops;
    vector<CoordType> coord_types;
    void RandomCluster(RMCluster& rmc);   /// three algorithm in abc
    void Cross3PCluster(RMCluster& rmc) const;
    void Cross5PCluster(RMCluster& rmc) const;

    void DoPos(int nmols, const double* real_params, double* eCoords) const;
    void DoInBox(int nmols, const double* real_params, double* eCoords) const;
    void DoOutBox(int nmols, const double* real_params, double* eCoords) const;
    void DoInSphere(int nmols, const double* real_params, double* eCoords) const;
    void DoOnSurface(int nmols, const Vec3D& a0, const Vec3D& da, const vector<double>& points, double* eCoords) const;
    void DoMicelle(int nmols, const Vec3D& a0, const Vec3D& da, const double* real_params, double* eCoords) const;
    void DoLayer(int nmols, bool do_random, const Vec3D& a0, const Vec3D& da, const double* real_params, double* eCoords) const;
///  void Docircle()

private:    ////basic math cak\l
    void CalcConvexHull(const vector<double>& points, vector<tuple<int, int, int, Vec3D> >& faces) const;
    void CalcExtCoord(const Vec3D& a0, const Vec3D& an, const Vec3D& u0, const Vec3D& un, double* eCoords) const;
    void GenerateUniformSpherePoints(int n, vector<Vec3D>& points) const;
    void CalcMeanPlane(const vector<Vec3D>& points, vector<Vec3D>& proj_points) const;

    static int progress_CalcExtCoord(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
            const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);
    static lbfgsfloatval_t evaluate_CalcExtCoord(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
            const lbfgsfloatval_t step);

    static int progress_GenerateUniformSpherePoints(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
            const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);
    static lbfgsfloatval_t evaluate_GenerateUniformSpherePoints(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
            const lbfgsfloatval_t step);

    static int progress_CalcMeanPlane(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
            const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);
    static lbfgsfloatval_t evaluate_CalcMeanPlane(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
            const lbfgsfloatval_t step);

private:
    vector<pair<double, vector<double> > > training_set;
};

}

#endif //  __STRUCTUREGENERATOR_H__
