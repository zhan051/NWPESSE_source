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

#include "RigidMolecule.h"
#include "lbfgs.h"
#include "tinyfun.h"
#include "RunStatus.h"
#include <boost/spirit/include/classic.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>

using namespace abcluster;
using namespace boost::spirit::classic;
using namespace std;

RigidMolecule::RigidMolecule(void)
{
    eCoord.phiZ1 = 0.0000;
    eCoord.phiY2 = 0.0000;
    eCoord.phiZ3 = 0.0000;
    eCoord.mX = 0.0000;
    eCoord.mY = 0.0000;
    eCoord.mZ = 0.0000;
}

RigidMolecule::~RigidMolecule(void)
{
    /* Nothing to do */
}

void RigidMolecule::LoadMolecule(const string& fn)
{
   // Open and load coordinates
   ifstream fd(fn.c_str(), ios::in);
   if(!fd.is_open()) 
   {   
       RunStatus::ErrorTermination("Fail to open file [ %s ] to load internal coordinates.", fn.c_str());
   }   
   string filecontent((istreambuf_iterator<char>(fd)), istreambuf_iterator<char>());
   fd.close();

   string tstr;
   istringstream iss(filecontent);
   
   // Number of atoms
   int tn;    
   const rule<> RULE_int = (*blank_p)>>(int_p)[assign(tn)]>>(*blank_p);
   getline(iss, tstr, '\n');
   parse_info<> parser = parse(tstr.c_str(), RULE_int);
   if(!parser.full)
   {	
       RunStatus::ErrorTermination("Incorrect XYZ file [ %s ] format near: \" %s \". Correct example: \" 12 \".", fn.c_str(), tstr.c_str());
   }
   if(tn <= 0)
   {
       RunStatus::ErrorTermination("Invalid number of atoms (%d): A positive integer is required.", tn);
   }
   natoms = tn;
   // Comments
   getline(iss, tstr, '\n');
   comments = tstr;
   // Coordinates
   string te;
   double tx, ty, tz;
   const rule<> RULE_xyz = (*blank_p)>>(*alnum_p)[assign(te)]>>
       (*blank_p)>>(real_p)[assign(tx)]>>
       (*blank_p)>>(real_p)[assign(ty)]>>
       (*blank_p)>>(real_p)[assign(tz)]>>
       (*blank_p);
   iCoord.clear(); iCoord.reserve(tn*6);
   symbols.clear(); symbols.reserve(tn*2);
   for(int i = 0; i < tn; ++i)
   {
       getline(iss, tstr, '\n');
       if(iss.fail()) { RunStatus::ErrorTermination("Insufficient number of atoms in XYZ file [ %s ]. It should have %d.", fn.c_str(), tn); }
       parser = parse(tstr.c_str(), RULE_xyz);
       if(!parser.full)
       {
	   RunStatus::ErrorTermination("Incorrect XYZ file [ %s ] format near: \" %s \". Correct example: \" C 3.4 5.2 7.9 \".", fn.c_str(), tstr.c_str());
       }
       iCoord.push_back(tx);
       iCoord.push_back(ty);
       iCoord.push_back(tz);
       symbols.push_back(te);
   }
   MoveToCenter();
}

int RigidMolecule::GetNAtoms(void) const
{
    return natoms;
}

const string& RigidMolecule::GetComments(void) const
{
    return comments;
}

const vector<double>& RigidMolecule::GetInternalCoordinates(void) const
{
    return iCoord;
}

const vector<string>& RigidMolecule::GetSymbols(void) const
{
    return symbols;
}

void RigidMolecule::SetExternalCoordinates(const ExtCoord& teCoord)
{
    eCoord = teCoord;
}

void RigidMolecule::GetExternalCoordinates(ExtCoord& teCoord) const
{
    teCoord = eCoord;
}

RigidMolecule& RigidMolecule::operator = (const RigidMolecule& rm)
{
    natoms = rm.natoms;
    comments = rm.comments;
    iCoord = rm.iCoord;       
    symbols = rm.symbols;
    eCoord = rm.eCoord;
    return *this;
}

void RigidMolecule::GetLabQuantity(vector<double>& labCoord, vector<double>& labDerivative) const
{
    labCoord.clear(); labCoord.assign(natoms*3, 0.);
    labDerivative.clear(); labDerivative.assign(natoms*18, 0.);

    const double c1 = cos(eCoord.phiZ1);
    const double s1 = sin(eCoord.phiZ1);
    const double c2 = cos(eCoord.phiY2);
    const double s2 = sin(eCoord.phiY2);
    const double c3 = cos(eCoord.phiZ3);
    const double s3 = sin(eCoord.phiZ3);
    const double RM[18] = {
         c1*c2*c3-s1*s3,  s1*c2*c3+c1*s3,  s2*c3, 
        -s1*c3-c1*c2*s3,  c1*c3-s1*c2*s3, -s2*s3, 
                 -c1*s2,          -s1*s2,     c2,
              -c1*s2*c3,       -s1*s2*c3,  c2*c3,
               c1*s2*s3,        s1*s2*s3, -c2*s3,
                 -c1*c2,          -s1*c2,    -s2
    };
    double tmp2[18] = {
        0., 0., 0., 1., 0., 0.,
        0., 0., 0., 0., 1., 0.,
        0., 0., 0., 0., 0., 1.,
    };

    for(int i = 0, i3 = 0, i18 = 0; i < natoms; ++i, i3 += 3, i18 += 18)
    {
	const double x0 = iCoord[i3];
	const double y0 = iCoord[i3+1];
	const double z0 = iCoord[i3+2];

        const double t1 = RM[3]*x0+RM[4]*y0+RM[5]*z0;
        const double t2 = RM[0]*x0+RM[1]*y0+RM[2]*z0;

	// Coordinates
	labCoord[i3]   = t2+eCoord.mX;
	labCoord[i3+1] = t1+eCoord.mY;
	labCoord[i3+2] = RM[6]*x0+RM[7]*y0+RM[8]*z0+eCoord.mZ;
	
	// Derivatives
        tmp2[0]  = -RM[1]*x0+RM[0]*y0;
	tmp2[1]  = RM[9]*x0+RM[10]*y0+RM[11]*z0;
	tmp2[2]  = t1;
        tmp2[6]  = -RM[4]*x0+RM[3]*y0;
	tmp2[7]  = RM[12]*x0+RM[13]*y0+RM[14]*z0;
        tmp2[8]  = -t2;
        tmp2[12]  = -RM[7]*x0+RM[6]*y0;
        tmp2[13]  = RM[15]*x0+RM[16]*y0+RM[17]*z0;
        //tmp2[14] = 0.0000;

        memcpy(&(labDerivative[i18]), tmp2, sizeof(double)*18);
    }
}

void RigidMolecule::GetLabQuantity(const double* x, vector<double>& labCoord,  vector<double>& labDerivative) const
{
    labCoord.clear(); labCoord.assign(natoms*3, 0.);
    labDerivative.clear(); labDerivative.assign(natoms*18, 0.);

    const double c1 = cos(x[0]);
    const double s1 = sin(x[0]);
    const double c2 = cos(x[1]);
    const double s2 = sin(x[1]);
    const double c3 = cos(x[2]);
    const double s3 = sin(x[2]);
    const double RM[18] = {
         c1*c2*c3-s1*s3,  s1*c2*c3+c1*s3,  s2*c3, 
        -s1*c3-c1*c2*s3,  c1*c3-s1*c2*s3, -s2*s3, 
                 -c1*s2,          -s1*s2,     c2,
              -c1*s2*c3,       -s1*s2*c3,  c2*c3,
               c1*s2*s3,        s1*s2*s3, -c2*s3,
                 -c1*c2,          -s1*c2,    -s2
    };
    double tmp2[18] = {
        0., 0., 0., 1., 0., 0.,
        0., 0., 0., 0., 1., 0.,
        0., 0., 0., 0., 0., 1.,
    };

    for(int i = 0, i3 = 0, i18 = 0; i < natoms; ++i, i3 += 3, i18 += 18)
    {
  const double x0 = iCoord[i3];
  const double y0 = iCoord[i3+1];
  const double z0 = iCoord[i3+2];

        const double t1 = RM[3]*x0+RM[4]*y0+RM[5]*z0;
        const double t2 = RM[0]*x0+RM[1]*y0+RM[2]*z0;

  // Coordinates
  labCoord[i3]   = t2+x[3];
  labCoord[i3+1] = t1+x[4];
  labCoord[i3+2] = RM[6]*x0+RM[7]*y0+RM[8]*z0+x[5];
  
  // Derivatives
        tmp2[0]  = -RM[1]*x0+RM[0]*y0;
  tmp2[1]  = RM[9]*x0+RM[10]*y0+RM[11]*z0;
  tmp2[2]  = t1;
        tmp2[6]  = -RM[4]*x0+RM[3]*y0;
  tmp2[7]  = RM[12]*x0+RM[13]*y0+RM[14]*z0;
        tmp2[8]  = -t2;
        tmp2[12]  = -RM[7]*x0+RM[6]*y0;
        tmp2[13]  = RM[15]*x0+RM[16]*y0+RM[17]*z0;
        //tmp2[14] = 0.0000;

        memcpy(&(labDerivative[i18]), tmp2, sizeof(double)*18);
    }
}

void RigidMolecule::MoveToCenter(void)
{
    double xc = 0.0000;
    double yc = 0.0000;
    double zc = 0.0000;

    for(int i = 0; i < natoms; ++i)
    {
        xc += iCoord[i*3];
        yc += iCoord[i*3+1];
        zc += iCoord[i*3+2];
    }
    xc /= natoms;
    yc /= natoms;
    zc /= natoms;
    
    for(int i = 0; i < natoms; ++i)
    {
        iCoord[i*3] -= xc;
        iCoord[i*3+1] -= yc;
        iCoord[i*3+2] -= zc;
    }
}

struct OptData {
    const RigidMolecule* mol;
    const double* labcoords;
};

void RigidMolecule::SolveExternalCoordinates(const double* labCoords, ExtCoord& eCoords) const
{
    // Initialization.    
    double x[6];
    for(int i = 0; i < 6; ++i) { x[i] = 0.; }
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    // Search solutions.
    double delta;
    OptData od; od.mol = this; od.labcoords = labCoords;
    const int ret = lbfgs(6, x, &delta, evaluate, progress, (void*)(&od), &param);
    eCoords = ExtCoord(x[0], x[1], x[2], x[3], x[4], x[5]);
}

int RigidMolecule::progress(void *instance, const lbfgsfloatval_t *x, 
        const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step, int n, int k, int ls)
{
#if 0
    printf("step: %5d energy: %15.8f gradient: %15.8f change: %15.8f\n", 
            k, fx, gnorm, step); 
#endif    
    return 0;
}

lbfgsfloatval_t RigidMolecule::evaluate(void *instance, 
        const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
        const lbfgsfloatval_t step)
{
    OptData* od = (OptData*)instance;
    const int num_atoms = od->mol->natoms;
    vector<double> coords;
    vector<double> derivatives;
    od->mol->GetLabQuantity(x, coords, derivatives);
    double delta = 0.;
    for(int i = 0; i < n; ++i) { g[i] = 0.; }
    for(int i = 0, i3 = 0, i18 = 0; i < num_atoms; ++i, i3 += 3, i18 += 18)
    {
        const double dx = coords[i3]-(od->labcoords[i3]);
        const double dy = coords[i3+1]-(od->labcoords[i3+1]);
        const double dz = coords[i3+2]-(od->labcoords[i3+2]);
        delta += dx*dx+dy*dy+dz*dz;
        for(int j = 0; j < 6; ++j) { g[j] += dx*derivatives[i18+j]; }
        for(int j = 0; j < 6; ++j) { g[j] += dy*derivatives[i18+6+j]; }
        for(int j = 0; j < 6; ++j) { g[j] += dz*derivatives[i18+12+j]; }
    }
    for(int i = 0; i < n; ++i) { g[i] *= 2; }
    return delta;
}