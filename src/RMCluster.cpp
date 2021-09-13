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

#include "RMCluster.h"
#include "RigidMolecule.h"
#include "RunStatus.h"
#include "tinyfun.h"
#include "lbfgs.h"
#include <boost/spirit/include/classic.hpp>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <cmath>
#include <cstring>
#include <cstdio>

using namespace abcluster;
using namespace boost::spirit::classic;
using namespace std;

RMCluster::RMCluster(void) : energy(0.0000)
{
    /* Nothing to do */ 
}

RMCluster::~RMCluster(void) 
{
    /* Nothing to do */ 
}

void RMCluster::LoadRMCluster(const string& fn)
{
   // Open and load coordinates
   ifstream fd(fn.c_str(), ios::in);
   if(!fd.is_open()) 
   {   
       RunStatus::ErrorTermination("Fail to open file [ %s ] to load Rigid molecular cluster.", fn.c_str());
   }   
   string filecontent((istreambuf_iterator<char>(fd)), istreambuf_iterator<char>());
   fd.close();

   string tstr;
   istringstream iss(filecontent);
   
   // Number of components
   int nc;    
   const rule<> RULE_int = (*blank_p)>>(int_p)[assign(nc)]>>(*blank_p);
   getline(iss, tstr, '\n');
   parse_info<> parser = parse(tstr.c_str(), RULE_int);
   if(!parser.full)
   {	
       RunStatus::ErrorTermination("Incorrect file [ %s ] format near: \" %s \". Correct example: \" 13 \"", fn.c_str(), tstr.c_str());
   }
   if(nc <= 0)
   {	
       RunStatus::ErrorTermination("Invalid number of components (%d): A positive integer is required.", nc);
   }
   // Each components
   string tfn;
   int tnc;
   const rule<> RULE_compt = (*blank_p)>>(*graph_p)[assign(tfn)]>>
       (*blank_p)>>(int_p)[assign(tnc)]>>
       (*blank_p);
   nrms = 0;
   components.clear(); components.reserve(nc*2);
   for(int i = 0; i < nc; ++i)
   {
       getline(iss, tstr, '\n');
       if(iss.fail()) { RunStatus::ErrorTermination("Insufficient number of clusters in  file [ %s ]. It should have %d.", fn.c_str(), nc); }
       parser = parse(tstr.c_str(), RULE_compt);
       if(!parser.full)
       {
	   RunStatus::ErrorTermination("Incorrect file [ %s ] format near: \" %s \". Correct example: \" tip4p.xyz  12 \".", fn.c_str(), tstr.c_str());
       }
       if(nc <= 0)
       {	
	   RunStatus::ErrorTermination("Invalid component number (%d): A positive integer is required.", tnc);
       }
       components.push_back(pair<string, int>(tfn, tnc));
       nrms += tnc;
   }
   // Load each components
   rms.clear(); rms.reserve(nrms*2);
   for(int i = 0; i < nc; ++i)
   {
       RigidMolecule rm;
       rm.LoadMolecule(components[i].first);
       for(int j = 0; j < components[i].second; ++j)
       {
           rms.push_back(rm);
       }
   } 
}

int RMCluster::CoarseOpt(const vector<int>& tfix_indices)
{
    // Initialization.
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    fix_indices = tfix_indices;
    // Optimization.
    int ret;
    const int dof = nrms*6;
    double* x = (double*)malloc(sizeof(double)*dof);
    GeteCoordsArray(x);
    lbfgsfloatval_t tenergy;
    ret = lbfgs(dof, x, &tenergy, evaluate, progress, (RMCluster*)this, &param);
    free(x);
    // Finish.
    return ret; 
}

double RMCluster::GetEnergy(void) const
{
    return energy;
}

void RMCluster::SetEnergy(double tenergy)
{
    energy = tenergy;
}


const vector<RigidMolecule>& RMCluster::GetRigidMolecules(void) const
{
    return rms;
}

int RMCluster::GetNRigidMolecules(void) const
{
    return nrms;
}

const vector<pair<string, int> >& RMCluster::GetComponents(void) const
{
    return components;
}

RMCluster& RMCluster::operator = (const RMCluster& rmc)
{
    nrms = rmc.nrms;
    rms = rmc.rms;
    energy = rmc.energy;
    components = rmc.components;
    return (*this);
}

void RMCluster::SaveAsXYZ(const string& fn) const
{
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL)
    {
        RunStatus::ErrorTermination("Cannot open file [ %s ] for saving the cluster.", fn.c_str());
    }
    int natoms = 0;
    for(int i = 0; i < nrms; ++i)
    {
        natoms += rms[i].GetNAtoms();
    }
    fprintf(fd, "%d\n", natoms);
    fprintf(fd, "Energy: %15.8f\n", energy);
    for(int i = 0; i < nrms; ++i)
    {
        const vector<string>& symbs = rms[i].GetSymbols();
        vector<double> labCoord;
        vector<double> labDerivative;
        rms[i].GetLabQuantity(labCoord, labDerivative);
        int nic = rms[i].GetNAtoms();
        for(int j = 0; j < nic; ++j)
        {
            fprintf(fd, "%5s %15.8f %15.8f %15.8f\n", symbs[j].c_str(), labCoord[3*j], labCoord[3*j+1], labCoord[3*j+2]);
        }
    }
    fclose(fd);
}

void RMCluster::GeteCoordsArray(double* eCoord) const
{
    for(int i = 0; i < nrms; ++i)
    {
        int i6 = i*6;
        ExtCoord ec;
        rms[i].GetExternalCoordinates(ec);
        eCoord[i6+0] = ec.phiZ1;
        eCoord[i6+1] = ec.phiY2;
        eCoord[i6+2] = ec.phiZ3;
        eCoord[i6+3] = ec.mX;
        eCoord[i6+4] = ec.mY;
        eCoord[i6+5] = ec.mZ;
    }
}

void RMCluster::SeteCoordsArray(const double* eCoord)
{
    for(int i = 0; i < nrms; ++i)
    {
        int i6 = i*6;
        ExtCoord ec(eCoord[i6+0], eCoord[i6+1], eCoord[i6+2], eCoord[i6+3], eCoord[i6+4], eCoord[i6+5]);
        rms[i].SetExternalCoordinates(ec);
    }
}


// Compute energy and grad.
void RMCluster::ComputEnergyGrad(double& tenergy, vector<double>& tgrad)
{
    const int nrms6 = nrms*6;

    tenergy = 0.;
    tgrad.clear(); tgrad.assign(nrms*6, 0.);
  
    for(int I = 0; I < nrms; ++I)
    {
	const int naI = rms[I].GetNAtoms();
	vector<double> CartI;
	vector<double> DerivI; 
	rms[I].GetLabQuantity(CartI, DerivI);
        const int I6 = I*6;
        double dEdIphiZ1, dEdIphiY2, dEdIphiZ3, dEdImX, dEdImY, dEdImZ;

        // Two-body interactions.
        for(int J = I+1; J < nrms; ++J)
	{
	    const int naJ = rms[J].GetNAtoms();
	    vector<double> CartJ;
   	    vector<double> DerivJ;
	    rms[J].GetLabQuantity(CartJ, DerivJ);
            const int J6 = J*6;
            double dEdJphiZ1, dEdJphiY2, dEdJphiZ3, dEdJmX, dEdJmY, dEdJmZ;

            double energyIJ;
            
            ComputeEnergyGradCore(
                    naI, CartI, DerivI, naJ, CartJ, DerivJ,
                    energyIJ,
                    dEdIphiZ1, dEdIphiY2, dEdIphiZ3, dEdImX, dEdImY, dEdImZ, 
                    dEdJphiZ1, dEdJphiY2, dEdJphiZ3, dEdJmX, dEdJmY, dEdJmZ
                    );

            tenergy += energyIJ;	
            tgrad[I6+0] += dEdIphiZ1;
            tgrad[I6+1] += dEdIphiY2;
            tgrad[I6+2] += dEdIphiZ3;
            tgrad[I6+3] += dEdImX;
            tgrad[I6+4] += dEdImY;
            tgrad[I6+5] += dEdImZ;
            tgrad[J6+0] += dEdJphiZ1;
            tgrad[J6+1] += dEdJphiY2;
            tgrad[J6+2] += dEdJphiZ3;
            tgrad[J6+3] += dEdJmX;
            tgrad[J6+4] += dEdJmY;
            tgrad[J6+5] += dEdJmZ;
	}
    }

    // One-body contraint terms.

}

void RMCluster::ComputeEnergyGradCore(
        int naI, const vector<double>& CartI, const vector<double>& DerivI,
        int naJ, const vector<double>& CartJ, const vector<double>& DerivJ,
        double& energyIJ,
        double& dEdIphiZ1, double& dEdIphiY2, double& dEdIphiZ3, double& dEdImX, double& dEdImY, double& dEdImZ,
        double& dEdJphiZ1, double& dEdJphiY2, double& dEdJphiZ3, double& dEdJmX, double& dEdJmY, double& dEdJmZ
        ) const
{
    energyIJ = 0.;
    dEdIphiZ1 = dEdIphiY2 = dEdIphiZ3 = dEdImX = dEdImY = dEdImZ = 0.;
    dEdJphiZ1 = dEdJphiY2 = dEdJphiZ3 = dEdJmX = dEdJmY = dEdJmZ = 0.;   

    for(int i = 0; i < naI; ++i)
    {
	const int i3 = i*3;
        const int i18 = i*18;
        for(int j = 0; j < naJ; ++j)
	{
	    const int j3 = j*3;
	    const int j18 = j*18;
	    const double xij = CartI[i3]-CartJ[j3];
	    const double yij = CartI[i3+1]-CartJ[j3+1];
	    const double zij = CartI[i3+2]-CartJ[j3+2];
	    const double rij2 = xij*xij+yij*yij+zij*zij;
	    const double rij  = sqrt(rij2);

            const double r0 = 2.0;
            const double k = 500;
	 
            // Energy
            energyIJ += (rij < r0) ? k*(rij-r0)*(rij-r0) : 0. ;
            // Gradient
            const double dedrij = (rij < r0) ? 2*k*(rij-r0) : 0. ;
            const double drijdxi = xij/rij; const double drijdxj = -drijdxi;
            const double drijdyi = yij/rij; const double drijdyj = -drijdyi;
            const double drijdzi = zij/rij; const double drijdzj = -drijdzi;

            dEdIphiZ1 += dedrij*(drijdxi*DerivI[i18+0]+drijdyi*DerivI[i18+6] +drijdzi*DerivI[i18+12]);
            dEdIphiY2 += dedrij*(drijdxi*DerivI[i18+1]+drijdyi*DerivI[i18+7] +drijdzi*DerivI[i18+13]);
            dEdIphiZ3 += dedrij*(drijdxi*DerivI[i18+2]+drijdyi*DerivI[i18+8] +drijdzi*DerivI[i18+14]);
            dEdImX    += dedrij*(drijdxi*DerivI[i18+3]+drijdyi*DerivI[i18+9] +drijdzi*DerivI[i18+15]);
            dEdImY    += dedrij*(drijdxi*DerivI[i18+4]+drijdyi*DerivI[i18+10]+drijdzi*DerivI[i18+16]);
            dEdImZ    += dedrij*(drijdxi*DerivI[i18+5]+drijdyi*DerivI[i18+11]+drijdzi*DerivI[i18+17]);

            dEdJphiZ1 += dedrij*(drijdxj*DerivJ[j18+0]+drijdyj*DerivJ[j18+6] +drijdzj*DerivJ[j18+12]);
            dEdJphiY2 += dedrij*(drijdxj*DerivJ[j18+1]+drijdyj*DerivJ[j18+7] +drijdzj*DerivJ[j18+13]);
            dEdJphiZ3 += dedrij*(drijdxj*DerivJ[j18+2]+drijdyj*DerivJ[j18+8] +drijdzj*DerivJ[j18+14]);
            dEdJmX    += dedrij*(drijdxj*DerivJ[j18+3]+drijdyj*DerivJ[j18+9] +drijdzj*DerivJ[j18+15]);
            dEdJmY    += dedrij*(drijdxj*DerivJ[j18+4]+drijdyj*DerivJ[j18+10]+drijdzj*DerivJ[j18+16]);
            dEdJmZ    += dedrij*(drijdxj*DerivJ[j18+5]+drijdyj*DerivJ[j18+11]+drijdzj*DerivJ[j18+17]);
	}
    }
}

int RMCluster::progress(void *instance, const lbfgsfloatval_t *x, 
        const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step, int n, int k, int ls)
{
#if 0
    printf("step: %5d energy: %10.7f gradient: %10.7f change: %10.7f\n", 
            k, fx, gnorm, step); 
#endif    
    return 0;
}

lbfgsfloatval_t RMCluster::evaluate(void *instance, 
        const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
        const lbfgsfloatval_t step)
{
    double tenergy;
    vector<double> tgrad;
    RMCluster* rmc = (RMCluster*)instance;
    rmc->SeteCoordsArray(x);
    rmc->ComputEnergyGrad(tenergy, tgrad);
    
    // Treat the fixed coordinates.
    for(auto i : rmc->fix_indices)
    {
        tgrad[i*6+0] = 0.;
        tgrad[i*6+1] = 0.;
        tgrad[i*6+2] = 0.;
        tgrad[i*6+3] = 0.;
        tgrad[i*6+4] = 0.;
        tgrad[i*6+5] = 0.;
    }

    memcpy(g, &(tgrad[0]), sizeof(double)*n);
    return tenergy;
}
