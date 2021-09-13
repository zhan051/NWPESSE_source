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

#include "CEResult.h"
#include "RMCluster.h"
#include "RigidMolecule.h"
#include "RunStatus.h"
#include <boost/spirit/include/classic.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>

using namespace abcluster;
using namespace boost::spirit::classic;
using namespace std;

CEResult::CEResult(void)
{
    n = 0;
    energy = 0.00;
    x.clear();
    symbols.clear();
}

CEResult::~CEResult(void)
{
    /* Nothing to do */
}

int CEResult::Save(const string& fn) const
{
    SaveGJF(fn+".gjf"); 
    SaveXYZ(fn+".xyz");
    return RunStatus::Succeed;
}

int CEResult::Show(void) const
{
    for(int i = 0; i < n; ++i)
    {
        printf("%3s %15.8f %15.8f %15.8f\n", symbols[i].c_str(), x[i*3], x[i*3+1], x[i*3+2]);
    }
    return 0;
}

void CEResult::SetN(int tn)
{
    n = tn;
    if(tn <= 0)
    {
        RunStatus::ErrorTermination("Incorrect dimension of X (%d): A positive number is required.", tn);
    }    
    dof = n*3;
    x.clear();
    symbols.clear();

}

void CEResult::SetEnergy(double tenergy)
{
    energy = tenergy;
}

void CEResult::SetX(const vector<double>& tx)
{
    x.clear();
    if(tx.size() != dof)
    {
	RunStatus::ErrorTermination("The dimension of X (%d) is inconsistent with the cluster (%d).", tx.size(), dof);
    }
    else
    {
        x.assign(tx.begin(), tx.end());
    }
}

void CEResult::SetSymbols(const vector<string>& tsymbols)
{
    symbols.clear();
    if(tsymbols.size() != n)
    {
	RunStatus::ErrorTermination("The dimension of symbols (%d) is inconsistent with the cluster (%d).", 
		tsymbols.size(), n);
    }
    else
    {
        symbols.assign(tsymbols.begin(), tsymbols.end());
    }
}

double CEResult::GetEnergy(void) const
{
    return energy;
}

void CEResult::MoveSC(void)
{
    double xc = 0.0000;
    double yc = 0.0000;
    double zc = 0.0000;
    for(unsigned int i = 0; i < n; ++i)
    {
	xc += x[i*3];
	yc += x[i*3+1];
	zc += x[i*3+2];
    }
    xc /= n;
    yc /= n;
    zc /= n;
    for(unsigned int i = 0; i < n; ++i)
    {
	x[i*3] -= xc;
	x[i*3+1] -= yc;
	x[i*3+2] -= zc;
    }
}

int CEResult::SaveGJF(const string& fn) const
{
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL) 
    { 
	return RunStatus::ErrorTermination("Fail to write GJF file [ %s ]", fn.c_str());
    }

    fprintf(fd, "#B3LYP/6-31+G(d) OPT FREQ\n");
    fprintf(fd, "\n");
    fprintf(fd, "%d, Energy = %15.8f au\n", n, energy);
    fprintf(fd, "\n");
    fprintf(fd, "0 1\n");
    for(int i = 0; i < n; ++i)
    {
        fprintf(fd, "%3s %15.8f %15.8f %15.8f\n", symbols[i].c_str(), x[i*3], x[i*3+1], x[i*3+2]);	
    }
    fprintf(fd, "\n");
    fclose(fd);
    return RunStatus::Succeed;
}

int CEResult::SaveXYZ(const string& fn) const
{
    FILE* fd = fopen(fn.c_str(), "w");
    if(fd == NULL) 
    { 
	return RunStatus::ErrorTermination("Fail to write XYZ file [ %s ]", fn.c_str());
    }
    
    fprintf(fd, "%d\n", n);
    fprintf(fd, "Energy = %15.8f au\n", energy);
    for(int i = 0; i < n; ++i)
    {
        fprintf(fd, "%3s %15.8f %15.8f %15.8f\n", symbols[i].c_str(), x[i*3], x[i*3+1], x[i*3+2]);	
    }
    fclose(fd);
    return RunStatus::Succeed;
}

void CEResult::SolveExtCoords(const RMCluster& rmc, vector<ExtCoord>& eCoords) const
{        
    const vector<RigidMolecule>& mols = rmc.GetRigidMolecules();
    const int nmols = rmc.GetNRigidMolecules();
    eCoords.clear(); eCoords.assign(nmols, ExtCoord());
    for(int i = 0, p = 0; i < nmols; ++i)
    {
        mols[i].SolveExternalCoordinates(x.data()+3*p, eCoords[i]);
        p += mols[i].GetNAtoms();
    }
}



