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
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"
#include "StructureGenerator.h"
#include "RMCluster.h"
#include "RigidMolecule.h"
#include "Vec3D.h"
#include "RunStatus.h"
#include "tinyfun.h"
#include "lbfgs.h"
#include <tuple>
#include <utility>
#include <algorithm>
#include <cmath>
#include <boost/spirit/include/classic.hpp>
#include <string>
#include <cstdio>

using namespace abcluster;
using namespace geom;
using namespace boost::spirit::classic;
using namespace std;

StructureGenerator::StructureGenerator(void) : initflag(false)
{
    training_set.clear();
    coord_types.clear();    
}

StructureGenerator::~StructureGenerator(void)
{
    /* Nothing to do */
}

void StructureGenerator::Initialization(const vector<string>& structuretypes, vector<string>& descriptions)  ////define grammar of feature
{
    if(initflag)
    {
        RunStatus::ErrorTermination("Class \"StructureGenerator\" has been initialized.");
    }
    initflag = true;


    // Parse the structure types.
    double td1, td2, td3, td4, td5, td6, td7, td8, td9, td10, td11, td12, td13, td14, td15;
    int tn1, tn2, tn3, tn4, tn5, tn6;
    const rule<> RULE_pos = str_p("pos")>>(*blank_p)>>(real_p)[assign(td1)]>>(*blank_p)>>(real_p)[assign(td2)]
        >>(*blank_p)>>(real_p)[assign(td3)]>>(*blank_p)>>(real_p)[assign(td4)]
        >>(*blank_p)>>(real_p)[assign(td5)]>>(*blank_p)>>(real_p)[assign(td6)];
    const rule<> RULE_inbox = str_p("inbox")>>(*blank_p)>>(real_p)[assign(td1)]>>(*blank_p)>>(real_p)[assign(td2)]
        >>(*blank_p)>>(real_p)[assign(td3)]>>(*blank_p)>>(real_p)[assign(td4)]
        >>(*blank_p)>>(real_p)[assign(td5)]>>(*blank_p)>>(real_p)[assign(td6)];
    const rule<> RULE_outbox = str_p("outbox")>>(*blank_p)>>(real_p)[assign(td1)]>>(*blank_p)>>(real_p)[assign(td2)]
        >>(*blank_p)>>(real_p)[assign(td3)]>>(*blank_p)>>(real_p)[assign(td4)]
        >>(*blank_p)>>(real_p)[assign(td5)]>>(*blank_p)>>(real_p)[assign(td6)]
        >>(*blank_p)>>(real_p)[assign(td7)]>>(*blank_p)>>(real_p)[assign(td8)]
        >>(*blank_p)>>(real_p)[assign(td9)]>>(*blank_p)>>(real_p)[assign(td10)]
        >>(*blank_p)>>(real_p)[assign(td11)]>>(*blank_p)>>(real_p)[assign(td12)];
    const rule<> RULE_insphere = str_p("insphere")>>(*blank_p)>>(real_p)[assign(td1)]>>(*blank_p)>>(real_p)[assign(td2)]
        >>(*blank_p)>>(real_p)[assign(td3)]>>(*blank_p)>>(real_p)[assign(td4)]>>(*blank_p)>>(real_p)[assign(td5)]>>(*blank_p)>>(real_p)[assign(td6)];
    // Surface treatment.
    const rule<> RULE_onsurface = str_p("onsurface")>>(*blank_p)>>(int_p)[assign(tn1)]
        >>(*blank_p)>>(int_p)[assign(tn2)]>>(*blank_p)>>(int_p)[assign(tn3)];

    const rule<> RULE_micelle = str_p("micelle")>>(*blank_p)>>(int_p)[assign(tn1)]>>(*blank_p)>>(int_p)[assign(tn2)]
        >>(*blank_p)>>(real_p)[assign(td1)]>>(*blank_p)>>(real_p)[assign(td2)]>>(*blank_p)>>(real_p)[assign(td3)]
        >>(*blank_p)>>(real_p)[assign(td4)]>>(*blank_p)>>(real_p)[assign(td5)]>>(*blank_p)>>(real_p)[assign(td6)];
    const rule<> RULE_micelle4P = str_p("micelle")>>(*blank_p)>>(int_p)[assign(tn1)]>>(*blank_p)>>(int_p)[assign(tn2)]
        >>(*blank_p)>>(int_p)[assign(tn3)]>>(*blank_p)>>(int_p)[assign(tn4)]
        >>(*blank_p)>>(real_p)[assign(td1)]>>(*blank_p)>>(real_p)[assign(td2)]>>(*blank_p)>>(real_p)[assign(td3)]
        >>(*blank_p)>>(real_p)[assign(td4)]>>(*blank_p)>>(real_p)[assign(td5)]>>(*blank_p)>>(real_p)[assign(td6)];
    const rule<> RULE_micelle6P = str_p("micelle")>>(*blank_p)>>(int_p)[assign(tn1)]>>(*blank_p)>>(int_p)[assign(tn2)]
        >>(*blank_p)>>(int_p)[assign(tn3)]>>(*blank_p)>>(int_p)[assign(tn4)]
        >>(*blank_p)>>(int_p)[assign(tn5)]>>(*blank_p)>>(int_p)[assign(tn6)]
        >>(*blank_p)>>(real_p)[assign(td1)]>>(*blank_p)>>(real_p)[assign(td2)]>>(*blank_p)>>(real_p)[assign(td3)]
        >>(*blank_p)>>(real_p)[assign(td4)]>>(*blank_p)>>(real_p)[assign(td5)]>>(*blank_p)>>(real_p)[assign(td6)];                

    const rule<> RULE_layer = str_p("layer")>>(*blank_p)>>(int_p)[assign(tn1)]>>(*blank_p)>>(int_p)[assign(tn2)]>>(*blank_p)>>(int_p)[assign(tn3)]
        >>(*blank_p)>>(real_p)[assign(td1)]>>(*blank_p)>>(real_p)[assign(td2)]>>(*blank_p)>>(real_p)[assign(td3)]
        >>(*blank_p)>>(real_p)[assign(td4)]>>(*blank_p)>>(real_p)[assign(td5)]>>(*blank_p)>>(real_p)[assign(td6)]
        >>(*blank_p)>>(real_p)[assign(td7)]>>(*blank_p)>>(real_p)[assign(td8)]>>(*blank_p)>>(real_p)[assign(td9)]
        >>(*blank_p)>>(real_p)[assign(td10)]>>(*blank_p)>>(real_p)[assign(td11)]>>(*blank_p)>>(real_p)[assign(td12)];        
// const rule<> Rule_circle = ...
    ops.clear();
    descriptions.clear();
    for(int i = 0; i < structuretypes.size(); ++i)
    {
        char buffer[1024];
        Op op0;
        if(parse(structuretypes[i].c_str(), RULE_pos).full)
        {
            op0.type = Pos;
            op0.real_params.clear();
            op0.real_params.push_back(td1); op0.real_params.push_back(td2); op0.real_params.push_back(td3);
            op0.real_params.push_back(td4); op0.real_params.push_back(td5); op0.real_params.push_back(td6);
            op0.int_params.clear();
            ops.push_back(op0);
            sprintf(buffer, "Center fixed at position (%.3f, %.3f, %.3f) with rotation (%.3f, %.3f, %.3f)",
                td1, td2, td3, td4, td5, td6);
            descriptions.push_back(buffer);            
            continue;
        }
        if(parse(structuretypes[i].c_str(), RULE_inbox).full)
        {
            op0.type = InBox;
            op0.real_params.clear();
            op0.real_params.push_back(td1); op0.real_params.push_back(td2); op0.real_params.push_back(td3);
            op0.real_params.push_back(td4); op0.real_params.push_back(td5); op0.real_params.push_back(td6);
            op0.int_params.clear();
            ops.push_back(op0);
            sprintf(buffer, "Distributed inside the box (%.3f, %.3f, %.3f) - (%.3f, %.3f, %.3f)",
                td1, td2, td3, td4, td5, td6);
            descriptions.push_back(buffer);
            continue;
        }
        if(parse(structuretypes[i].c_str(), RULE_outbox).full)
        {
            op0.type = OutBox;
            op0.real_params.clear();
            op0.real_params.push_back(td1); op0.real_params.push_back(td2); op0.real_params.push_back(td3);
            op0.real_params.push_back(td4); op0.real_params.push_back(td5); op0.real_params.push_back(td6);
            op0.real_params.push_back(td7); op0.real_params.push_back(td8); op0.real_params.push_back(td9);
            op0.real_params.push_back(td10); op0.real_params.push_back(td11); op0.real_params.push_back(td12);
            op0.int_params.clear();
            ops.push_back(op0);
            sprintf(buffer, "Distributed inside the box (%.3f, %.3f, %.3f) - (%.3f, %.3f, %.3f) but outside the box (%.3f, %.3f, %.3f) - (%.3f, %.3f, %.3f)",
                td7, td8, td9, td10, td11, td12, td1, td2, td3, td4, td5, td6);
            descriptions.push_back(buffer);
            continue;
        }
        if(parse(structuretypes[i].c_str(), RULE_insphere).full)
        {
            op0.type = InSphere;
            op0.real_params.clear();
            op0.real_params.push_back(td1); op0.real_params.push_back(td2); op0.real_params.push_back(td3);
            op0.real_params.push_back(td4); op0.real_params.push_back(td5); op0.real_params.push_back(td6);
            op0.int_params.clear();
            ops.push_back(op0);
            sprintf(buffer, "Distributed inside the sphere centered at (%.3f, %.3f, %.3f) with (ra, rb, rc) = (%.3f, %.3f %.3f)",
                td1, td2, td3, td4, td5, td6);
            descriptions.push_back(buffer);
            continue;
        }        
        if(parse(structuretypes[i].c_str(), RULE_onsurface).full)
        {
            op0.type = OnSurface;
            op0.real_params.clear();
            op0.int_params.clear();
            op0.int_params.push_back(tn1);
            op0.int_params.push_back(tn2);
            op0.int_params.push_back(tn3);
            ops.push_back(op0);
            sprintf(buffer, "Distributed on the surface of the %dth molecule, with axis %d - %d along the normal, %d attached on it",
                tn1, tn2, tn3, tn2);
            descriptions.push_back(buffer);            
            continue;
        }
        if(parse(structuretypes[i].c_str(), RULE_micelle).full)
        {
            op0.type = Micelle;
            op0.real_params.clear();
            op0.int_params.clear();
            op0.int_params.push_back(tn1);
            op0.int_params.push_back(tn2);
            op0.real_params.push_back(td1);
            op0.real_params.push_back(td2);
            op0.real_params.push_back(td3);
            op0.real_params.push_back(td4);
            op0.real_params.push_back(td5);
            op0.real_params.push_back(td6);
            ops.push_back(op0);
            sprintf(buffer, "Form micelle-like structure with axis %d - %d along the normal centered at (%.3f, %.3f, %.3f), core radius (%.3f, %.3f, %.3f)",
                tn1, tn2, td1, td2, td3, td4, td5, td6);
            descriptions.push_back(buffer);
            continue;
        }
        if(parse(structuretypes[i].c_str(), RULE_micelle4P).full)
        {
            op0.type = Micelle4P;
            op0.real_params.clear();
            op0.int_params.clear();
            op0.int_params.push_back(tn1);
            op0.int_params.push_back(tn2);
            op0.int_params.push_back(tn3);
            op0.int_params.push_back(tn4);
            op0.real_params.push_back(td1);
            op0.real_params.push_back(td2);
            op0.real_params.push_back(td3);
            op0.real_params.push_back(td4);
            op0.real_params.push_back(td5);
            op0.real_params.push_back(td6);
            ops.push_back(op0);
            sprintf(buffer, "Form micelle-like structure with axis (%d, %d) - (%d, %d) along the normal centered at (%.3f, %.3f, %.3f), core radius (%.3f, %.3f, %.3f)",
                tn1, tn2, tn3, tn4, td1, td2, td3, td4, td5, td6);
            descriptions.push_back(buffer);
            continue;
        }    
        if(parse(structuretypes[i].c_str(), RULE_micelle6P).full)
        {
            op0.type = Micelle6P;
            op0.real_params.clear();
            op0.int_params.clear();
            op0.int_params.push_back(tn1);
            op0.int_params.push_back(tn2);
            op0.int_params.push_back(tn3);
            op0.int_params.push_back(tn4);
            op0.int_params.push_back(tn5);
            op0.int_params.push_back(tn6);
            op0.real_params.push_back(td1);
            op0.real_params.push_back(td2);
            op0.real_params.push_back(td3);
            op0.real_params.push_back(td4);
            op0.real_params.push_back(td5);
            op0.real_params.push_back(td6);
            ops.push_back(op0);
            sprintf(buffer, "Form micelle-like structure with axis (%d, %d, %d) - (%d, %d, %d) along the normal centered at (%.3f, %.3f, %.3f), core radius (%.3f, %.3f, %.3f)",
                tn1, tn2, tn3, tn4, tn5, tn6, td1, td2, td3, td4, td5, td6);
            descriptions.push_back(buffer);
            continue;
        }            
        if(parse(structuretypes[i].c_str(), RULE_layer).full)
        {
            op0.type = Layer;
            op0.real_params.clear();
            op0.int_params.clear();
            op0.int_params.push_back(tn1);
            op0.int_params.push_back(tn2);
            op0.int_params.push_back(tn3);
            op0.real_params.push_back(td1);
            op0.real_params.push_back(td2);
            op0.real_params.push_back(td3);
            op0.real_params.push_back(td4);
            op0.real_params.push_back(td5);
            op0.real_params.push_back(td6);
            op0.real_params.push_back(td7);
            op0.real_params.push_back(td8);
            op0.real_params.push_back(td9);
            op0.real_params.push_back(td10);
            op0.real_params.push_back(td11);
            op0.real_params.push_back(td12);
            ops.push_back(op0);
            sprintf(buffer, "Form a %s layer with axis %d - %d along (%.3f, %.3f, %.3f) on the plane region (%.3f, %.3f, %.3f) to (%.3f, %.3f, %.3f) and (%.3f, %.3f, %.3f)",
                (tn3 == 0) ? "regular" : "random", tn1, tn2, 
                td1, td2, td3, td4, td5, td6, td7, td8, td9, td10, td11, td12);
            descriptions.push_back(buffer);
            continue;
        }
//// if()......
        RunStatus::ErrorTermination("Cannot parse \" %s \".", structuretypes[i].c_str());
    }
}

void StructureGenerator::OneCluster(bool docoarseopt, RMCluster& rmc)
{
static int scout_idx = 0;

    const int Nmax3P = 4;
    const double Pc3p = 0.4;
    const double dice = rand01();

    RandomCluster(rmc);
    if(training_set.size() >= Nmax3P && dice < Pc3p/2)
    {
        Cross3PCluster(rmc);
    }
    else
    {
        if(training_set.size() >= Nmax3P && dice < Pc3p)
        {
            Cross5PCluster(rmc);
        }
else
{
char buffer[1024];
sprintf(buffer, "scout%d.xyz",scout_idx);
rmc.SaveAsXYZ(string(buffer));
++scout_idx;
}
    }

    vector<int> fix_indices;
    const int nmols = coord_types.size()/6;
    for(int i = 0; i < nmols; ++i) { if(coord_types[i*6] == Fix) { fix_indices.push_back(i); } }
    if(docoarseopt) { rmc.CoarseOpt(fix_indices); }
}

void StructureGenerator::Cross3PCluster(RMCluster& rmc) const
{
    const vector<pair<string, int> >& components = rmc.GetComponents();
    const int nc = components.size();
    int tot_nmols = 0; for(int i = 0; i < nc; ++i) { tot_nmols += components[i].second; }
    vector<double> eCoord(tot_nmols*6, 0.);
    rmc.GeteCoordsArray(eCoord.data());

    const int npop = training_set.size();
    const int k1 = npop*rand01();
    const int k2 = npop*rand01();
    const int k3 = npop*rand01();
    const double p1 = training_set[k1].first;
    const double p2 = training_set[k2].first;
    const double p3 = training_set[k3].first;
    const double p = p1+p2+p3;
    const double d21 = (p2-p1)/p;
    const double d32 = (p3-p2)/p;
    const double d13 = (p1-p3)/p;
    for(int i = 0; i < tot_nmols*6; ++i)
    {
        if(coord_types[i] != Random) { continue; }
        eCoord[i] = (training_set[k1].second[i]+training_set[k2].second[i]+training_set[k3].second[i])/3
            +d21*(training_set[k1].second[i]-training_set[k2].second[i])
            +d32*(training_set[k2].second[i]-training_set[k3].second[i])
            +d13*(training_set[k3].second[i]-training_set[k1].second[i]);            
    }
    rmc.SeteCoordsArray(&(eCoord[0]));
//printf("3P: %d %d %d ", k1, k2, k3);
}

void StructureGenerator::Cross5PCluster(RMCluster& rmc) const
{
    const vector<pair<string, int> >& components = rmc.GetComponents();
    const int nc = components.size();
    int tot_nmols = 0; for(int i = 0; i < nc; ++i) { tot_nmols += components[i].second; }
    vector<double> eCoord(tot_nmols*6, 0.);
    rmc.GeteCoordsArray(eCoord.data());

    const int npop = training_set.size();
    const int k1 = npop*rand01();
    const int k2 = npop*rand01();
    const int k3 = npop*rand01();
    const int k4 = npop*rand01();
    const double length = rand01();
    for(int i = 0; i < tot_nmols*6; ++i)
    {
        if(coord_types[i] != Random) { continue; }
        eCoord[i] = training_set[0].second[i]+
            length*(training_set[k1].second[i]-training_set[k3].second[i]+training_set[k2].second[i]-training_set[k4].second[i]);
    }
    rmc.SeteCoordsArray(&(eCoord[0]));
//printf("5P: %d %d %d %d ", k1, k2, k3, k4);
}

void StructureGenerator::RandomCluster(RMCluster& rmc)
{
    const vector<pair<string, int> >& components = rmc.GetComponents();
    const int nc = components.size();
    int tot_nmols = 0; for(int i = 0; i < nc; ++i) { tot_nmols += components[i].second; }
    vector<double> eCoord(tot_nmols*6, 0.);
    coord_types.clear(); coord_types.assign(tot_nmols*6, Random);

    // Struture distribution for each component.
    for(int i = 0, p = 0; i < nc; ++i)
    {
        const int nmols = components[i].second;
        if(ops[i].type == Pos)
        {
            for(int j = 0; j < nmols*6; ++j) { coord_types[p+j] = Fix; } // This is one must be fixed!
            DoPos(nmols, ops[i].real_params.data(), eCoord.data()+p);
            p += nmols*6;
            continue;
        }
        if(ops[i].type == InBox)
        {
            DoInBox(nmols, ops[i].real_params.data(), eCoord.data()+p);
            p += nmols*6;
            continue;
        }
        if(ops[i].type == OutBox)
        {
            DoOutBox(nmols, ops[i].real_params.data(), eCoord.data()+p);
            p += nmols*6;
            continue;
        }
        if(ops[i].type == InSphere)
        {
            DoInSphere(nmols, ops[i].real_params.data(), eCoord.data()+p);
            p += nmols*6;
            continue;
        }        
        if(ops[i].type == OnSurface)
        {
            for(int j = 0; j < nmols*6; ++j) { coord_types[p+j] = Fix; } // This is one must be fixed!
            // Get the absolute coordinates of te host molecule.
            const int idx = ops[i].int_params[0]-1;
            if(idx >= p/6)
            {
                RunStatus::ErrorTermination("The absorbed molecules (%d) must come after the host molecule (%d).", p/6+1, idx+1);
            }            
            RigidMolecule rm(rmc.GetRigidMolecules()[idx]);
            rm.SetExternalCoordinates(ExtCoord(eCoord[idx*6], eCoord[idx*6+1], eCoord[idx*6+2], eCoord[idx*6+3], eCoord[idx*6+4], eCoord[idx*6+5]));
            vector<double> points, labDerivative;
            rm.GetLabQuantity(points, labDerivative);
            Vec3D a0(0., 0., 0.), da(0., 0., 0.);
            // Construct the unit axis.
            const vector<double>& ligand = rmc.GetRigidMolecules()[p/6].GetInternalCoordinates();
            const int iatom0 = ops[i].int_params[1]-1;
            const int iatom1 = ops[i].int_params[2]-1;                
            if(iatom0 < 0 || iatom0 >= ligand.size()/3 || iatom1 < 0 || iatom1 >= ligand.size()/3)
            {
                RunStatus::ErrorTermination("The atomic indices %d and %d are smaller than 1 or larger than the total number of atoms of this moleucle (%d).",
                    iatom0+1, iatom1+1, ligand.size()/3);
            }
            // TODO: extract coordinates.
            a0 = Vec3D(ligand[iatom0*3], ligand[iatom0*3+1], ligand[iatom0*3+2]);
            const Vec3D a1(ligand[iatom1*3], ligand[iatom1*3+1], ligand[iatom1*3+2]);
            da = a1-a0;

            DoOnSurface(nmols, a0, da, points, eCoord.data()+p);
            p += nmols*6;
            continue;
        }
        if(ops[i].type == Micelle)
        {
            for(int j = 0; j < nmols*6; ++j) { coord_types[p+j] = Fix; } // This is one must be fixed!
            Vec3D a0(0., 0., 0.), da(0., 0., 0.);
            // Construct the unit axis.
            const vector<double>& ligand = rmc.GetRigidMolecules()[p/6].GetInternalCoordinates();
            const int iatom0 = ops[i].int_params[0]-1;
            const int iatom1 = ops[i].int_params[1]-1;                
            if(iatom0 < 0 || iatom0 >= ligand.size()/3 || iatom1 < 0 || iatom1 >= ligand.size()/3)
            {
                RunStatus::ErrorTermination("The atomic indices %d and %d are smaller than 1 or larger than the total number of atoms of this moleucle (%d).",
                    iatom0+1, iatom1+1, ligand.size()/3);
            }
            a0 = Vec3D(ligand[iatom0*3], ligand[iatom0*3+1], ligand[iatom0*3+2]);
            const Vec3D a1(ligand[iatom1*3], ligand[iatom1*3+1], ligand[iatom1*3+2]);
            da = a1-a0;
            if(ops[i].real_params[3] <= 0 || ops[i].real_params[4] <= 0 || ops[i].real_params[5] <= 0)
            {
                RunStatus::ErrorTermination("A positive real number is needed for the internal radius of micelle but not %.3f, %.3f, %.3f..",
                    ops[i].real_params[3], ops[i].real_params[4], ops[i].real_params[5]);
            }
            DoMicelle(nmols, a0, da, ops[i].real_params.data(), eCoord.data()+p);
            p += nmols*6;
            continue;
        }
        if(ops[i].type == Micelle4P)
        {
            for(int j = 0; j < nmols*6; ++j) { coord_types[p+j] = Fix; } // This is one must be fixed!
            Vec3D a0(0., 0., 0.), da(0., 0., 0.);
            // Construct the unit axis.
            const vector<double>& ligand = rmc.GetRigidMolecules()[p/6].GetInternalCoordinates();
            const int iatom00 = ops[i].int_params[0]-1;
            const int iatom01 = ops[i].int_params[1]-1;
            const int iatom10 = ops[i].int_params[2]-1;                
            const int iatom11 = ops[i].int_params[3]-1;                
            if(iatom00 < 0 || iatom00 >= ligand.size()/3 || iatom01 < 0 || iatom01 >= ligand.size()/3
                 || iatom10 < 0 || iatom10 >= ligand.size()/3  || iatom11 < 0 || iatom11 >= ligand.size()/3)
            {
                RunStatus::ErrorTermination("The atomic indices %d, %d, %d, and %d are smaller than 1 or larger than the total number of atoms of this moleucle (%d).",
                    iatom00+1, iatom01+1, iatom10+1, iatom11+1, ligand.size()/3);
            }
            a0 = Vec3D((ligand[iatom00*3]+ligand[iatom01*3])/2, (ligand[iatom00*3+1]+ligand[iatom01*3+1])/2, (ligand[iatom00*3+2]+ligand[iatom01*3+2])/2);
            const Vec3D a1((ligand[iatom10*3]+ligand[iatom11*3])/2, (ligand[iatom10*3+1]+ligand[iatom11*3+1])/2, (ligand[iatom10*3+2]+ligand[iatom11*3+2])/2);
            da = a1-a0;
            if(ops[i].real_params[3] <= 0 || ops[i].real_params[4] <= 0 || ops[i].real_params[5] <= 0)
            {
                RunStatus::ErrorTermination("A positive real number is needed for the internal radius of micelle but not %.3f, %.3f, %.3f..",
                    ops[i].real_params[3], ops[i].real_params[4], ops[i].real_params[5]);
            }
            DoMicelle(nmols, a0, da, ops[i].real_params.data(), eCoord.data()+p);
            p += nmols*6;
            continue;
        }
        if(ops[i].type == Micelle6P)
        {
            for(int j = 0; j < nmols*6; ++j) { coord_types[p+j] = Fix; } // This is one must be fixed!
            Vec3D a0(0., 0., 0.), da(0., 0., 0.);
            // Construct the unit axis.
            const vector<double>& ligand = rmc.GetRigidMolecules()[p/6].GetInternalCoordinates();
            const int iatom00 = ops[i].int_params[0]-1;
            const int iatom01 = ops[i].int_params[1]-1;
            const int iatom02 = ops[i].int_params[2]-1;
            const int iatom10 = ops[i].int_params[3]-1;                
            const int iatom11 = ops[i].int_params[4]-1;                
            const int iatom12 = ops[i].int_params[5]-1;                
            if(iatom00 < 0 || iatom00 >= ligand.size()/3 || iatom01 < 0 || iatom01 >= ligand.size()/3 || iatom02 < 0 || iatom02 >= ligand.size()/3
                 || iatom10 < 0 || iatom10 >= ligand.size()/3  || iatom11 < 0 || iatom11 >= ligand.size()/3 || iatom12 < 0 || iatom12 >= ligand.size()/3)
            {
                RunStatus::ErrorTermination("The atomic indices %d, %d, %d, %d, %d, and %d are smaller than 1 or larger than the total number of atoms of this moleucle (%d).",
                    iatom00+1, iatom01+1, iatom02+1, iatom10+1, iatom11+1, iatom12+1, ligand.size()/3);
            }
            a0 = Vec3D((ligand[iatom00*3]+ligand[iatom01*3]+ligand[iatom02*3])/3, (ligand[iatom00*3+1]+ligand[iatom01*3+1]+ligand[iatom02*3+1])/3, (ligand[iatom00*3+2]+ligand[iatom00*3+2]+ligand[iatom02*3+2])/3);
            const Vec3D a1((ligand[iatom10*3]+ligand[iatom11*3]+ligand[iatom12*3])/3, (ligand[iatom10*3+1]+ligand[iatom11*3+1]+ligand[iatom12*3+1])/3, (ligand[iatom10*3+2]+ligand[iatom11*3+2]+ligand[iatom12*3+2])/3);
            da = a1-a0;
            if(ops[i].real_params[3] <= 0 || ops[i].real_params[4] <= 0 || ops[i].real_params[5] <= 0)
            {
                RunStatus::ErrorTermination("A positive real number is needed for the internal radius of micelle but not %.3f, %.3f, %.3f..",
                    ops[i].real_params[3], ops[i].real_params[4], ops[i].real_params[5]);
            }
            DoMicelle(nmols, a0, da, ops[i].real_params.data(), eCoord.data()+p);
            p += nmols*6;
            continue;
        }              
        if(ops[i].type == Layer)
        {
            for(int j = 0; j < nmols*6; ++j) { coord_types[p+j] = Fix; } // This is one must be fixed!            
            Vec3D a0(0., 0., 0.), da(0., 0., 0.);
            // Construct the unit axis.
            const vector<double>& ligand = rmc.GetRigidMolecules()[p/6].GetInternalCoordinates();
            const int iatom0 = ops[i].int_params[0]-1;
            const int iatom1 = ops[i].int_params[1]-1;                
            if(iatom0 < 0 || iatom0 >= ligand.size()/3 || iatom1 < 0 || iatom1 >= ligand.size()/3)
            {
                RunStatus::ErrorTermination("The atomic indices %d and %d are smaller than 1 or larger than the total number of atoms of this moleucle (%d).",
                    iatom0+1, iatom1+1, ligand.size()/3);
            }
            a0 = Vec3D(ligand[iatom0*3], ligand[iatom0*3+1], ligand[iatom0*3+2]);
            const Vec3D a1(ligand[iatom1*3], ligand[iatom1*3+1], ligand[iatom1*3+2]);
            da = a1-a0;
            if(ops[i].int_params[2] == 0) // Regular
            {            
                DoLayer(nmols, false, a0, da, ops[i].real_params.data(), eCoord.data()+p);
            }
            else // Random
            {
                DoLayer(nmols, true, a0, da, ops[i].real_params.data(), eCoord.data()+p);
            }
            p += nmols*6;
            continue;
        }
////   if(op[i].type = circle )   ...call do circles..
    }

    // Set Clusters
    rmc.SeteCoordsArray(&(eCoord[0]));
    rmc.SetEnergy(0.0000);
}

void StructureGenerator::AddClusterForTraining(double energy, const vector<ExtCoord>& eCoords)
{
    vector<double> coords;
    for(int i = 0; i < eCoords.size(); ++i)
    {
        coords.push_back(eCoords[i].phiZ1);
        coords.push_back(eCoords[i].phiY2);
        coords.push_back(eCoords[i].phiZ3);
        coords.push_back(eCoords[i].mX);
        coords.push_back(eCoords[i].mY);
        coords.push_back(eCoords[i].mZ);
    }
    const int MaxTrainingSetSize = 10;
    if(training_set.size() < MaxTrainingSetSize)
    {
        training_set.push_back(make_pair(energy, coords));
    }
    else
    {
        sort(training_set.begin(), training_set.end(), 
            [](const pair<double, vector<double> >& a, const pair<double, vector<double> >& b){ return a.first < b.first; });
        const int i = training_set.size()-1;
        if(energy < training_set[i].first)
        {
            training_set[i].first = energy;
            training_set[i].second = coords;
        }

//for(auto x: training_set) { printf("%15.8f ", x.first); } printf("\n");
    }
}

void StructureGenerator::DoPos(int nmols, const double* real_params, double* eCoords) const
{
    const double x0 = real_params[0];
    const double y0 = real_params[1];
    const double z0 = real_params[2];
    const double a0 = real_params[3];
    const double b0 = real_params[4];
    const double c0 = real_params[5];

    const double ep = -1.E-6;
    for(int i = 0, i6 = 0; i < nmols; ++i, i6 += 6)
    {
        eCoords[i6] = (a0 >= ep) ? a0*Ang2Rad : rand01()*Pi2;
        eCoords[i6+1] = (b0 >= ep) ? b0*Ang2Rad : rand01()*Pi2;
        eCoords[i6+2] = (c0 >= ep) ? c0*Ang2Rad : rand01()*Pi2;
        eCoords[i6+3] = x0;
        eCoords[i6+4] = y0;
        eCoords[i6+5] = z0;
    }
}

void StructureGenerator::DoInBox(int nmols, const double* real_params, double* eCoords) const
{
    const auto bx = minmax(real_params[0], real_params[3]);
    const auto by = minmax(real_params[1], real_params[4]);
    const auto bz = minmax(real_params[2], real_params[5]);
    const double x0 = bx.first; const double dx = bx.second-bx.first;
    const double y0 = by.first; const double dy = by.second-by.first;
    const double z0 = bz.first; const double dz = bz.second-bz.first;

    for(int i = 0, i6 = 0; i < nmols; ++i, i6 += 6)
    {
        eCoords[i6] = rand01()*Pi2;
        eCoords[i6+1] = rand01()*Pi2;
        eCoords[i6+2] = rand01()*Pi2;
        eCoords[i6+3] = x0+rand01()*dx;
        eCoords[i6+4] = y0+rand01()*dy;
        eCoords[i6+5] = z0+rand01()*dz;
    }
}

void StructureGenerator::DoOutBox(int nmols, const double* real_params, double* eCoords) const
{
    const auto box = minmax(real_params[0], real_params[3]);
    const auto boy = minmax(real_params[1], real_params[4]);
    const auto boz = minmax(real_params[2], real_params[5]);
    const auto bix = minmax(real_params[6], real_params[9]);
    const auto biy = minmax(real_params[7], real_params[10]);
    const auto biz = minmax(real_params[8], real_params[11]);
    
    const double x1 = bix.first;  const double dx1 = bix.second-bix.first;
    const double y1 = biy.first;  const double dy1 = biy.second-biy.first;
    const double z1 = biz.first;  const double dz1 = biz.second-biz.first;

    const double x2 = box.first;  const double dx2 = box.second-box.first;
    const double y2 = boy.first;  const double dy2 = boy.second-boy.first;
    const double z2 = boz.first;  const double dz2 = boz.second-boz.first;

    for(int i = 0, i6 = 0; i < nmols; ++i, i6 += 6)
    {
        eCoords[i6] = rand01()*Pi2;
        eCoords[i6+1] = rand01()*Pi2;
        eCoords[i6+2] = rand01()*Pi2;

        double x, y, z;
        const int Ntry = 1000;
        for(int j = 0; j < Ntry; ++j)
        {
            x = x1+rand01()*dx1;
            y = y1+rand01()*dy1;
            z = z1+rand01()*dz1;

            const double dlx = (x-x2)/dx2;
            const double dly = (y-y2)/dy2;
            const double dlz = (z-z2)/dz2;
            if(!(dlx >= 0. && dlx <= 1. && dly >= 0. && dly <= 1. && dlz >= 0. && dlz <= 1.)) { break; }
        }

        eCoords[i6+3] = x;
        eCoords[i6+4] = y;
        eCoords[i6+5] = z;
    }
}


void StructureGenerator::DoInSphere(int nmols, const double* real_params, double* eCoords) const
{
    const double x0 = real_params[0];
    const double y0 = real_params[1];
    const double z0 = real_params[2];
    const double ra = real_params[3];
    const double rb = real_params[4];
    const double rc = real_params[5];    
    for(int i = 0, i6 = 0; i < nmols; ++i, i6 += 6)
    {
        eCoords[i6] = rand01()*Pi2;
        eCoords[i6+1] = rand01()*Pi2;
        eCoords[i6+2] = rand01()*Pi2;
        
        const double phi = rand01()*Pi2;
        const double costheta = rand01()*2-1;
        const double u = pow(rand01(), 1./3.);
        const double theta = acos(costheta);
        const double dx = ra*u*sin(theta)*cos(phi);
        const double dy = rb*u*sin(theta)*sin(phi);
        const double dz = rc*u*cos(theta);

        eCoords[i6+3] = x0+dx;
        eCoords[i6+4] = y0+dy;
        eCoords[i6+5] = z0+dz;
    }
}

void StructureGenerator::DoOnSurface(int nmols, const Vec3D& a0, const Vec3D& da, 
    const vector<double>& points, double* eCoords) const
{
    // Get the convex hull of points.
    vector<tuple<int, int, int, Vec3D> > faces;
    CalcConvexHull(points, faces);

    const int nFaces = faces.size();
    for(int i = 0, i6 = 0; i < nmols; ++i, i6 += 6)
    {
        // Pick up a face randomly.
        auto face = faces[nFaces*rand01()];
        const int i0 = get<0>(face);
        const int i1 = get<1>(face);
        const int i2 = get<2>(face);            
        const Vec3D v0(points[i0*3], points[i0*3+1], points[i0*3+2]);
        const Vec3D v1(points[i1*3], points[i1*3+1], points[i1*3+2]);
        const Vec3D v2(points[i2*3], points[i2*3+1], points[i2*3+2]);        
        const Vec3D v01 = v0+(v1-v0)*rand01();
        const Vec3D v012 = v01+(v2-v01)*rand01();
        const double ds = 2.5;
        
        const Vec3D un(get<3>(face));                
        const Vec3D u0 = v012+un*ds;//+a0;

        if(da.norm2() < 1.E-6) // No well-defined norm.
        {
            eCoords[i6] = 0.;
            eCoords[i6+1] = 0.;
            eCoords[i6+2] = 0.;
            eCoords[i6+3] = u0.x;
            eCoords[i6+4] = u0.y;
            eCoords[i6+5] = u0.z;            
        }
        else // Rotate.
        {
            const Vec3D an = da.normI();
            CalcExtCoord(a0, an, u0, un, eCoords+i6);
        }
    }
}

void StructureGenerator::DoMicelle(int nmols, const Vec3D& a0, const Vec3D& da, const double* real_params, double* eCoords) const
{
    // Evenly points on unit sphere.
    vector<Vec3D> points;
    GenerateUniformSpherePoints(nmols, points);

    const Vec3D center(real_params[0], real_params[1], real_params[2]);
    const double ra = real_params[3];    
    const double rb = real_params[4];    
    const double rc = real_params[5];    
    for(int i = 0, i6 = 0; i < nmols; ++i, i6 += 6)
    {        
        const double s = 25.;
        const Vec3D un = points[i];
        const Vec3D u0 = center+Vec3D(points[i].x*ra, points[i].y*rb, points[i].z*rc);

        if(da.norm2() < 1.E-6) // No well-defined norm.
        {
            eCoords[i6] = 0.;
            eCoords[i6+1] = 0.;
            eCoords[i6+2] = 0.;
            eCoords[i6+3] = u0.x;
            eCoords[i6+4] = u0.y;
            eCoords[i6+5] = u0.z;            
        }
        else // Rotate.
        {
            const Vec3D an = da.normI();
            CalcExtCoord(a0, an, u0, un, eCoords+i6);
        }
    }
}

void StructureGenerator::DoLayer(int nmols, bool do_random, const Vec3D& a0, const Vec3D& da, const double* real_params, double* eCoords) const
{
    // Define the plane and norm.
    const Vec3D un(real_params[0], real_params[1], real_params[2]);
    const Vec3D v0(real_params[3], real_params[4], real_params[5]);
    const Vec3D v1(real_params[6], real_params[7], real_params[8]);
    const Vec3D v2(real_params[9], real_params[10], real_params[11]);
    const Vec3D v10 = v1-v0;
    const Vec3D v20 = v2-v0;
    vector<Vec3D> vs;
    const double p = v10.norm()/v20.norm();
    int n1 = sqrt(nmols*p);
    int n2 = sqrt(nmols/p);
    if(n1*n2 < nmols) { ++n1; ++n2; }
    for(int i = 0; i < n1; ++i)
    {
        for(int j = 0; j < n2; ++j)
        {
            vs.push_back(v10/n1*i+v20/n2*j);
        }
    }

    for(int i = 0, i6 = 0; i < nmols; ++i, i6 += 6)
    {        
        const Vec3D u0 = (do_random) ? v10*rand01()+v20*rand01()+v0 : vs[i]+v0;
        if(da.norm2() < 1.E-6) // No well-defined norm.
        {
            eCoords[i6] = 0.;
            eCoords[i6+1] = 0.;
            eCoords[i6+2] = 0.;
            eCoords[i6+3] = u0.x;
            eCoords[i6+4] = u0.y;
            eCoords[i6+5] = u0.z;            
        }
        else // Rotate.
        {
            const Vec3D an = da.normI();            
            CalcExtCoord(a0, an, u0, un, eCoords+i6);
        }
    }
}
/// + voild docircle ....

void StructureGenerator::CalcConvexHull(const vector<double>& points, vector<tuple<int, int, int, Vec3D> >& faces) const
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Ugly codes.
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    const int num_points = points.size()/3;
    ch_vertex* vertices = (ch_vertex*)malloc(num_points*sizeof(ch_vertex));
    for(int i = 0; i < num_points; ++i)
    {
        vertices[i].x = points[i*3];
        vertices[i].y = points[i*3+1];
        vertices[i].z = points[i*3+2];
    }    
    int* faceIndices = NULL;
    int nFaces;
    convhull_3d_build(vertices, num_points, &faceIndices, &nFaces);

    faces.clear();
    for(int i = 0; i < nFaces; ++i)
    {
        const int idx0 = faceIndices[i*3];
        const int idx1 = faceIndices[i*3+1];
        const int idx2 = faceIndices[i*3+2];
        ch_vec3 v0, v1, v2, normal;
        v0 = vertices[idx0];
        v1 = vertices[idx1];
        v2 = vertices[idx2];
        v1.x -= v0.x;
        v1.y -= v0.y;
        v1.z -= v0.z;
        v2.x -= v0.x;
        v2.y -= v0.y;
        v2.z -= v0.z;
        normal.x = v1.y*v2.z-v1.z*v2.y;
        normal.y = v1.z*v2.x-v1.x*v2.z;
        normal.z = v1.x*v2.y-v1.y*v2.x;
        /* normalise to unit length */
        double scale = 1.0/(sqrt(pow(normal.x, 2.0)+pow(normal.y, 2.0)+pow(normal.z, 2.0))+2.23e-9);
        normal.x *= scale;
        normal.y *= scale;
        normal.z *= scale;
        faces.push_back(make_tuple(idx0, idx1, idx2, Vec3D(normal.x, normal.y, normal.z)));
    }
    free(vertices);
    free(faceIndices);
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////    
}

struct OptArgs {
    Vec3D a0;
    Vec3D a1;
    Vec3D u0;
    Vec3D u1;
};

void StructureGenerator::CalcExtCoord(const Vec3D& a0, const Vec3D& an, const Vec3D& u0, const Vec3D& un, double* eCoords) const
{
    // Initialization.
    const int dim = 6;
    for(int i = 0; i < dim; ++i) { eCoords[i] = rand01(); }
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    // Search solutions.
    double delta;
    OptArgs args = { a0, a0+an.normI(), u0, u0+un.normI() };
    const int ret = lbfgs(dim, eCoords, &delta, evaluate_CalcExtCoord, progress_CalcExtCoord, (void*)(&args), &param);
}


int StructureGenerator::progress_CalcExtCoord(void *instance, const lbfgsfloatval_t *x, 
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

lbfgsfloatval_t StructureGenerator::evaluate_CalcExtCoord(void *instance, 
        const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
        const lbfgsfloatval_t step)
{
    OptArgs* oa = (OptArgs*)instance;
    const int num_atoms = 2;
    vector<double> iCoords(num_atoms*3, 0.);
    iCoords[0] = oa->a0.x; iCoords[1] = oa->a0.y; iCoords[2] = oa->a0.z;
    iCoords[3] = oa->a1.x; iCoords[4] = oa->a1.y; iCoords[5] = oa->a1.z;
    vector<double> eCoords(num_atoms*3, 0.);
    eCoords[0] = oa->u0.x; eCoords[1] = oa->u0.y; eCoords[2] = oa->u0.z;
    eCoords[3] = oa->u1.x; eCoords[4] = oa->u1.y; eCoords[5] = oa->u1.z;

    vector<double> coords(num_atoms*3, 0.);
    vector<double> derivatives(num_atoms*18, 0.);
    
    // Calculate ua0 and ua1.
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
    const double RX = x[3];
    const double RY = x[4];
    const double RZ = x[5];
    for(int i = 0, i3 = 0, i18 = 0; i < num_atoms; ++i, i3 += 3, i18 += 18)
    {
        const double x0 = iCoords[i3];
        const double y0 = iCoords[i3+1];
        const double z0 = iCoords[i3+2];
        const double t1 = RM[3]*x0+RM[4]*y0+RM[5]*z0;
        const double t2 = RM[0]*x0+RM[1]*y0+RM[2]*z0;
        // Coordinates
        coords[i3]   = t2+RX;
        coords[i3+1] = t1+RY;
        coords[i3+2] = RM[6]*x0+RM[7]*y0+RM[8]*z0+RZ;
        // Derivatives
            derivatives[i18]  = -RM[1]*x0+RM[0]*y0;
        derivatives[i18+1]  = RM[9]*x0+RM[10]*y0+RM[11]*z0;
        derivatives[i18+2]  = t1;
        derivatives[i18+3]  = 1.;
        derivatives[i18+4]  = 0.;
        derivatives[i18+5]  = 0.;
            derivatives[i18+6]  = -RM[4]*x0+RM[3]*y0;
        derivatives[i18+7]  = RM[12]*x0+RM[13]*y0+RM[14]*z0;
            derivatives[i18+8]  = -t2;
        derivatives[i18+9]  = 0.;
        derivatives[i18+10]  = 1.;
        derivatives[i18+11]  = 0.;        
            derivatives[i18+12]  = -RM[7]*x0+RM[6]*y0;
            derivatives[i18+13]  = RM[15]*x0+RM[16]*y0+RM[17]*z0;
            derivatives[i18+14] = 0.;
        derivatives[i18+15]  = 0.;
        derivatives[i18+16]  = 0.;
        derivatives[i18+17]  = 1.;
    }

    // Calculate delta and its derivatives.
    double delta = 0.;
    for(int i = 0; i < n; ++i) { g[i] = 0.; }
    for(int i = 0, i3 = 0, i18 = 0; i < num_atoms; ++i, i3 += 3, i18 += 18)
    {
        const double dx = coords[i3]-eCoords[i3];
        const double dy = coords[i3+1]-eCoords[i3+1];
        const double dz = coords[i3+2]-eCoords[i3+2];
        delta += dx*dx+dy*dy+dz*dz;
        for(int j = 0; j < 6; ++j) { g[j] += dx*derivatives[i18+j]; }
        for(int j = 0; j < 6; ++j) { g[j] += dy*derivatives[i18+6+j]; }
        for(int j = 0; j < 6; ++j) { g[j] += dz*derivatives[i18+12+j]; }
    }
    for(int i = 0; i < n; ++i) { g[i] *= 2; }
    return delta;
}

void StructureGenerator::GenerateUniformSpherePoints(int n, vector<Vec3D>& points) const
{
     // Initialization.
    const int dim = 2*n;
    vector<double> x(dim, 0.);
    for(int i = 0; i < dim; ++i) { x[i] = rand01()*Pi2; }
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    // Search solutions.
    double delta;    
    const int ret = lbfgs(dim, x.data(), &delta, evaluate_GenerateUniformSpherePoints, progress_GenerateUniformSpherePoints, NULL, &param);

    // To unit shpere.
    points.clear();
    for(int i = 0; i < dim; i += 2)
    {
        const double cthetai = cos(x[i]);
        const double sthetai = sin(x[i]);
        const double cphii = cos(x[i+1]);
        const double sphii = sin(x[i+1]);
        const double xi = cthetai*cphii;
        const double yi = cthetai*sphii;
        const double zi = sthetai;
        points.push_back(Vec3D(xi, yi, zi));
    }
}

int StructureGenerator::progress_GenerateUniformSpherePoints(void *instance, const lbfgsfloatval_t *x, 
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

lbfgsfloatval_t StructureGenerator::evaluate_GenerateUniformSpherePoints(void *instance, 
        const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
        const lbfgsfloatval_t step)
{
    const double nShift = 1.E-8;
    double delta = 0.;
    for(int i = 0; i < n; ++i) { g[i] = 0.; }
    for(int i = 0; i < n; i += 2)
    {        
        const double cthetai = cos(x[i]);
        const double sthetai = sin(x[i]);
        const double cphii = cos(x[i+1]);
        const double sphii = sin(x[i+1]);
        const double xi = cthetai*cphii;
        const double yi = cthetai*sphii;
        const double zi = sthetai;
        for(int j = 0; j < i; j += 2)
        {
            const double cthetaj = cos(x[j]);
            const double sthetaj = sin(x[j]);
            const double cphij = cos(x[j+1]);
            const double sphij = sin(x[j+1]);
            const double xj = cthetaj*cphij;
            const double yj = cthetaj*sphij;
            const double zj = sthetaj;

            const double xij = xi-xj;
            const double yij = yi-yj;
            const double zij = zi-zj;

            const double r2 = 1./(xij*xij+yij*yij+zij*zij+nShift);
            delta += r2;

            const double r4 = -2*r2*r2;
            const double dAdxi = r4*xij; const double dAdxj = -dAdxi;
            const double dAdyi = r4*yij; const double dAdyj = -dAdyi;
            const double dAdzi = r4*zij; const double dAdzj = -dAdzi;

            g[i]   += dAdxi*(-sthetai*cphii)+dAdyi*(-sthetai*sphii)+dAdzi*(cthetai);
            g[i+1] += dAdxi*(-cthetai*sphii)+dAdyi*(cthetai*cphii);
            g[j]   += dAdxj*(-sthetaj*cphij)+dAdyj*(-sthetaj*sphij)+dAdzj*(cthetaj);
            g[j+1] += dAdxj*(-cthetaj*sphij)+dAdyj*(cthetaj*cphij);
        }
    }
    return delta;
}

void StructureGenerator::CalcMeanPlane(const vector<Vec3D>& points, vector<Vec3D>& proj_points) const
{
    // Initialization.
    const int dim = 3;
    vector<double> x(dim, 0.); // theta, phi, d
    for(int i = 0; i < dim; ++i) { x[i] = rand01(); }
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    // Search solutions.
    double dist;    
    const int ret = lbfgs(dim, x.data(), &dist, evaluate_CalcMeanPlane, progress_CalcMeanPlane, (void*)(&points), &param);
}

int StructureGenerator::progress_CalcMeanPlane(void *instance, const lbfgsfloatval_t *x, 
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

lbfgsfloatval_t StructureGenerator::evaluate_CalcMeanPlane(void *instance, 
        const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, 
        const lbfgsfloatval_t step)
{
    const vector<Vec3D>& points = *((const vector<Vec3D>*)instance);
    const int num_points = points.size();

    double dist = 0.;
    for(int i = 0; i < n; ++i) { g[i] = 0.; }
    for(int i = 0; i < num_points; ++i)
    {
        const double ct = cos(x[0]);
        const double st = sin(x[0]);
        const double cp = cos(x[1]);
        const double sp = sin(x[1]);
        const double d = x[2];
        const double di = ct*cp*points[i].x+ct*sp*points[i].y+st*points[i].z+d;
        dist += di*di;
        g[0] += di*2*(-st*cp*points[i].x-st*sp*points[i].y+ct*points[i].z);
        g[1] += di*2*(-ct*sp*points[i].x+ct*cp*points[i].y);
        g[2] += di*2;
    }
    dist /= num_points;
    g[0] /= num_points;
    g[1] /= num_points;
    g[2] /= num_points;
    return dist;
}
