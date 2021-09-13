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

#ifndef    __ABCLUSTERINPUT_H__
#define    __ABCLUSTERINPUT_H__

#include <vector>
#include <string>

namespace abcluster {

using namespace std;

class ABClusterInput {
public:
    ABClusterInput(const string& tinpfn);
    ~ABClusterInput(void); 

    string           inpfn;           // Input file name     
    string           resultfn;        // result file name
    string           clusterfn;       // Cluster file name
    int              maxncals;        // The maximal number of calculations
    vector<string>   structuretypes;  // Structure types    
    vector<string>   optcmds;         // External program excuting commands    
    bool             docoarseopt;     // Flag that if we should do the coarse optimization.
    double           energy_ub;       // Energy upper bound
    double           energy_lb;       // Energy lower bound

    static const string NoNodes;    
private:    

};

}

#endif //  __ABCLUSTERINPUT_H__
