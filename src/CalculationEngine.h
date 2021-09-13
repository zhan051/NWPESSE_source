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

#ifndef    __CALCULATIONENGINE_H__
#define    __CALCULATIONENGINE_H__

#include "CEResult.h"
#include "RMCluster.h"
#include <string>
#include <vector>

namespace abcluster {

using namespace std;

// For external program 
class CalculationEngine {
public:
    // Constructor.
    CalculationEngine(void);
    // Destructor.    
    ~CalculationEngine(void);

    // Perform a local optimization.
    // Return: double Energy, optimized X
    double LocalOptimization(RMCluster& rmc, string& state, CEResult& cer, const string& label, const vector<string>& optcmds, const string& node);

private:
    void RunExternalProgram(const vector<string>& optcmds, const string& inpfn, const string& outfn, const string& label, const string& node) const;
    void WriteInputXYZFile(const RMCluster& rmc, const string& fn) const;
    void ReadOutputXYZFile(RMCluster& rmc, string& state, CEResult& cer, const string& fn) const;
    void DeleteTmpFiles(const string& inpfn, const string& outfn) const;
};

}

#endif //  __CALCULATIONENGINE_H__
