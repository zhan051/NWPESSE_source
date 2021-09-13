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

#ifndef    __PARALLELIZATION_H__
#define    __PARALLELIZATION_H__

#include <string>
#include <vector>
#include <pthread.h>

namespace abcluster {

using namespace std;

// For external program 
class Parallelization {
public:
    // Constructor.
    Parallelization(void);
    // Destructor.    
    ~Parallelization(void);
    
    void Initialize(const string& tfn);

    void PrintNames(void) const;
    string GetName(int i) const;
    int GetNumNodes(void) const;
    bool IsAvailable(int i) const;

    void SetAvailable(int i);
    void SetUnavailable(int i);
    void AddThread(pthread_t pt);    
    void Barrier(void) const;

private:
    string fn;
    vector<string>    node_names;
    vector<bool>      is_available;
    vector<pthread_t> thread_pool;;
};

}

#endif //  __PARALLELIZATION_H__
