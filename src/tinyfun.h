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

#ifndef    __TINYFUN_H__
#define    __TINYFUN_H__

#include <string>

namespace abcluster {

using namespace std;

void srand(void);
// Generate a random real number in [0,1).
// Warning: Before using it, one has to set up a random number seed manually.
double rand01(void);

// Generate a random integer in [min,max) that is not n!
int randnotn(int min, int max, int n);

// Generating step function by nit iteration of 
//     y = ax^3+bx^2
// where the crosspoint with y = x is at shift...
double continuestep01(double x, double shift, int nit);

// Replace all oldstr by newstr occured in str.
string replaceall(const string& str, const string& oldstr, const string& newstr);

// Create a new directory.
bool makedir(const string& str);

// Trim by '#'
string triminput(const string& str);
}

#endif //  __TINYFUN_H__
