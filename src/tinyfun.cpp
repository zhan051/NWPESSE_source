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

#include "tinyfun.h"
#if defined(CODEWINDOWS)
#include <windows.h>
#elif defined(CODELINUX) || defined(CODEMACOSX)
#include <sys/stat.h>
#else
#error You should assign the operating system: CODEWINDOWS or CODELINUX or CODEMACOSX
#endif
#include <string>
#include <cstdlib>
#include <ctime>
#include <boost/random.hpp>

using namespace boost;
using namespace std;

typedef lagged_fibonacci19937 abclusterrand;
//typedef boost::mt19937 abclusterrand;
abclusterrand rng(time(NULL));
uniform_01<abclusterrand&> u01(rng);  

void abcluster::srand(void)
{
    std::srand(time(NULL));   
}

double abcluster::rand01(void)
{
    return u01();
}

int abcluster::randnotn(int min, int max, int n)
{
    int d = max-min;
    int k;
    for(int i = 0; i < 100; ++i)
    {
        k = rand01()*d+min;
	if(k != n) return k;
    }
    return k;
}

double abcluster::continuestep01(double x, double shift, int nit)
{
    if(shift == 0.0000) return 0;
    double a = -1.0000/shift;
    double b = 1.0000-a;
    double y = x;
    for(int i = 0; i < nit; ++i)
    {
        y = y*y*(a*y+b);
    }
    return y;
}

string abcluster::replaceall(const string& str, const string& oldstr, const string& newstr)
{
     string str0 = str;
     while(true)   
     {   
	 size_t pos = str0.find(oldstr);				       
	 if(pos != string::npos)   
	 {
	     str0.replace(pos, oldstr.length(),newstr);
	 }
	 else
	 {
	     break;   
	 }
     }   
     return str0;  
}

bool abcluster::makedir(const string& str)
{
#if defined(CODEWINDOWS)
    return (CreateDirectory(str.c_str(), NULL) == 1) ? true : false;
#elif defined(CODELINUX) || defined(CODEMACOSX)
    return (mkdir(str.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0) ? true : false;     
#else
#error You should assign the operating system: CODEWINDOWS or CODELINUX or CODEMACOSX
    return false;
#endif
}

string abcluster::triminput(const string& str) 
{
    // Delete the comment.
    string tstr = str.substr(0, str.find('#'));
    // Delete the spaces.
    tstr = tstr.substr(0, tstr.find_last_not_of(" \n\r\t")+1);
    if(!tstr.empty())
    {   
	tstr = tstr.substr(tstr.find_first_not_of(" \n\r\t"));
    }   
    return tstr;
}




