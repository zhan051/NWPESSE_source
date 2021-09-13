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

#include "Parallelization.h"
#include "RunStatus.h"
#include <boost/spirit/include/classic.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <pthread.h>

using namespace abcluster;
using namespace boost::spirit::classic;
using namespace std;

Parallelization::Parallelization(void) 
{
    /* Nothing to do */
}

Parallelization::~Parallelization(void) 
{
    /* Nothing to do */
}

void Parallelization::Initialize(const string& tfn)
{
    fn = tfn;
    // Open nodes file
    ifstream fd(fn.c_str(), ios::in);
    if(!fd.is_open())
    {
	RunStatus::ErrorTermination("Cannot open input file [ %s ].", fn.c_str());
    }
    string filecontent((istreambuf_iterator<char>(fd)), istreambuf_iterator<char>());
    fd.close();

    string te;
    int tn;
    const rule<> RULE_str = (+graph_p)[assign(te)];
    const rule<> RULE_int = (int_p)[assign(tn)];
    parse_info<> parser;

    string tstr;
    istringstream iss(filecontent);
    getline(iss, tstr, '\n');
    parser = parse(tstr.c_str(), RULE_int);
    if(!parser.full || tn <= 0)
    {
	RunStatus::ErrorTermination("Incorrect node file [ %s ] format near \" %s \": A positive number is required.", fn.c_str(), tstr.c_str());
    }
    const int num_nodes = tn;
    node_names.clear(); node_names.reserve(num_nodes);
    is_available.clear(); is_available.reserve(num_nodes);
    for(unsigned int i = 0; i < num_nodes; ++i)
    {
	getline(iss, tstr, '\n');
	parser = parse(tstr.c_str(), RULE_str);
        if(!parser.full)
        { 
	    RunStatus::ErrorTermination("Incorrect node file [ %s ] format near: \" %s \". One host name or IP is needed.", fn.c_str(), tstr.c_str());
        }
        node_names.push_back(te);
        is_available.push_back(true);
    }

    thread_pool.clear();
}

void Parallelization::PrintNames(void) const
{
    for(int i = 0; i < node_names.size(); ++i)
    {
        printf("%3d: %s\n", i, node_names[i].c_str());
    }
}

string Parallelization::GetName(int i) const
{
    return node_names[i];
}

int Parallelization::GetNumNodes(void) const
{
    return node_names.size();
}

void Parallelization::AddThread(pthread_t pt) 
{ 
    thread_pool.push_back(pt); 
}

bool Parallelization::IsAvailable(int i) const
{
    return is_available[i];
}

void Parallelization::SetAvailable(int i)
{
    is_available[i] = true;
}

void Parallelization::SetUnavailable(int i)
{
    is_available[i] = false;
}

void Parallelization::Barrier(void) const
{ 
    for(int i = 0; i < thread_pool.size(); ++i)
    {
        pthread_join(thread_pool[i], NULL);
    }
}
