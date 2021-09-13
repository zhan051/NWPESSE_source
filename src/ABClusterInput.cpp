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

#include "ABClusterInput.h"
#include "RunStatus.h"
#include "tinyfun.h"
#include <string>
#include <vector>
#include <boost/spirit/include/classic.hpp>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace abcluster;
using namespace boost::spirit::classic;
using namespace std;

const string ABClusterInput::NoNodes = "*";

ABClusterInput::ABClusterInput(const string& tinpfn)
{
    inpfn.assign(tinpfn);

    string te;
    int tn;
    const rule<> RULE_str = (+graph_p)[assign(te)];
    const rule<> RULE_int = (int_p)[assign(tn)];
    parse_info<> parser;

    // Open the input file.
    ifstream fd(inpfn.c_str(), ios::in);
    if(!fd.is_open())
    {
	RunStatus::ErrorTermination("Cannot open input file [ %s ].", inpfn.c_str());
    }
    string filecontent((istreambuf_iterator<char>(fd)), istreambuf_iterator<char>());
    fd.close();

    string tstr;    
    istringstream iss(filecontent);
  
    // Get result file name
    getline(iss, tstr, '\n');
    parser = parse((triminput(tstr)).c_str(), RULE_str);
    if(!parser.full) 
    {
	RunStatus::ErrorTermination("Error input in [ %s ] near \n    %s \nCannot read the result file name.", inpfn.c_str(), tstr.c_str());
    }
    resultfn.assign(te);

    // Get symbols
    getline(iss, tstr, '\n');
    clusterfn.assign(triminput(tstr));

    // Get maximal number of calculations
    getline(iss, tstr, '\n');
    parser = parse((triminput(tstr)).c_str(), RULE_int);
    if(!parser.full || tn <= 0) 
    {
	RunStatus::ErrorTermination("Error input in [ %s ] near \n    %s \nA positive number is required.", inpfn.c_str(), tstr.c_str());
    }
    maxncals = tn;

    // Get external program excuting commands
    const char* Separator = ">>>>";
    const rule<> RULE_separator = (*blank_p)>>str_p(Separator)>>(*blank_p);
    getline(iss, tstr, '\n');
    parser = parse((triminput(tstr)).c_str(), RULE_separator);
    if(!parser.full)
    {
        RunStatus::ErrorTermination("Structure types should start by \"%s\".", Separator);        
    }
    structuretypes.clear();
    while(1)
    {
        if(getline(iss, tstr, '\n').eof())
        {
            RunStatus::ErrorTermination("Commands should terminate by \"%s\".", Separator);
        }
        parser = parse((triminput(tstr)).c_str(), RULE_separator);
        if(!parser.full)
        {
            structuretypes.push_back(triminput(tstr));
        }
        else
        {
            break;
        }
    }
    optcmds.clear();
    while(1)
    {
    	if(getline(iss, tstr, '\n').eof())
    	{
            if(parse((triminput(tstr)).c_str(), RULE_separator).full) { break; }
    	    RunStatus::ErrorTermination("Commands should terminate by \"%s\".", Separator);
    	}
        parser = parse((triminput(tstr)).c_str(), RULE_separator); 
    	if(!parser.full)
    	{
    	    optcmds.push_back(triminput(tstr));
    	}
    	else
    	{
    	    break;
    	}
    }

    // Get optimization strategy.    
    docoarseopt = true;
    energy_ub = +1.E+10;
    energy_lb = -1.E+10;
    tstr = "";

    getline(iss, tstr, '\n');
    const rule<> RULE_optstrategy = (*blank_p)>>(+graph_p)[assign(te)]>>(*blank_p)>>real_p[assign(energy_lb)]>>(*blank_p)>>real_p[assign(energy_ub)]>>(*blank_p);
    if(parse((triminput(tstr)).c_str(), RULE_optstrategy).full)
    {
        if(te == "n" || te == "N") { docoarseopt = false; }
        else
        {
            if(te == "y" || te == "Y") { docoarseopt = true; }
            else { RunStatus::ErrorTermination("Cannot understand \"%s\". ", te.c_str()); }
        }
        if(energy_ub < energy_lb) { double td = energy_lb; energy_lb = energy_ub; energy_ub = td; }
    }
    else
    {
        if(tstr == "n" || tstr == "N") { docoarseopt = false; }
        else
        {
            if(tstr == "y" || tstr == "Y") { docoarseopt = true; }
            else
            { 
                if(!triminput(tstr).empty()) { RunStatus::ErrorTermination("Cannot understand \"%s\". ", tstr.c_str()); }
            }
        }
    }    
}

ABClusterInput::~ABClusterInput(void)
{
    // Nothing to do
}

