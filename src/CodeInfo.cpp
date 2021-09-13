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

#include "CodeInfo.h"
#include <vector>
#include <cstdio>

using namespace abcluster;
using namespace std;

const string CodeInfo::name = "NWPEsSe";
const string CodeInfo::flag = CODEFLAG;
const string CodeInfo::majorversion = CODEMAJORVER;
const string CodeInfo::minorversion = CODEMINORVER;
Author tauthors[] = {
    Author("Jun Zhang, Difan Zhang, Vanda Glezakou ", "Pacific Northwest National Laboratory", ""),

};
const vector<Author> CodeInfo::authors(tauthors, tauthors+1); 
const CompilingInfo CodeInfo::compilerinfo(
        CODECOMPILER,
        CODECOMPFLAG,
        CODEUSER,
	string(__DATE__)+" "+__TIME__,
	CODEMACHINE
	);

CodeInfo::CodeInfo(void)
{
    /* Nothing to do */ 
}

CodeInfo::~CodeInfo(void) 
{ 
    /* Nothing to do */ 
}

const string& CodeInfo::GetName(void)
{
    return name;
}

const string& CodeInfo::GetFlag(void)
{
    return flag;
}

const string& CodeInfo::GetMajorVersion(void)
{
    return majorversion;
}


const string& CodeInfo::GetMinorVersion(void)
{
    return minorversion;
}

const vector<Author>& CodeInfo::GetAuthors(void)
{
    return authors;
}

const CompilingInfo& CodeInfo::GetCompilerinfo(void)
{
    return compilerinfo;
}

void CodeInfo::PrintCodeInfo(void)
{
    PrintCodeLogo();
    printf("\n");
    PrintCodeVersion();
    printf("\n");
    PrintCodeAuthors();
    printf("\n");
    PrintPrologue();
}

void CodeInfo::PrintCodeLogo(void)
{
    // Print logo.
printf("  _   ___          _______  ______      _____       \n");
printf(" | \\ | \\ \\        / /  __ \\|  ____|    / ____|      \n");
printf(" |  \\| |\\ \\  /\\  / /| |__) | |__   ___| (___   ___  \n");
printf(" | . ` | \\ \\/  \\/ / |  ___/|  __| / __|\\___ \\ / _ \\ \n");
printf(" | |\\  |  \\  /\\  /  | |    | |____\\__ \\____) |  __/ \n");
printf(" |_| \\_|   \\/  \\/   |_|    |______|___/_____/ \\___| \n");
printf("                                                    \n");
  
}

void CodeInfo::PrintCodeVersion(void)
{
    // Print version.
    printf("%s %s.%s\n", name.c_str(), majorversion.c_str(), minorversion.c_str());
    printf("Specific flag:  %s\n", flag.c_str());

    // Print compiling informaiton.
    printf("Platform:       %s\n", compilerinfo.platform.c_str());
    printf("Compiled by:    %s@%s\n",compilerinfo.user.c_str(), 
	    compilerinfo.machine.c_str());
    printf("Compiling date: %s\n", compilerinfo.time.c_str());
    printf("C++ compiler:   %s\n", compilerinfo.compiler.c_str());
    printf("C++ options:    %s\n", compilerinfo.compflag.c_str());
}

void CodeInfo::PrintCodeAuthors(void)
{
    printf("Contributors:\n");
    for(vector<Author>::const_iterator iter = authors.begin();
            iter != authors.end(); ++iter)
    {
        printf("%s(%s): %s\n", iter->name.c_str(), iter->college.c_str(), iter->contribution.c_str());
    }
}

void CodeInfo::PrintPrologue(void)
{
    printf("For questions on the theory and ");
    printf("implementation of %s, please contact:\n", name.c_str());
    printf("    Vanda Glezakou, vanda.glezakou@pnnl.gov\n");
    printf("    Difan Zhang, difan.zhang@pnnl.gov\n");
    printf("\n");
    printf("NWPEsSe references:\n");
    printf("[1] Zhang, J.; Glezakou, V.-A.; Rousseau, R.; Nguyen, M.-T.; NWPEsSe: An Adaptive-Learning\n");
    printf("    Global Optimization Algorithm for Nanosized Cluster Systems. J. Chem. Theory Comput. 2020,\n");
    printf("    16, 3947-3958\n");
    printf("\n");
    printf("Other references:\n");
    printf("[1] Zhang, J.; Dolg, M.; ABCluster: The Artificial Bee Colony Algorithm for\n");
    printf("    Cluster Global Optimization. Phys. Chem. Chem. Phys. 2015, 17, 24173-24181.\n");
    printf("[2] Zhang, J.; Dolg, M.; Global Optimization of Clusters of Rigid Molecules by the\n");
    printf("    Artificial Bee Colony Algorithm. Phys. Chem. Chem. Phys. 2016, 18, 3003-3010.\n");

    printf("\n");
    //printf("For more efficient global optimization of atomic and molecular clusters with force fields,\n");
    //printf("please use ABCluster from \n");
    //printf("\n");
    //printf("http://www.zhjun-sci.com/software-abcluster-download.php\n");
    //printf("\n");
}
