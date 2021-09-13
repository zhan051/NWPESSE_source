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

#ifndef    __CODEINFO_H__
#define    __CODEINFO_H__

#include <string>
#include <vector>

namespace abcluster {

using namespace std;

class Author {
public:
    Author(const string& tname, const string& tcollege, const string& tcontribution) :
	name(tname), college(tcollege), contribution(tcontribution)
    {
        /* Nohting to do. */
    }
    string name;
    string college;
    string contribution;
};

class CompilingInfo {
public:
    CompilingInfo(const string& tcompiler, const string& tcompflag,
            const string& tuser, const string& ttime, const string& tmachine) :
	compiler(tcompiler), compflag(tcompflag),
        user(tuser), time(ttime), machine(tmachine)
    {
#if defined(CODEWINDOWS)
        platform.assign("WINDOWS");
#elif defined(CODELINUX)
        platform.assign("LINUX/UNIX");
#elif defined(CODEMACOSX)
        platform.assign("Mac OS X");
#else
        #error You should assign the operating system: CODEWINDOWS or CODELINUX or CODEMACOSX
        platform.assign("Unknown");
#endif
    }
    string compiler;
    string compflag;
    string user;
    string time;
    string machine;
    string platform;
};

class CodeInfo {
public:
    CodeInfo(void);
    ~CodeInfo(void);

    static const string& GetName(void);
    static const string& GetFlag(void);
    static const string& GetMajorVersion(void);
    static const string& GetMinorVersion(void);
    static const vector<Author>& GetAuthors(void);
    static const CompilingInfo& GetCompilerinfo(void);
    
    static void PrintCodeInfo(void);

private:
    static const string name;
    static const string flag;
    static const string majorversion;
    static const string minorversion;
    static const vector<Author> authors;
    static const CompilingInfo compilerinfo;    

    static void PrintCodeLogo(void);
    static void PrintCodeVersion(void);
    static void PrintCodeAuthors(void);    
    static void PrintPrologue(void);    
};

}

#endif //  __CODEINFO_H__
