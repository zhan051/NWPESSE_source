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

#ifndef    __RUNSTATUS_H__
#define    __RUNSTATUS_H__

#include <string>
#include <cstdarg>
#include <ctime>

namespace abcluster {

using namespace std;

class RunStatus {
public:
    RunStatus(void);
    ~RunStatus(void);

    // Runtime status control.
    static void Start(int argc, char** argv);
    static unsigned int NormalTermination(void);
    static unsigned int ErrorTermination(const char* fmt, ...);
    static void Warning(const char* fmt, ...);
    static unsigned int GetPid(void);
    static unsigned int GetTid(void);
    static string GetUser(void);
    static string GetHostname(void);
    static string GetTime(time_t& lt);
    static string GetRunningTime(const time_t& latet, const time_t& earlyt, unsigned int& duration);
    static string Get_OMP_NUM_THREADS(void);
    static string GetCurrentPath(void);

    static const unsigned int Succeed;
    static const unsigned int ErrorExit;

private:

    static time_t starttime;
    static time_t terminationtime;
    static unsigned int pid;
    static string user;
    static string hostname;

    static int GetNCPUCores(void);
    static void GetUserInfo(void);
    static void ProcessSignal(void);
};

void BlockTermSignal(int signal);

}

#endif //  __RUNSTATUS_H__
