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

#include "RunStatus.h"
#include "CodeInfo.h"
#if defined(CODEWINDOWS)
#include <windows.h>
#include <direct.h>
#elif defined(CODELINUX)
#include <sys/sysinfo.h>
#include <pwd.h>
#include <unistd.h>
#include <pthread.h>
#elif defined(CODEMACOSX)
#include <sys/types.h>
#include <sys/sysctl.h>
#include <pwd.h>
#include <unistd.h>
#include <pthread.h>
#else
#error You should assign the operating system: CODEWINDOWS or CODELINUX or CODEMACOSX
#endif
#include <string>
#include <csignal>
#include <cstdarg>
#include <ctime>
#include <cstdlib>
#include <cstdio>

using namespace abcluster;
using namespace std;

const unsigned int RunStatus::Succeed = 0x00000000;
const unsigned int RunStatus::ErrorExit = 0x00000001;

time_t RunStatus::starttime;
time_t RunStatus::terminationtime;
unsigned int RunStatus::pid;
string RunStatus::user;
string RunStatus::hostname;

RunStatus::RunStatus(void) 
{ 
    /* Nothing to do */ 
}

RunStatus::~RunStatus(void)
{
    /* Nothing to do */ 
}

void RunStatus::Start(int argc, char** argv)
{
    setvbuf(stdout, NULL, _IONBF, 0); // No buffer!
    ProcessSignal(); // Blocking signals
    // Begin!
    CodeInfo::PrintCodeInfo();
    printf("\n");
    printf("Program starts at %s.\n", GetTime(starttime).c_str());
    printf("Command line: \" ");
    for(int i = 0; i < argc; ++i)
    {
        printf("%s ", argv[i]);
    }
    printf("\"\n");
    GetUserInfo();
    printf("Process ID:   %u\n", pid);
    printf("Run by:       %s@%s\n", user.c_str(), hostname.c_str());
    printf("Available CPU cores:    %d\n", GetNCPUCores());
    printf("Read \"OMP_NUM_THREADS\": %s\n", Get_OMP_NUM_THREADS().c_str());
    printf("\n");
}

unsigned int RunStatus::NormalTermination(void)
{
    printf("\n");
    printf("===== (^ _ ^) =====\n");
    printf("Normal termination at %s.\n", GetTime(terminationtime).c_str());
    
    unsigned int duration;
    string runtimestr;
    runtimestr = GetRunningTime(terminationtime, starttime, duration);
    printf("Total run time: %u seconds = %s.\n", duration, runtimestr.c_str());
    
    return Succeed;
}

unsigned int RunStatus::ErrorTermination(const char* fmt, ...)
{
    va_list args;
    char buffer[2048];
    va_start(args, fmt);
    vsprintf(buffer, fmt, args);
    va_end(args);
    printf("\n");
    printf("===== (@ o @) =====\n");
    printf("Error occurs: %s\n", buffer);
    printf("Error termination at %s.\n", GetTime(terminationtime).c_str());

    unsigned int duration;
    string runtimestr;
    runtimestr = GetRunningTime(terminationtime, starttime, duration);
    printf("Total run time: %u seconds = %s.\n", duration, runtimestr.c_str());

    exit(ErrorExit);
    return ErrorExit;
}

void RunStatus::Warning(const char* fmt, ...)
{
    va_list args;
    char buffer[2048];
    va_start(args, fmt);
    vsprintf(buffer, fmt, args);
    va_end(args);
    printf(" ! ! ! WARNING: %s\n", buffer);
}

unsigned int RunStatus::GetPid(void)
{
    return pid;
}

unsigned int RunStatus::GetTid(void)
{
#if defined(CODEWINDOWS)
    return GetCurrentThreadId();
#elif defined(CODELINUX)
    return pthread_self();
#elif defined(CODEMACOSX)
    return pthread_self()->__sig;    
#else
#error You should assign the operating system: CODEWINDOWS or CODELINUX or CODEMACOSX
    return 1;
#endif    
}

string RunStatus::GetUser(void)
{
    return user;
}

string RunStatus::GetHostname(void)
{
    return hostname;
}

string RunStatus::Get_OMP_NUM_THREADS(void)
{
    char* p = getenv("OMP_NUM_THREADS");
    return (p == NULL)?string("Unset"):string(p);
}

string RunStatus::GetCurrentPath(void)
{
    string path;
    char* buffer = getcwd(NULL, 0);
    path.assign(buffer);
    free(buffer);
    return path;
}

string RunStatus::GetTime(time_t& lt)
{
    time(&lt);
    string timestr(asctime(localtime(&lt)));
    return timestr.substr(0, timestr.find('\n'));
}

string RunStatus::GetRunningTime(const time_t& latet, const time_t& earlyt,
	unsigned int& duration)
{
    duration = (unsigned int)difftime(latet, earlyt);
    
    const unsigned int s2d = 86400;
    const unsigned int s2h = 3600;
    const unsigned int s2m = 60;
    unsigned int day, hour, min, second;
    
    day = duration/s2d;
    hour = duration%s2d/s2h;
    min = duration%s2d%s2h/s2m;
    second = duration%s2d%s2h%s2m;

    char buffer[2048];
    sprintf(buffer, "%u days %u hours %u minutes %u seconds", day, hour,
	    min, second);
    return string(buffer);
}

int RunStatus::GetNCPUCores(void)
{
#if defined(CODEWINDOWS)
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors;
#elif defined(CODELINUX)
    return get_nprocs();
#elif defined(CODEMACOSX)
    const size_t buffersize = 32;
    size_t strsize = 32;
    char buffer[buffersize];
    sysctlbyname("hw.physicalcpu", &buffer, &strsize, NULL, NULL);
    return atoi(buffer);    
#else
#error You should assign the operating system: CODEWINDOWS or CODELINUX or CODEMACOSX
    return 0;
#endif    
}

void RunStatus::GetUserInfo(void)
{
#if defined(CODEWINDOWS)
    pid = GetCurrentProcessId();
    const DWORD buffersize = 2048;
    char buffer[buffersize];
    DWORD strsize;
    strsize = buffersize;
    GetUserName(buffer, &strsize);
    user.assign(buffer);
    strsize = buffersize;    
    GetComputerName(buffer, &strsize);
    hostname.assign(buffer);
#elif defined(CODELINUX) || defined(CODEMACOSX)
    pid = (unsigned int)getpid();
    user.assign(getpwuid(geteuid())->pw_name);
    const unsigned int buffersize = 2048;
    char buffer[buffersize];
    gethostname(buffer, buffersize);
    hostname.assign(buffer);
#else
    #error You should assign the operating system: CODEWINDOWS or CODELINUX or CODEMACOSX
    pid = 1;
    user.assign("?");
    hostname.assign("?");
#endif
}

void RunStatus::ProcessSignal(void)
{
    signal(SIGINT, abcluster::BlockTermSignal);
    signal(SIGILL, abcluster::BlockTermSignal);
    signal(SIGFPE, abcluster::BlockTermSignal);
    signal(SIGSEGV, abcluster::BlockTermSignal);
    signal(SIGTERM, abcluster::BlockTermSignal);
    signal(SIGABRT, abcluster::BlockTermSignal);
}

void abcluster::BlockTermSignal(int signal)
{
    string siginfo;
    switch (signal)
    {
	case SIGINT:
	{
	    siginfo.assign("SIGINT: Interrupt");
	    break;
	}
	case SIGILL:
	{
	    siginfo.assign("SIGILL: Illegal instruction");
	    break;
	}
	case SIGFPE:
	{
	    siginfo.assign("SIGFPE: Floating point exception");
	    break;
	}
	case SIGSEGV:
	{
	    siginfo.assign("SIGSEGV: Segment violation");
	    break;
	}
	case SIGTERM:
	{
	    siginfo.assign("SIGTERM: Software termination from kill");
	    break;
	}
	case SIGABRT:
	{
	    siginfo.assign("SIGABRT: Termination by abort call");
	    break;
	}
	default:
	{
	    siginfo.assign("Unknown signal");
	    break;
	}
    }
    RunStatus::ErrorTermination("Recieve signal: %s (Value: %d).", siginfo.c_str(), signal);
}

