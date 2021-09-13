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

#include "EnergyOrder.h"
#include "Parallelization.h"
#include "CalculationEngine.h"
#include "StructureGenerator.h"
#include "CEResult.h"
#include "RMCluster.h"
#include "ABClusterInput.h"
#include "tinyfun.h"
#include "RunStatus.h"
#include <algorithm>
#include <utility>
#include <string>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <pthread.h>
#include <unistd.h>


using namespace abcluster;
using namespace std;

class TaskArgs {
public:
    TaskArgs(Parallelization* tpdispatcher, int tnode_idx,
            int ti, const string& tlmdir, StructureGenerator* tpgeompool,
            ABClusterInput* tpabcinp, EnergyOrder* tppool,
            CalculationEngine* tpeg, RMCluster* tprmc) : 
        pdispatcher(tpdispatcher), node_idx(tnode_idx),
        i(ti), lmdir(tlmdir), pgeompool(tpgeompool), 
        pabcinp(tpabcinp), ppool(tppool), peg(tpeg), prmc(tprmc)
    {}
    Parallelization*    pdispatcher;
    int                 node_idx;

    int                 i;          // i.
    string              lmdir;      // Directory
    StructureGenerator* pgeompool;  // Structure.
    ABClusterInput*     pabcinp;    // ABC input.
    EnergyOrder*        ppool;      // Energy pool.
    CalculationEngine*  peg;        // Energy eigine.
    RMCluster*          prmc;       // Seed cluster.
};

bool ChangeLine = true;

void DoTask(TaskArgs& args)
{
    Parallelization&    dispatcher = *(args.pdispatcher);    
    const int           node_idx   = args.node_idx;
    const int           i          = args.i;
    const string        lmdir      = args.lmdir;
    StructureGenerator& geompool   = *(args.pgeompool);
    ABClusterInput&     abcinp     = *(args.pabcinp);
    EnergyOrder&        pool       = *(args.ppool);
    CalculationEngine&  eg         = *(args.peg);
    RMCluster&          rmc        = *(args.prmc);
// -------------------------------------------------------------------------
    char buffer[1024];
    sprintf(buffer, "abcluster%s%s%u-%d", RunStatus::GetUser().c_str(), RunStatus::GetHostname().c_str(), RunStatus::GetPid(), i);
    geompool.OneCluster(abcinp.docoarseopt, rmc);  //////// first initial, update with bees
    
    string node;
    if(dispatcher.GetNumNodes() > 0)
    { 
        node = dispatcher.GetName(node_idx);
        if(ChangeLine) { printf("Doing "); ChangeLine = false; }
        printf("%d@%s, ", i, node.c_str());
    }
    else
    {
        printf("%5d ", i);
    }
    
    string CurrentDirectory(RunStatus::GetCurrentPath());
    string state;
    CEResult cer;    ///////////////// contins str and ene
    time_t t1; RunStatus::GetTime(t1);
    eg.LocalOptimization(rmc, state, cer, string(CurrentDirectory)+"/"+string(buffer), abcinp.optcmds, node); // <- Calculation
    time_t t2; RunStatus::GetTime(t2);
    unsigned int duration; RunStatus::GetRunningTime(t2, t1, duration);
    sprintf(buffer, "%s/%d", lmdir.c_str(), i);                
    cer.Save(string(buffer));   //////////// save str to disk

    if(state == "Succeed"){
        vector<ExtCoord> eCoords;
        cer.SolveExtCoords(rmc, eCoords); // Solve its external coordinates.
        geompool.AddClusterForTraining(cer.GetEnergy(), eCoords);  // Add to the training pool.
    }

    if(dispatcher.GetNumNodes() > 0)
    { 
        if(!ChangeLine) { printf("\n"); }
        printf("%5d %15.8f %10u   %s\n", i, cer.GetEnergy(), duration, state.c_str());
        ChangeLine = true;
    }
    else
    {
        printf("%15.8f %10u   %s\n", cer.GetEnergy(), duration, state.c_str());
    }
    if(cer.GetEnergy() >= abcinp.energy_lb && cer.GetEnergy() <= abcinp.energy_ub)
    {
        pool.AddEI(i, cer.GetEnergy());
    }
// -------------------------------------------------------------------------
}

void* DoTaskThread(void* pargs)
{
    TaskArgs* pt = (TaskArgs*)(pargs);
    DoTask(*pt);
    pt->pdispatcher->SetAvailable(pt->node_idx);
    delete pt;
    pthread_exit(NULL);
}


int main(int argc, char** argv)
{
    RunStatus::Start(argc, argv);
    
    abcluster::srand();
    printf("Random number testing (integers): %d...%d\n", rand(), rand());
    printf("Random number testing (fractals): %.8f...%.8f\n", rand01(), rand01());
    printf("\n");

    if(argc < 2)
    {
        return RunStatus::ErrorTermination("Need one input file name.");
    }
    
    ABClusterInput abcinp(argv[1]);   //  input files
    string lmdir = abcinp.resultfn+"-LM";
    if(!makedir(lmdir))
    {
	RunStatus::ErrorTermination("Cannot create the directory [ %s ] to store the local minima. If this directory already exists, please delete it.", lmdir.c_str());
    }        

    printf(" -- Input/Output Parameters --\n");    
    printf("Input file name:       %s\n", abcinp.inpfn.c_str());
    printf("Cluster file name:     %s\n", abcinp.clusterfn.c_str());
    printf("Result file directory: %s\n", lmdir.c_str());
    printf("\n");

    RMCluster rmc;  // read xyz, .cluster .....
    rmc.LoadRMCluster(abcinp.clusterfn);
    const vector<pair<string, int> >& components = rmc.GetComponents();
    printf(" -- Cluster Components: from [ %s ] --\n", abcinp.clusterfn.c_str());
    printf(" -------------------------------------------------\n");
    printf("   %-40s %5s\n", "Source", "#");
    printf(" -------------------------------------------------\n");
    int nc = components.size();
    for(int i = 0; i < nc; ++i)
    {
        printf("   %-40s %5d\n", components[i].first.c_str(), components[i].second);
    }
    printf(" -------------------------------------------------\n");
    printf("\n");

    // Generating Clsuters
    printf(" -- Structure Generation --\n");    
    if(abcinp.structuretypes.size() != nc)
    {
        RunStatus::ErrorTermination("The number of cluster components (%d) must be equal to the number of structure descripions (%d).\n",
            nc, abcinp.structuretypes.size());
    }
    StructureGenerator geompool;   ////////////////////////////////// add struc feature
    vector<string> descriptions;
    geompool.Initialization(abcinp.structuretypes, descriptions);
    for(int i = 0; i < nc; ++i)
    {
        printf(" %3d %s: %s.\n", i, components[i].first.c_str(), descriptions[i].c_str());
    }
    printf("\n");

    printf(" -- External Commands --\n");    ///////////////// extenal cmd like cp2k
    int ncmds = abcinp.optcmds.size();
    for(int i = 0; i < ncmds; ++i)
    {
        printf("%s\n", abcinp.optcmds[i].c_str());
    }
    printf("\n");

    Parallelization dispatcher;

    EnergyOrder pool; pool.Clear();
    printf(" -- NWPEsSe Optimization --\n");
    printf("Will perform at most %d calculations.\n", abcinp.maxncals);
    if(abcinp.docoarseopt)
    {
        printf("A coarse optimization will be done for each generated structure.\n"); 
    }
    else
    {
        printf("NO coarse optimization will be done for each generated structure.\n"); 
    }
    printf("Cluster energy boundary for training: (%.8f, %.8f)\n", abcinp.energy_lb, abcinp.energy_ub);
    CalculationEngine eg;

    if(dispatcher.GetNumNodes() > 0)   ////// start opt
    {
        printf("===============================================================\n");
        printf("%5s %15s %10s   %s\n", "#", "Energy", "Time", "State");    
        printf("===============================================================\n");
        for(int i = 0; i < abcinp.maxncals; ++i)
        {        
            bool is_running = false;
            while(!is_running)
            {
                for(int node_idx = 0; node_idx < dispatcher.GetNumNodes(); ++node_idx)
                {
                    if(dispatcher.IsAvailable(node_idx))
                    {
                        dispatcher.SetUnavailable(node_idx);		
                        pthread_t pt;
                        
                        TaskArgs* pargs = new TaskArgs(&dispatcher, node_idx, 
                                i, lmdir, &geompool, &abcinp, &pool, &eg, &rmc);
                        pthread_create(&pt, NULL, DoTaskThread, (void*)pargs);

                        dispatcher.AddThread(pt);
                        is_running = true;     
                        break;
                    }
                }
            }
            const int WaitTime = 1; // 1 Second
            sleep(WaitTime);
        }
        dispatcher.Barrier();                
        printf("===============================================================\n");
    }
    else
    {
        printf("===============================================================\n");
        printf("%5s %15s %10s   %s\n", "#", "Energy", "Time", "State");
        printf("===============================================================\n");
        for(int i = 0; i < abcinp.maxncals; ++i)   //////////////////// now in cycles
        {
            TaskArgs args(&dispatcher, 0, i, lmdir, &geompool, &abcinp, &pool, &eg, &rmc);
            DoTask(args);
        }
        printf("===============================================================\n");
    }
    printf("All the calculations are finished!\n");
    printf("\n");
    
    printf("Reordered from low to high energy:\n");
    pool.ListOrderedEIs(abcinp.energy_lb, abcinp.energy_ub);

    return RunStatus::NormalTermination();
}
