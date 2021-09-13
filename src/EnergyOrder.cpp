#include "EnergyOrder.h"
#include <algorithm>
#include <cstdio>

using namespace abcluster;
using namespace std;

EnergyOrder::EnergyOrder(void) 
{
    eis.clear();
}

EnergyOrder::~EnergyOrder(void) 
{
    /* Nothing to do */ 
}

void EnergyOrder::Clear(void)
{
    eis.clear();
}

void EnergyOrder::AddEI(int index, double energy)
{
    eis.push_back(EnergyOrder::EnergyIndex(index, energy));
}

void EnergyOrder::ListOrderedEIs(double energy_lb, double energy_ub)
{
    sort(eis.begin(), eis.end());
    eis.erase(unique(eis.begin(), eis.end()), eis.end());
    printf("===============================================================\n");
    printf("%5s %15s\n", "#", "Energy");   
    printf("===============================================================\n");
    int neis = eis.size();
    for(int i = 0; i < neis; ++i)
    {
        if(eis[i].energy >= energy_lb && eis[i].energy <= energy_ub)
        {
            printf("%5d %15.8f\n", eis[i].index, eis[i].energy);
        }
    }
    printf("===============================================================\n");        
}
