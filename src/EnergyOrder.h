#ifndef    __ENERGYORDER_H__
#define    __ENERGYORDER_H__

#include <vector>
#include <cmath>

namespace abcluster {

using namespace std;

class EnergyOrder {
public:
    EnergyOrder(void);
    ~EnergyOrder(void);

    void Clear(void);
    void AddEI(int index, double energy);
    void ListOrderedEIs(double energy_lb, double energy_ub);

private:
    class EnergyIndex {
    public:
        EnergyIndex(int tindex, double tenergy) : index(tindex), energy(tenergy) {}
        ~EnergyIndex(void) {}
        bool operator < (const EnergyIndex& ei0) const { return ((energy < ei0.energy) && (abs(energy-ei0.energy) > 1.E-4)); }
        bool operator == (const EnergyIndex& ei0) const { return abs(energy-ei0.energy) < 1.E-4; }
        int index;
        double energy;
    };

    vector<EnergyIndex> eis;
};

}

#endif //  __ENERGYORDER_H__
