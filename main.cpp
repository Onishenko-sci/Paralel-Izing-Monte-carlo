#include "Izing.hpp"
#include <iostream>

int main()
{
    Izing A(-1);
    A.init(48);
    A.MC_simulation(1e+6, 0.9, 4);
    A.show();
    std::cout << "Work time: " << A.how_long() << std::endl
              << "Magnetization: " << A.magnetization() << std::endl
              << "Energy: " << A.config_energy() << std::endl;

    int cores[] = {1,2,3,4};
    long int steps = 100e+6;
    std::cout << "\n---------- " << steps << " Metropolis steps compalition ----------\n";
    for (auto &&i : cores)
    {
        A.init(48);
        A.MC_simulation(steps, 0.9, i);
        std::cout << "\tWork time for " << i << " cores: " << A.how_long() << std::endl;
    }
    return 0;
}