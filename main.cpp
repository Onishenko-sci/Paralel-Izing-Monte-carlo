#include "Izing.hpp"
#include <iostream>

int main()
{
    Izing A(1);
    A.init(48);

    int cores[] = {1, 2, 3, 4};
    double fr_rate[] = { 0.5, 1.0, 2.5, 5.0, 10};
    for (auto &&T : fr_rate)
    {
        for (auto &&core : cores)
        {
            A.simulation_block(10e+6, T, 100, core);
            std::cout << A.how_long() << ";";
            A.simulation(10e+6, T, 100, core);
            std::cout << A.how_long() << ";\n";
        }
        std::cout << std::endl;
    }

    /*
        long int steps = 100e+6;
        std::cout << "\n---------- " << steps << " Metropolis steps compalition ----------\n";
        for (auto &&i : cores)
        {
            A.init(48);
            A.warming_up(steps, 0.9, i);
            std::cout << "\tWork time for " << i << " cores: " << A.how_long() << std::endl;
        }
        */
    return 0;
}