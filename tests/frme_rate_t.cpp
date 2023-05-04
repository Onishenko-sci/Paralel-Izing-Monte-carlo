
#include "../Izing.hpp"
#include <iostream>
#include <vector>

int main()
{
    double J = 1;
    int lattice_size = 12;
    long int steps = 1e7;
    double Temperature = 2.2;
    int Number_of_threads = 4;

    int cores[] = {1, 2, 3, 4};
    double fr_rate[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 60, 90, 150, 1000, 10000};

    Izing A(J);
    A.init(lattice_size);
    A.warming_up(steps / 10, Temperature, Number_of_threads);
    for (auto &&fr : fr_rate)
    {

        A.layered_rand(steps, Temperature, fr, Number_of_threads);
        std::cout << fr << ';' << A.how_long() << '\n';
    }

    return 0;
}
