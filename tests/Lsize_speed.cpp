#include "../Izing.hpp"
#include <iostream>
#include <vector>

int main()
{
    double J = 1;
    int lattice_size = 2;
    long int steps = 1e8;
    double Temperature = 2.2;
    int Number_of_threads = 4;
    int frame_rate = 1000;
    int core = 4;
    int cores[] = {1, 2, 3, 4};
    int l_sizes[] = {8,16,32,64,128,216,516,1024,2048,4096};

    Izing A(J);
    for (auto &&size : l_sizes)
    {
        A.init(size);
        A.layered_rand(steps, Temperature, 1000, core);
        std::cout << size << ';' << A.how_long();
        A.layered_order(steps, Temperature, 5, core);
        std::cout << ';' << A.how_long();
        A.layered_chess(steps, Temperature, core);
        std::cout << ';' << A.how_long() << '\n';
    }

    return 0;
}