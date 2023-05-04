#include "../Izing.hpp"
#include <iostream>
#include <vector>

int main()
{
    double J = 1;
    int lattice_size = 48;
    long int steps = 1e8;
    double Temperature = 2.2;
    int Number_of_threads = 4;
    int n_Experiments = 3;
    int frame_rate = 1000;
    int cores[] = {1, 2, 3, 4};
    int l_sizes[] = {4, 8, 16, 32, 64, 128, 216, 516, 1024, 2048};

    Izing A(J);
    A.init(lattice_size);

    A.warming_up(steps / 10, Temperature, Number_of_threads);
    for (auto &&core : cores)
    {
        double sum = 0;
        for (int i = 0; i < n_Experiments; i++)
        {
            A.layered_rand(steps, Temperature, 1000, core);
            sum += A.how_long();
        }
        std::cout << core << ';' << sum / n_Experiments;

        sum = 0;
        for (int i = 0; i < n_Experiments; i++)
        {
            A.layered_order(steps, Temperature, 1, core);
            sum += A.how_long();
        }
        std::cout << ';' << sum / n_Experiments;

        sum = 0;
        for (int i = 0; i < n_Experiments; i++)
        {
            A.layered_chess(steps, Temperature, core);
            sum += A.how_long();
        }
        std::cout << ';' << sum / n_Experiments;

        std::cout << ';' << '\n';
    }
    return 0;
}
