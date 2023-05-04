#include "../Izing.hpp"
#include <iostream>
#include <vector>

int main()
{
    using std::vector;

    int lattice_size = 48;
    int n_experiments = 10;
    double J = 1;
    long int steps = 1e7;
    int Number_of_threads = 4;
    int frame_rate = 1000;

    Izing *Systems[3];
    for (int j = 0; j < 3; j++)
        Systems[j] = new Izing[n_experiments];


    for (int i = 0; i < n_experiments; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Systems[j][i].init(lattice_size);
        }
    }

    double temperature = 5.1;
    double m = 0;
    vector<double> magnetizations;
    magnetizations.resize(3);
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 3; j++)
            magnetizations[j] = 0;

        for (int n = 0; n < n_experiments; n++)
        {
            Systems[0][n].layered_rand(steps, temperature, 1000, 4);
            Systems[1][n].layered_rand(steps, temperature, 1000, 1);
            Systems[2][n].layered_chess(steps, temperature, 1);
            for (int j = 0; j < 3; j++)
            {
                m = Systems[j][n].magnetization() / (lattice_size * lattice_size);
                magnetizations[j] += (m < 0) ? -m : m;
            }
        }

        std::cout << temperature << ';';
        for (int j = 0; j < 3; j++)
            std::cout << magnetizations[j] / n_experiments << ';';
        std::cout << '\n';
        temperature -= 0.1;
    }

    return 0;
}