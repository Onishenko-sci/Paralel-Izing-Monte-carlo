#include "../Izing.hpp"
#include <iostream>
#include <vector>
void report(std::vector<int> A);

int main()
{
    double J = 1;
    int lattice_size = 16;
    long int steps = 1e4;
    double Temperature = 0.5;
    int Number_of_threads = 1;
    int order_frame_rate = 1;
    int rand_frame_rate = 1000;
    int Number_of_experiments = 1000;
    int max_steps = 1000;
    std::vector<int> Rand_steps;
    std::vector<int> order_steps;

    Izing A(J);
    Izing B(J);

    for (int i = 0; i < Number_of_experiments; i++)
    {
        A.init(lattice_size, i + 123);
        B.init(lattice_size, i + 123);

        int ACount = 0;
        while ((A.config_energy() > -480) && (ACount < max_steps))
        {
            A.layered_chess(steps, Temperature, Number_of_threads);
            ACount++;
        }
        order_steps.push_back(ACount);

        int BCount = 0;
        while ((B.config_energy() > -480) && (BCount < max_steps))
        {
            B.layered_order(steps, Temperature, order_frame_rate, Number_of_threads);
            BCount++;
        }
        Rand_steps.push_back(BCount);
        std::cout << i << ' ' << ACount << ' ' << BCount << '\n';
    }
    std::cout << "Temperature: " << Temperature << '\n';
    std::cout << "Number of experiments: " << Number_of_experiments << "\nOrder Steps\n";
    report(Rand_steps);
    std::cout << "Chess Steps\n";
    report(order_steps);

    //   A.save_Energy("E_order4.txt");
    //   B.save_Energy("E_rand4.txt");

    return 0;
}

void report(std::vector<int> A)
{
    int neud = 0;
    int sum = 0;
    for (auto &&count : A)
    {
        if (count == 1000)
            neud++;
        else
        {
            sum += count;
        }
    }

    std::cout << "Neud: " << double(100 * neud) / A.size() << "%\n";
    std::cout << "Mean: " << double(sum) / (A.size() - neud) << '\n';
}