#include "Izing.hpp"
#include <iostream>
#include <thread>
#include <chrono>

double mean(int N_of_exp, double fr, int core)
{
    Izing A(2);
    double *experiments = new double[N_of_exp];
    double sum = 0;
    using namespace std::chrono_literals;
    using namespace std::this_thread;

    for (int i = 0; i < N_of_exp; i++)
    {
        A.init(4);
        A.layered_rand(10e+7, 2.0, fr, core);
        experiments[i] = A.how_long();
        sum += experiments[i];
        sleep_for(1s);
    }

    delete[] experiments;
    return sum / N_of_exp;
}

int main()
{
    Izing A(1);
    A.init(40);

    int cores[] = {1, 2, 3, 4};
    double Temperature[] = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 5.0};
    double fr_rate[] = {2, 4, 16, 64, 128, 0};
    double ll[] = {4,6,8,9};
    int i = 0000;
        using namespace std::chrono_literals;
    using namespace std::this_thread;

    int core = 2;
    for (int i = 0; i < 4; i++)
    {   
        A.init(40);
        std::cout << cores[i] << ";";
        A.layered_chess(10e+7, 1.8,cores[i]);
        std::cout << A.how_long() << ";\n";
      //  sleep_for(1s);
    }
    A.show();
    return 0;
}