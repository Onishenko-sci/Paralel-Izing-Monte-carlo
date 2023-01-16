#include "Izing.hpp"
#include <random>
#include <iostream>
#include <math.h>
#include <ctime>
#include <chrono>
#include <thread>

using std::cout;
using std::thread;
using std::vector;

void Izing::show() const
{
    std::cout << "---Izing config---\n";
    for (int i = 0; i < Lattice_size; i++)
    {
        for (int j = 0; j < Lattice_size; j++)
            if (config[i][j] == -1)
                std::cout << "0 ";
            else
                std::cout << ". ";
        std::cout << std::endl;
    }
}

double Izing::near_bounds_energy(int i, int j)
{
    double Energy = 0;
    int last = Lattice_size - 1;
    // Left
    if (i != 0)
        Energy -= J * config[i][j] * config[i - 1][j];

    // Right
    if (i != last)
        Energy -= J * config[i][j] * config[i + 1][j];

    // Top
    if (j != 0)
        Energy -= J * config[i][j] * config[i][j - 1];

    // Bottom
    if (j != last)
        Energy -= J * config[i][j] * config[i][j + 1];

    return Energy;
}

Izing::Izing(double init_J, double init_h)
{
    J = init_J;
    h = init_h;
}

void Izing::init(size_t Lattice_Size)
{
    Lattice_size = Lattice_Size;
    config.resize(Lattice_size);

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> distr(0, 1);
    for (int i = 0; i < Lattice_size; i++)
    {
        config[i].reserve(Lattice_size);
        for (int j = 0; j < Lattice_size; j++)
            config[i][j] = 2 * distr(generator) - 1;
    }
}

void work(size_t thread_id, Izing &izing, unsigned int steps)
{
    std::random_device seed_init;
    std::mt19937 gen(seed_init());

    std::uniform_real_distribution<double> distr(0, 1);
    std::uniform_int_distribution<int> random_i(0, izing.Lattice_size - 1);
    std::uniform_int_distribution<int> random_j(0, izing.Lattice_size - 1);
    cout << "Thread " << thread_id << " started.\n";

    int Spin_i;
    int Spin_j;
    double energy_before_change;
    double dE;
    double R;
    bool change_spin;
    for (unsigned long int step = 0; step < steps; step++)
    {
        Spin_i = random_i(gen);
        Spin_j = random_j(gen);

        izing.synch.lock();
        energy_before_change = izing.near_bounds_energy(Spin_i, Spin_j);
        izing.synch.unlock();

        dE = -2 * energy_before_change;

        change_spin = false;
        if (dE > 0)
        {
            R = std::exp(-dE / izing.kT);
            if (R > distr(gen))
                change_spin = true;
        }
        else
            change_spin = true;

        izing.synch.lock();
        if (change_spin)
            izing.config[Spin_i][Spin_j] *= -1;
        izing.synch.unlock();
    }

    cout << "Done " << thread_id << "\n";
    return;
}

void Izing::MC_simulation(unsigned long int steps, double kT_in, int number_of_threads, int warm_steps)
{
    unsigned int start_time = clock();
    cout << "real max is " << thread::hardware_concurrency() << '\n';
    int max_threads = thread::hardware_concurrency();
    if (number_of_threads > max_threads && max_threads)
    {
        cout << "Too much threads!\n"
             << "Maximum number of threads is " << max_threads << '\n';
        work_time = 0;
        return;
    }

    Number_of_threads = number_of_threads;
    kT = kT_in;
    steps = steps / number_of_threads;

    vector<thread> threads;

    for (int i = 1; i < Number_of_threads; i++)
        threads.push_back(thread(work, i, std::ref(*this), steps));

    work(0, *this, steps);

    for (int i = 0; i < Number_of_threads - 1; i++)
        threads[i].join();

    work_time = (clock() - start_time) / (double)CLOCKS_PER_SEC;
    return;
}

Izing::~Izing()
{
}
