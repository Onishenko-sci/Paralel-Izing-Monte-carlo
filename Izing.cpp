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

void Izing::init(size_t lattice_size)
{
    Lattice_size = lattice_size;
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

double Izing::config_energy() const
{
    double Energy = 0;
    int last = Lattice_size - 1;
    for (int i = 0; i < Lattice_size; i++)
        for (int j = 0; j < Lattice_size; j++)
        {
            // Right
            if (i != last)
                Energy -= J * config[i][j] * config[i + 1][j];

            // Bottom
            if (j != last)
                Energy -= J * config[i][j] * config[i][j + 1];
        }
    return Energy;
}

double Izing::magnetization() const
{
    double Sum = 0;
    for (int i = 0; i < Lattice_size; i++)
        for (int j = 0; j < Lattice_size; j++)
            Sum += config[i][j];
    return Sum;
}

void simulation(size_t thread_id, Izing &izing, unsigned int steps)
{
    // Init distributions for spins
    std::random_device seed_init;
    std::mt19937 gen(seed_init());

    std::uniform_real_distribution<double> distr(0, 1);
    size_t thread_first_layer = thread_id * izing.layer_hight;
    size_t thread_last_layer = thread_id * izing.layer_hight + izing.layer_hight - 1;

    // Rows numbers generated inside of strip for current thread.
    std::uniform_int_distribution<int> random_i(thread_first_layer, thread_last_layer);
    std::uniform_int_distribution<int> random_j(0, izing.Lattice_size - 1);
    /*
        izing.synch.lock();
        cout << "Thread " << thread_id
             << " for Layers " << thread_first_layer + 1 << "-" << thread_last_layer + 1 << " "
             << "started.\n";
        izing.synch.unlock();
    */
    int Spin_i;
    int Spin_j;
    double energy_before_change;
    double dE;
    double R;
    bool change_spin;
    // Metropolis
    for (unsigned long int step = 0; step < steps; step++)
    {
        // Random spin
        Spin_i = random_i(gen);
        Spin_j = random_j(gen);
        // Progress bar
        if (step % (steps / 20) == 0 && thread_id == 0)
            cout << "=" << std::flush;

        // Calculating dE
        // For contiguous rows using Mutex
        if (Spin_i == thread_first_layer || Spin_i == thread_last_layer)
        {
            std::lock_guard<std::mutex> locker(izing.synch);
            energy_before_change = izing.near_bounds_energy(Spin_i, Spin_j);
        }
        else
            energy_before_change = izing.near_bounds_energy(Spin_i, Spin_j);
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

        if (change_spin && (Spin_i == thread_first_layer || Spin_i == thread_last_layer))
        {
            std::lock_guard<std::mutex> locker(izing.synch);
            izing.config[Spin_i][Spin_j] *= -1;
        }
        else if (change_spin)
            izing.config[Spin_i][Spin_j] *= -1;
    }
    // cout << thread_id << " end his work.\n";
    return;
}

void Izing::MC_simulation(unsigned long int steps, double kT_in, int number_of_threads)
{
    unsigned int start_time = clock();
    // Check number_of_threads
    int max_threads = thread::hardware_concurrency();
    if (number_of_threads > max_threads && max_threads)
    {
        cout << "Too much threads!\n"
             << "Maximum number of threads is " << max_threads << '\n';
        work_time = 0;
        return;
    }

    // Initialization Metropolis parameters
    Number_of_threads = number_of_threads;
    kT = kT_in;
    steps = steps / Number_of_threads;
    layer_hight = Lattice_size / Number_of_threads;

    // Threads creation
    vector<thread> threads;
    for (int i = 1; i < Number_of_threads; i++)
        threads.push_back(thread(simulation, i, std::ref(*this), steps));

    // Use main as 0 thread
    simulation(0, *this, steps);

    // Join threads
    for (int i = 0; i < Number_of_threads - 1; i++)
        threads[i].join();

    // Calculate time
    work_time = (clock() - start_time) / ((double)CLOCKS_PER_SEC * Number_of_threads);
    return;
}

Izing::~Izing()
{
}
