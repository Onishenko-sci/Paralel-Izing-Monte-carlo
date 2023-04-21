#include "Izing.hpp"
#include <random>
#include <iostream>
#include <math.h>
#include <ctime>
#include <chrono>
#include <thread>

#include <fstream>

using std::cout;
using std::mutex;
using std::thread;
using std::vector;

void Izing::report()
{
    cout << "-----Report----- T:" << kT << "-- fr:" << Frame_rate << "---Latice:" << Lattice_size << "---Proc:" << Number_of_threads << '\n';
    std::cout << "Work time: " << work_time << std::endl;
    std::cout << "Steps with block: " << steps_with_block << std::endl;
    std::cout << "\% of blocked: " << steps_with_block / double(Steps * Number_of_threads) << std::endl;
    cout << "----------------------------------------\n";
}

double Izing::block()
{
    return steps_with_block / double(Steps * Number_of_threads);
}

void Izing::safe_data(const char *name)
{
    std::ofstream A(name);
    for (int i = 0; i < Steps; i++)
    {
        A << i << ";";
        for (int j = 0; j < Number_of_threads; j++)
            A << Magnetizations[j][i] << ";";
        A << "\n";
    }
}

void Izing::write_relative_step()
{
    std::ofstream A("relative_steps.txt");
    for (int i = 0; i < Steps; i++)
    {
        for (int j = 0; j < Number_of_threads; j++)
            A << thread_relative_step[j][i] << ";";
        A << std::endl;
    }
}

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
    // Top
    if (i != 0)
        Energy -= J * config[i][j] * config[i - 1][j];

    // Bot
    if (i != last)
        Energy -= J * config[i][j] * config[i + 1][j];

    // Left
    if (j != 0)
        Energy -= J * config[i][j] * config[i][j - 1];

    // Right
    if (j != last)
        Energy -= J * config[i][j] * config[i][j + 1];

    return Energy;
}

Izing::Izing(double init_J, double init_h)
{
    J = init_J;
    h = init_h;
}

void Izing::init(size_t lattice_size, int seed)
{
    Lattice_size = lattice_size;
    config.resize(Lattice_size);

    if (!seed)
    {
        std::random_device rd;
        seed = rd();
    }
    std::mt19937 generator(seed);

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

void warm(size_t thread_id, Izing &izing, unsigned int steps)
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
    int Last_thread_id = izing.Number_of_threads - 1;
    bool change_spin;
    // Metropolis
    for (unsigned long int step = 0; step < steps; step++)
    {
        // Random spin
        Spin_i = random_i(gen);
        Spin_j = random_j(gen);
        // Progress bar
        if (step % (steps / 5) == 0)
            cout << thread_id << std::flush;

        // Calculating dE
        // For contiguous rows using Mutex
        if (Spin_i == thread_first_layer && thread_id != 0)
        {
            std::lock_guard<mutex> locker(izing.Mutexi[thread_id - 1]);
            energy_before_change = izing.near_bounds_energy(Spin_i, Spin_j);
        }
        else if (Spin_i == thread_last_layer && thread_id != Last_thread_id)
        {
            std::lock_guard<mutex> locker(izing.Mutexi[thread_id]);
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

        if (change_spin)
        {
            if (Spin_i == thread_first_layer && thread_id != 0)
            {
                std::lock_guard<mutex> locker(izing.Mutexi[thread_id - 1]);
                izing.config[Spin_i][Spin_j] *= -1;
            }
            else if (Spin_i == thread_last_layer && thread_id != Last_thread_id)
            {
                std::lock_guard<mutex> locker(izing.Mutexi[thread_id]);
                izing.config[Spin_i][Spin_j] *= -1;
            }
            else
                izing.config[Spin_i][Spin_j] *= -1;
        }
    }
    // cout << thread_id << " end his work.\n";
    return;
}

void Izing::warming_up(unsigned long int steps, double kT_in, int number_of_threads)
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
    Mutexi = new mutex[Number_of_threads - 1];
    for (int i = 1; i < Number_of_threads; i++)
    {
        threads.push_back(thread(warm, i, std::ref(*this), steps));
    }

    // Use main as 0 thread
    warm(0, *this, steps);

    // Join threads
    for (int i = 0; i < Number_of_threads - 1; i++)
        threads[i].join();

    cout << std::endl;
    // Calculate time
    work_time = (clock() - start_time) / ((double)CLOCKS_PER_SEC * Number_of_threads);
    delete[] Mutexi;
    return;
}

void rand_spin(size_t thread_id, Izing &izing, unsigned int steps, Barrier &bar)
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

    int Spin_i;
    int Spin_j;
    double energy_before_change;
    double dE;
    double R;
    int Last_thread_id = izing.Number_of_threads - 1;
    bool change_spin;

    // cout << thread_id << ": " << thread_first_layer << " to  " << thread_last_layer << " started\n";

    int Current_M = 0;
    int Current_E = 0;
    for (int i = thread_first_layer; i <= thread_last_layer; i++)
    {
        for (int j = 0; j < izing.Lattice_size; j++)
        {
            Current_M += izing.config[i][j];
            Current_E += izing.near_bounds_energy(i, j);
            // не учитываем нижнюю связь нижнего слоя. она учитывается как верхняя в следующем слое.
            if (i == thread_last_layer && thread_id != Last_thread_id)
                Current_E += izing.J * izing.config[i][j] * izing.config[i + 1][j];
        }
    }

    // Metropolis
    for (unsigned long int step = 0; step < steps; step++)
    {
        // Random spin
        Spin_i = random_i(gen);
        Spin_j = random_j(gen);
        // Progress bar

        //  if (step % (steps / 20) == 0)
        //     cout << thread_id << std::flush;

        if (izing.Frame_rate && step % izing.Frame_rate == 0)
            bar.Wait();
        // Calculating dE
        // For contiguous rows using Mutex

        if (Spin_i == thread_first_layer && thread_id != 0)
        {
            izing.Mutexi[thread_id - 1].lock();
            energy_before_change = izing.near_bounds_energy(Spin_i, Spin_j);
            izing.steps_with_block++;
        }
        else if (Spin_i == thread_last_layer && thread_id != Last_thread_id)
        {
            izing.Mutexi[thread_id].lock();
            energy_before_change = izing.near_bounds_energy(Spin_i, Spin_j);
            izing.steps_with_block++;
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

        if (change_spin)
        {
            izing.config[Spin_i][Spin_j] *= -1;
            Current_M += 2 * izing.config[Spin_i][Spin_j];
            Current_E += dE;
        }

        if (Spin_i == thread_first_layer && thread_id != 0)
            izing.Mutexi[thread_id - 1].unlock();
        else if (Spin_i == thread_last_layer && thread_id != Last_thread_id)
            izing.Mutexi[thread_id].unlock();

        izing.Magnetizations[thread_id][step] = Current_M;
        izing.Energy[thread_id][step] = Current_E;
    }

    // cout << thread_id << " end his work.\n";
    return;
}

void Izing::layered_rand(unsigned long long int steps, double kT_in, int frame_rate, int number_of_threads)
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
    Frame_rate = frame_rate;
    kT = kT_in;
    steps = steps / Number_of_threads;
    Steps = steps;
    steps_with_block = 0;
    layer_hight = Lattice_size / Number_of_threads;

    Magnetizations.resize(Number_of_threads);
    Energy.resize(Number_of_threads);
    my_step.resize(Number_of_threads);
    thread_relative_step.resize(Number_of_threads);

    for (int i = 0; i < Number_of_threads; i++)
    {
        Magnetizations[i].resize(steps);
        Energy[i].resize(steps);
        //        thread_relative_step[i].resize(steps);
    }

    // Threads creation
    vector<thread> threads;
    Barrier Bar(Number_of_threads);
    Mutexi = new mutex[Number_of_threads - 1];
    for (int i = 1; i < Number_of_threads; i++)
        threads.push_back(thread(rand_spin, i, std::ref(*this), steps, std::ref(Bar)));
    // Use main as 0 thread
    rand_spin(0, *this, steps, Bar);
    // Join threads
    for (int i = 0; i < Number_of_threads - 1; i++)
        threads[i].join();
    // cout << std::endl;
    // Calculate time
    work_time = (clock() - start_time) / ((double)CLOCKS_PER_SEC * Number_of_threads);
    delete[] Mutexi;
    return;
}

void chess(size_t thread_id, Izing &izing, unsigned int steps, Barrier &bar)
{
    std::random_device seed_init;
    std::mt19937 gen(seed_init());
    std::uniform_real_distribution<double> distr(0, 1);
    size_t thread_first_layer = thread_id * izing.layer_hight;
    size_t thread_last_layer = thread_id * izing.layer_hight + izing.layer_hight - 1;
    double energy_before_change;
    double dE;
    double R;
    int Last_thread_id = izing.Number_of_threads - 1;
    bool change_spin;
    int last_j = izing.Lattice_size - 1;

    int Current_M = 0;
    int Current_E = 0;
    for (int i = thread_first_layer; i <= thread_last_layer; i++)
    {
        for (int j = 0; j < izing.Lattice_size; j++)
        {
            Current_M += izing.config[i][j];
            Current_E += izing.near_bounds_energy(i, j);
            if (i == thread_last_layer && thread_id != Last_thread_id)
                Current_E += izing.J * izing.config[i][j] * izing.config[i + 1][j];
        }
    }

    bool black = false;

    for (unsigned long int step = 0; step < steps; step++)
    {
        // Progress bar
        // if (step % (steps / 20) == 0)
        //      cout << thread_id << std::flush;

        // Black
        for (int i = thread_first_layer; i <= thread_last_layer; i++)
        {
            black = i % 2;
            for (int j = black; j <= last_j; j += 2)
            {
                energy_before_change = izing.near_bounds_energy(i, j);
                dE = -2 * energy_before_change;
                change_spin = false;
                if (dE > 0)
                {
                    R = std::exp(-dE / izing.kT);
                    if (R > distr(gen))
                    {
                        izing.config[i][j] *= -1;
                        Current_M += 2 * izing.config[i][j];
                        Current_E += dE;
                    }
                }
                else
                {
                    izing.config[i][j] *= -1;
                    Current_M += 2 * izing.config[i][j];
                    Current_E += dE;
                }
            }
        }

        bar.Wait();

        // White
        for (int i = thread_first_layer; i <= thread_last_layer; i++)
        {
            black = i % 2;
            for (int j = !black; j <= last_j; j += 2)
            {
                energy_before_change = izing.near_bounds_energy(i, j);
                dE = -2 * energy_before_change;
                change_spin = false;
                if (dE > 0)
                {
                    R = std::exp(-dE / izing.kT);
                    if (R > distr(gen))
                    {
                        izing.config[i][j] *= -1;
                        Current_M += 2 * izing.config[i][j];
                        Current_E += dE;
                    }
                }
                else
                {
                    izing.config[i][j] *= -1;
                    Current_M += 2 * izing.config[i][j];
                    Current_E += dE;
                }
            }
        }
        izing.Magnetizations[thread_id][step] = Current_M;
        izing.Energy[thread_id][step] = Current_E;
        bar.Wait();
    }
}

void Izing::layered_chess(unsigned long int steps, double kT_in, int number_of_threads)
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
    steps = steps / (Lattice_size * Lattice_size);
    Steps = steps;
    layer_hight = Lattice_size / Number_of_threads;
    Magnetizations.resize(Number_of_threads);
    Energy.resize(Number_of_threads);

    for (int i = 0; i < Number_of_threads; i++)
    {
        Magnetizations[i].resize(steps);
        Energy[i].resize(steps);
    }

    // Threads creation
    vector<thread> threads;
    Barrier Bar(Number_of_threads);
    for (int i = 1; i < Number_of_threads; i++)
        threads.push_back(thread(chess, i, std::ref(*this), steps, std::ref(Bar)));
    // Use main as 0 thread
    chess(0, *this, steps, Bar);
    // Join threads
    for (int i = 0; i < Number_of_threads - 1; i++)
        threads[i].join();
    // cout << std::endl;
    //  Calculate time
    work_time = (clock() - start_time) / ((double)CLOCKS_PER_SEC * Number_of_threads);
    return;
}

Izing::~Izing()
{
}
