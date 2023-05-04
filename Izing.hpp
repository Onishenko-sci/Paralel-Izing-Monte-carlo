#include <vector>
#include <mutex>
#include <iostream>
#include "Barrier.hpp"

#define NUM_OF_THREADS 4

class Izing
{
private:
    // Izing parameters
    std::vector<std::vector<int>> config;
    double J;
    double kT;
    double h;
    size_t Lattice_size;
    std::vector<std::vector<int>> Magnetizations;
    std::vector<std::vector<int>> Energy;
    std::vector<size_t> my_step;
    
    std::vector<std::vector<int>> thread_relative_step;

    // Threads parameters
    size_t Number_of_threads;
    std::mutex *Mutexi;
    std::mutex observed;
    int steps_with_block;
    size_t layer_hight;
    size_t Frame_rate;
    friend void warm(size_t thread_id, Izing &lattice, unsigned int steps);
    friend void rand_spin(size_t thread_id, Izing &lattice, unsigned int steps, Barrier &bar);

    friend void ordered_spin(size_t thread_id, Izing &lattice, unsigned int steps, Barrier &bar);

    friend void chess(size_t thread_id, Izing &lattice, unsigned int steps, Barrier& bar);
    long int Steps;
    double work_time;
    double near_bounds_energy(int i, int j);

public:
    Izing(double J = 1, double h = 0);
    void init(const size_t Lattice_Size, int seed = 0);
    void save_Energy(const char *name);
    void save_Magnetization(const char *name);

    void warming_up(unsigned long int steps, double kT, int Number_of_threads = 4);
    void layered_rand(unsigned long long int steps, double kT, int frame_rate = 1000, int Number_of_threads = 4);
    void layered_order(unsigned long long int steps, double kT, int frame_rate = 1, int Number_of_threads = 4);

    void layered_chess(unsigned long int steps, double kT, int Number_of_threads = 4);
    void write_relative_step();

    double config_energy() const;
    double magnetization() const;
    void show() const;
    double how_long() const { return work_time; };
    double block();
    void report();
    ~Izing();
};
