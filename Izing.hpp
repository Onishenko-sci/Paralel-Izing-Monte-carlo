#include <vector>
#include <random>
#include <mutex>

class Izing
{
private:
    // Izing parameters
    std::vector<std::vector<int>> config;
    double J;
    double kT;
    double h;
    int Lattice_size;
    // Threads parameters
    size_t Number_of_threads;
    std::mutex synch;
    double work_time;

    double near_bounds_energy(int i, int j);
    double full_config_energy();

    friend void work(size_t thread_id, Izing &lattice, unsigned int steps);

public:
    Izing(double J, double h = 0);
    void init(const size_t Lattice_Size);
    void MC_simulation(unsigned long int steps, double kT, int Number_of_threads = 4, int warm_steps = 0);
    void show() const;
    double how_long() const { return work_time; };
    ~Izing();
};
