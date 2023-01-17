#include <vector>
#include <mutex>

class Izing
{
private:
    // Izing parameters
    std::vector<std::vector<int>> config;
    double J;
    double kT;
    double h;
    size_t Lattice_size;

    // Threads parameters
    size_t Number_of_threads;
    std::mutex synch;
    size_t layer_hight;
    friend void simulation(size_t thread_id, Izing &lattice, unsigned int steps);

    double work_time;
    double near_bounds_energy(int i, int j);
public:
    Izing(double J, double h = 0);
    void init(const size_t Lattice_Size);
    void MC_simulation(unsigned long int steps, double kT, int Number_of_threads = 4);
    double config_energy() const;
    double magnetization() const;
    void show() const;
    double how_long() const { return work_time; };
    ~Izing();
};
