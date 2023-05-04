#include "Izing.hpp"

int main()
{
  double J = 1;
  int lattice_size = 48;
  long int steps = 1e8;
  double Temperature = 2.2;
  int Number_of_threads = 4;
  int frame_rate = 1000;

  Izing A(J);
  A.init(lattice_size);
  A.warming_up(steps / 10, Temperature, Number_of_threads);
  A.layered_rand(steps, Temperature, frame_rate, Number_of_threads);
  A.show();
  A.report();
  return 0;
}