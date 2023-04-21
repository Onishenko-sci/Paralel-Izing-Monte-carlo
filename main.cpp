#include "Izing.hpp"
#include <iostream>
#include <thread>
#include <chrono>

int main()
{
  Izing A(1);
  A.init(40);

  int cores[] = {1, 2, 3, 4};
  double Temperature[] = {1, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 3, 5};
  double fr_rate[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 60, 90, 150, 1000, 10000};

  using namespace std::chrono_literals;
  using namespace std::this_thread;

  for (auto &&Temp : Temperature)
  {
    A.init(48);
    double mag = 0;
    A.layered_rand(1e8,Temp,1000,4);
    mag = A.magnetization();
    A.init(48);
    A.layered_chess(1e8,Temp,4);
    std::cout << Temp << ';' << mag/(48*48) << ';' << A.magnetization()/(48*48) << ';' << std::endl;
  }
  

  return 0;
}