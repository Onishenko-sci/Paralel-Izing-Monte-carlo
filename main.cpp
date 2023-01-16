#include "Izing.hpp"
#include <iostream>

int main()
{
    Izing A(-1);

    A.init(20);
    A.show();
    A.MC_simulation(20e+6, 1.5, 1);
    std::cout << "Work time: " << A.how_long() << std::endl;

    A.show();
    return 0;
}