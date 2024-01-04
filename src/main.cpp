
#include "simulation.hpp"

int main() {
Simulation simulation(100,0.1,0.4);
simulation.run();
simulation.print_statistics();
}
