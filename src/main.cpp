#include "simulation.hpp"
#include <chrono>
#include <iostream>

int main() {
  // Start the clock
  auto start_time = std::chrono::high_resolution_clock::now();
  Simulation simulation(100, 25000, 0.4, 0.9);
  simulation.run();
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(
      end_time - start_time);
  std::cout << "Simulation took " << duration.count() << " seconds."
            << std::endl;
  simulation.print_statistics();
  return 0;
}
