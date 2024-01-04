#pragma once

#include "population.hpp"
#include <fstream>
#include <iostream>
#include <vector>


class Simulation {
public:
  Simulation(std::size_t max_gen, std::size_t population_size,  float mutation_probability,
             float recombination_probability)
      : max_generations(max_gen), population_size(population_size), mutation_probability(mutation_probability),
        recombination_probability(recombination_probability) {
    for (auto const &p : current_population.get_population()) {
    }
    history.push_back(current_population);
  }
  void run();
  void print_statistics() const;

private:
  std::vector<Population> history{};
  const std::size_t population_size;
  Population current_population{population_size};
  const std::size_t max_generations;
  const float mutation_probability;
  const float recombination_probability;
};