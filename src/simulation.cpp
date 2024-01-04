#include "simulation.hpp"
#include <iostream>
#include <random>

void Simulation::run() {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dis;
  for (std::size_t i = 0; i < max_generations; ++i) {
    auto parents = current_population.select_parents_roulette();
    std::vector<Gene> next_gen{};
    for (std::size_t j = 0; j < parents.size() - 1; j += 2) {
      if (j + 1 < parents.size()) {
        if (recombination_probability >= dis(gen)) {
          auto children =
              current_population.recombine({parents.at(j), parents.at(j + 1)});
          next_gen.push_back(children.first);
          next_gen.push_back(children.second);
          current_population.increment_recombination_count();
        } else {
          next_gen.push_back(parents.at(j));
          next_gen.push_back(parents.at(j + 1));
        }
      }
    }
    for (auto &child : next_gen) {
      if (mutation_probability >= dis(gen)) {
        child = current_population.mutate(child);
        current_population.increment_mutation_count();
      }
    }
    history.push_back(current_population);
    Population next_population(next_gen);
    current_population = next_population;
  }
}

void Simulation::print_statistics() const {
  std::ofstream fitness_statistics("fitness_statistics.csv");
  std::ofstream best_statistics("best_statistics.csv");
  std::ofstream mutation_count_statistics("mutation_count_statistics.csv");
  std::ofstream recombination_count_statistics(
      "recombination_count_statistics.csv");
  std::ofstream error_statistics("error_statistics.csv");
  if (fitness_statistics.is_open() && best_statistics.is_open() &&
      mutation_count_statistics.is_open() &&
      recombination_count_statistics.is_open() && error_statistics.is_open()) {
    for (std::size_t i = 1; i < history.size(); ++i) {
      auto const &p = history[i];
      fitness_statistics << i << "," << p.get_average_fitness() << '\n';
      best_statistics << i << "," << p.get_best_gene().get_fitness() << '\n';
      mutation_count_statistics << i << "," << p.get_mutation_count() << '\n';
      recombination_count_statistics << i << "," << p.get_recombination_count()
                                     << '\n';
      error_statistics << i << "," << p.get_mean_square_error() << '\n';
    }
    fitness_statistics.close();
    best_statistics.close();
    mutation_count_statistics.close();
    recombination_count_statistics.close();
    error_statistics.close();
  }
}
