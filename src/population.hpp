#pragma once

#include "gene.hpp"
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

class Population {
public:
  //  algorithm variables
  static constexpr std::size_t max_generations = 100;
  static constexpr std::size_t population_size = 1000;
  static constexpr float mutation_probability = 0.3f;
  static constexpr float recombination_probability = 0.8f;

  void printPopulationHistory() const;
  Population();
  Population(const Population& other);
  void run();

private:
  [[nodiscard]] std::pair<Gene, Gene>
  recombine(std::pair<Gene, Gene> const &parents) const;
  [[nodiscard]] Gene mutate(Gene const &gene) const;
  [[nodiscard]] std::vector<Gene> select_parents() const;
  [[nodiscard]] float calculate_average_fitness() const;
  [[nodiscard]] float calculate_average_error() const;
  [[nodiscard]] Gene get_best_gene() const;
  void set_population(std::vector<Gene> population);
  void replace_population(std::vector<Gene> children);

  // population tracking
  std::vector<Gene> population{};
  std::vector<Population> population_history{};

  // statistics
  std::size_t mutation_count{};
  std::size_t recombination_count{};
  float mean_square_error{};
  float average_fitness{};
  Gene best_gene{};
};