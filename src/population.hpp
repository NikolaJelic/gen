#pragma once

#include "gene.hpp"
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

class Population {
public:
  //  algorithm variables

  void print_statistics() const;
  Population(std::size_t population_size);
  Population(std::vector<Gene> const &new_population);
  [[nodiscard]] std::span<Gene> get_population() {
    return std::span<Gene>(population);
  }

  [[nodiscard]] std::pair<Gene, Gene>
  recombine(std::pair<Gene, Gene> const &parents);
  [[nodiscard]] Gene mutate(Gene const &gene);
  /// Uses the binary tournament selection for choosing parents
  [[nodiscard]] std::vector<Gene> select_parents_tournament() const;
  /// Uses a roulette wheel for choosing parents
  [[nodiscard]] std::vector<Gene> select_parents_roulette();
  [[nodiscard]] float calculate_average_fitness() const;
  [[nodiscard]] float calculate_average_error() const;
  [[nodiscard]] Gene get_best_gene() const;
  void set_population(std::vector<Gene> population);
  void replace_population(std::vector<Gene> children);

  constexpr void increment_mutation_count() { ++mutation_count; }
  constexpr void increment_recombination_count() { ++recombination_count; }

  [[nodiscard]] std::size_t get_mutation_count() const {
    return mutation_count;
  }
  [[nodiscard]] std::size_t get_recombination_count() const {
    return recombination_count;
  }
  [[nodiscard]] float get_mean_square_error() const {
    return mean_square_error;
  }
  [[nodiscard]] float get_average_fitness() const { return average_fitness; }

private:
  std::size_t population_size;
  // population tracking
  std::vector<Gene> population{};

  std::default_random_engine engine;

  // statistics
  std::size_t mutation_count{};
  std::size_t recombination_count{};
  float mean_square_error{};
  float average_fitness{};
  Gene best_gene{};
};