#include "population.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

Population::Population() {
  std::vector<Gene> new_population{};
  for (int i = 0; i < population_size; ++i) {
    new_population.push_back({});
  }
  set_population(new_population);
}

Population::Population(std::vector<Gene> const &new_population) {
  set_population(new_population);
}

void Population::set_population(std::vector<Gene> population) {
  this->population = std::move(population);
  average_fitness = calculate_average_fitness();
  mean_square_error = calculate_average_error();
  best_gene = get_best_gene();
  mutation_count = 0;
  recombination_count = 0;
}

Gene Population::get_best_gene() const {
  Gene best = population.at(0);
  for (auto const &gene : population) {
    best = best.get_fitness() > gene.get_fitness() ? best : gene;
  }
  return best;
}

/*
 std::vector<Gene> Population::select_parents_tournament() const {
  std::vector<Gene> parents{};

  // Create randomized indices for selection of random non-repeating pairs from
  // the population vector
  std::vector<int> indices(population.size());
  std::iota(indices.begin(), indices.end(), 0);

  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(indices.begin(), indices.end(), g);

  for (int i = 0; i < population.size(); i += 2) {
    if (i + 1 < population.size()) {
      auto candidate_1 = population[indices[i]];
      auto candidate_2 = population[indices[i + 1]];
      auto parent = (candidate_1.get_fitness() > candidate_2.get_fitness())
                        ? candidate_1
                        : candidate_2;
      parents.push_back(parent);
    }
  }

  return parents;
} */

std::vector<Gene> Population::select_parents_roulette() const {
  std::vector<Gene> parents;

  // Calculate total fitness
  float total_fitness = std::accumulate(
      population.begin(), population.end(), 0.0f,
      [](float sum, const Gene &gene) { return sum + gene.get_fitness(); });

  if (total_fitness == 0.0f) {
    // Handle division by zero or other appropriate error handling
    return parents;
  }

  // Perform Roulette Wheel Selection
  std::random_device rd;
  std::mt19937 g{rd()};
  std::uniform_real_distribution<float> dis(0.0f, total_fitness);

  // Adjust the loop condition to ensure pairs of parents can be selected
  for (size_t i = 0; i < population.size(); ++i) {
    float spin = dis(g);

    auto accumulate_fitness = [](float sum, const Gene &gene) {
      return sum + gene.get_fitness();
    };

    auto accumulator = 0.0f;
    auto parent = std::find_if(population.begin(), population.end(),
                               [&](const Gene &gene) {
                                 accumulator += gene.get_fitness();
                                 return accumulator >= spin;
                               });
    parents.push_back(*parent);
  }
  return parents;
}
std::pair<Gene, Gene>
Population::recombine(std::pair<Gene, Gene> const &parents) {
  std::random_device rd;
  std::default_random_engine engine(rd());
  std::uniform_int_distribution<int> distribution(1, parents.first.gene_length);
  std::size_t pivot = distribution(engine);

  std::pair<Gene::chromosome_array, Gene::chromosome_array> first =
      parents.first.get_chromosomes();
  std::pair<Gene::chromosome_array, Gene::chromosome_array> second =
      parents.second.get_chromosomes();

  for (int i = pivot; i < parents.first.gene_length; ++i) {
    first.first[i] = parents.second.get_chromosomes().first[i];
    first.second[i] = parents.second.get_chromosomes().second[i];

    second.first[i] = parents.first.get_chromosomes().first[i];
    second.second[i] = parents.first.get_chromosomes().second[i];
  }

  Gene first_child(first.first, first.second);
  Gene second_child(second.first, second.second);

  return {first_child, second_child};
}

Gene Population::mutate(Gene const &gene) {
  std::random_device rd;
  std::default_random_engine engine(rd());
  std::uniform_int_distribution<int> distribution(0, 1);
  Gene::chromosome_array ch_x = gene.get_chromosomes().first;
  Gene::chromosome_array ch_y = gene.get_chromosomes().second;
  for (int i = 0; i < ch_x.size(); ++i) {
    if (distribution(engine)) {
      ch_x[i] = !ch_x[i];
    }
  }
  for (int i = 0; i < ch_y.size(); ++i) {
    if (distribution(engine)) {
      ch_y[i] = !ch_y[i];
    }
  }
  Gene ret(ch_x, ch_y);
  return ret;
}

float Population::calculate_average_fitness() const {
  float sum = 0.0f;
  for (auto const &gene : population) {
    sum += gene.get_fitness();
  }
  return sum / population.size();
}

float Population::calculate_average_error() const {
  float max = 6.251637526363513f;
  float sum = 0.0f;
  for (auto const &gene : population) {
    sum += pow(max - gene.get_fitness(), 2);
  }
  return sum / population.size();
}
