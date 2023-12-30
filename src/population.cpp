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
  population_history.push_back(population);
}

void Population::set_population(std::vector<Gene> population) {
  this->population = population;
  // calculate the statistics
  this->average_fitness = calculate_average_fitness();
  this->mean_square_error = calculate_average_error();
  this->best_gene = get_best_gene();
}

[[nodiscard]] Gene Population::get_best_gene() const {
  Gene best = population.at(0);
  for (auto const &gene : population) {
    best = best.get_fitness() > gene.get_fitness() ? best : gene;
  }
  return best;
}

void Population::replacePopulation(std::vector<Gene> children) {
  std::sort(population.begin(), population.end(),
            [](Gene const &first, Gene const &second) {
              return first.get_fitness() < second.get_fitness();
            });

  for (std::size_t i = 0; i < children.size(); ++i) {
    population[i] = children[i];
  }

  set_population(population);
}

void Population::run() {
  // - combine
  // - mutate and crossover
  // - create new population
  // for (std::size_t i = 0; i <  max_generations; ++i) {
  auto parents = select_parents();
  //}
}

[[nodiscard]] std::vector<Gene> Population::select_parents() const {
  std::vector<Gene> parents{};

  // create randomized indices for selection of random non-repeating pairs from
  // the population vector
  std::vector<int> indices(population.size());
  std::iota(indices.begin(), indices.end(), 0);

  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(indices.begin(), indices.end(), g);

  for (int i = 0; i < population.size() - 1; i += 2) {
    auto candidate_1 = population[indices[i]];
    auto candidate_2 = population[indices[i + 1]];
    auto parent = (candidate_1.get_fitness() > candidate_2.get_fitness())
                      ? candidate_1
                      : candidate_2;
    parents.push_back(parent);
  }

  return parents;
}

[[nodiscard]] std::pair<Gene, Gene>
Population::recombine(std::pair<Gene, Gene> const &parents) const {
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

[[nodiscard]] Gene Population::mutate(Gene const &gene) const {
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

[[nodiscard]] float Population::calculate_average_fitness() const {
  float sum = 0.0f;
  for (auto const &gene : population) {
    sum += gene.get_fitness();
  }
  return sum / population.size();
}

[[nodiscard]] float Population::calculate_average_error() const {
  float max = 6.251637526363513f;
  float sum = 0.0f;
  // mean((f_max - fitness_values).^2);
  for (auto const &gene : population) {
    sum += pow(max - gene.get_fitness(), 2);
  }
  return sum / population.size();
}
