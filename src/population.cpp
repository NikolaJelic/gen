#include "population.hpp"
#include <algorithm>
#include <fstream>
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

Population::Population(const Population &other) {
  // Perform a deep copy
  population = other.population;
  mutation_count = other.mutation_count;
  recombination_count = other.recombination_count;
  mean_square_error = other.mean_square_error;
  average_fitness = other.average_fitness;
  best_gene = other.best_gene;
}

void Population::run() {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dis;
  for (std::size_t i = 0; i < max_generations; ++i) {
    auto parents = select_parents_roulette();
    std::vector<Gene> next_gen{};
    for (std::size_t j = 0; j < parents.size(); j += 2) {
      if (i + 1 < parents.size()) {
        if (recombination_probability >= dis(gen)) {
          auto children = recombine({parents[i], parents[i + 1]});
          next_gen.push_back(children.first);
          next_gen.push_back(children.second);
          ++recombination_count;
        } else {
          next_gen.push_back(parents[i]);
          next_gen.push_back(parents[i + 1]);
        }
      }
    }
    for (auto &child : next_gen) {
      if (mutation_probability >= dis(gen)) {
        child = mutate(child);
        ++mutation_count;
      }
    }
    replace_population(next_gen);
  }
}

void Population::print_statistics() const {
  std::ofstream fitness_statistics("fitness_statistics.csv");
  std::ofstream best_statistics("best_statistics.csv");
  std::ofstream mutation_count_statistics("mutation_count_statistics.csv");
  std::ofstream recombination_count_statistics(
      "recombination_count_statistics.csv");
  std::ofstream error_statistics("error_statistics.csv");
  if (fitness_statistics.is_open() && best_statistics.is_open() &&
      mutation_count_statistics.is_open() &&
      recombination_count_statistics.is_open() && error_statistics.is_open()) {
    for (std::size_t i = 0; i < population_history.size(); ++i) {
      auto const &p = population_history[i];
      fitness_statistics << i << "," << p.average_fitness << '\n';
      best_statistics << i << "," << p.best_gene.get_fitness() << '\n';
      mutation_count_statistics << i << "," << p.mutation_count << '\n';
      recombination_count_statistics << i << "," << p.recombination_count
                                     << '\n';
      error_statistics << i << "," << p.mean_square_error << '\n';
    }
    fitness_statistics.close();
    best_statistics.close();
    mutation_count_statistics.close();
    recombination_count_statistics.close();
    error_statistics.close();
  }
}

void Population::set_population(std::vector<Gene> population) {
  this->population = population;
  // calculate the statistics
  this->average_fitness = calculate_average_fitness();
  this->mean_square_error = calculate_average_error();
  this->best_gene = get_best_gene();
  population_history.push_back(*this);
}

[[nodiscard]] Gene Population::get_best_gene() const {
  Gene best = population.at(0);
  for (auto const &gene : population) {
    best = best.get_fitness() > gene.get_fitness() ? best : gene;
  }
  return best;
}

void Population::replace_population(std::vector<Gene> children) {
  // Sort the population to be replaced
  if (children.size() < population.size()) {
    std::sort(population.begin(), population.end(),
              [](const Gene &first, const Gene &second) {
                return first.get_fitness() < second.get_fitness();
              });
  }
  // Replace the worst individuals with the new children
  std::copy(children.begin(), children.end(), population.begin());
  set_population(population);
}

[[nodiscard]] std::vector<Gene> Population::select_parents_tournament() const {
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
}

[[nodiscard]] std::vector<Gene> Population::select_parents_roulette() const {
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
  for (size_t i = 0; i + 1 < population.size(); i += 2) {
    float spin1 = dis(g);
    float spin2 = dis(g);

    auto accumulate_fitness = [](float sum, const Gene &gene) {
      return sum + gene.get_fitness();
    };

    auto parent1 = std::find_if(population.begin(), population.end(),
                                [spin1, acc = 0.0f](const Gene &gene) mutable {
                                  acc += gene.get_fitness();
                                  return acc >= spin1;
                                });

    // Adjust the range for searching parent2 to avoid out-of-bounds access
    auto parent2 = std::find_if(parent1 + 1, population.end(),
                                [spin2, acc = 0.0f](const Gene &gene) mutable {
                                  acc += gene.get_fitness();
                                  return acc >= spin2;
                                });

    if (parent2 != population.end()) {
      parents.push_back(*parent1);
      parents.push_back(*parent2);
    }
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
  for (auto const &gene : population) {
    sum += pow(max - gene.get_fitness(), 2);
  }
  return sum / population.size();
}
