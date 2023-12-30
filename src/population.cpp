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
    auto parents = select_parents();
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

void Population::printPopulationHistory() const {
  std::cout << "Population count: " << population_history.size();
  for (const auto &pop : population_history) {
    std::cout << "Average fitness: " << pop.average_fitness << "\n";
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
  std::partial_sort(population.begin(), population.begin() + children.size(),
                    population.end(),
                    [](Gene const &first, Gene const &second) {
                      return first.get_fitness() < second.get_fitness();
                    });

  std::copy(children.begin(), children.end(), population.begin());

  set_population(population);
}

[[nodiscard]] std::vector<Gene> Population::select_parents() const {
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
