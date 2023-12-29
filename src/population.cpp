#include "population.hpp"
#include <random>
#include <utility>

Population::Population() {
    for(int i = 0; i < population_size; ++i){
        population.push_back({});
    }
    population_history.push_back(population);
}


void Population::run() {
  

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