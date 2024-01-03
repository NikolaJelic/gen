#include "gene.hpp"
#include <algorithm>
#include <cmath>
#include <random>

Gene::chromosome_array Gene::generate_random_chromosome() const {
  chromosome_array ret{};
  std::random_device rd;
  std::default_random_engine engine(rd());
  std::uniform_int_distribution<int> distribution(0, 1);
  for (int &c : ret) {
    c = distribution(engine);
  }
  return ret;
}

void Gene::set_phenotypes() {
  phenotype_x =
      calculate_phenotype(chromosome_to_decimal(chromosome_x), min_x, max_x);
  phenotype_y =
      calculate_phenotype(chromosome_to_decimal(chromosome_y), min_x, max_x);
}

float Gene::calculate_phenotype(float value, float min, float max) const {
  float in_min = 0.0f;
  float in_max = 1023.0f;
  return min + (max - min) * ((value - in_min) / (in_max - in_min));
}

std::size_t
Gene::chromosome_to_decimal(chromosome_array const &binary_array) const {
  std::size_t decimalValue = 0;
  std::size_t size = binary_array.size();

  for (size_t i = 0; i < size; ++i) {
    if (binary_array[i]) {
      decimalValue |= 1ULL << (size - 1 - i);
    }
  }
  return decimalValue;
}

void Gene::calculate_fitness() {

  fitness =
      1.25 * std::pow(1 - phenotype_x, 2) *
          std::exp(-std::pow(phenotype_x, 2) - std::pow(phenotype_y + 1, 2)) -
      10 * (phenotype_x - std::pow(phenotype_x, 5) - std::pow(phenotype_y, 5)) *
          std::exp(
              -((std::pow(phenotype_x, 2) + std::pow(phenotype_y, 2)) / 0.9)) -
      (1.0 / 5) *
          std::exp(-std::pow(phenotype_x + 1, 2) - std::pow(phenotype_y, 2));
  fitness = std::max(0.0f, fitness);
}