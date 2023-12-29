
#include "gene.hpp"
#include <iostream>

int main() {
  double max = 6.251637526363513;
  Gene gene{};
  std::cout << gene.get_phenotype().first << " | " << gene.get_phenotype().second;
}
