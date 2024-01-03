#pragma once

#include <array>
#include <span>
#include <utility>

class Gene {
public:
  constexpr static float min_x = -3.0f; 
  constexpr static float min_y = -3.0f; 
  constexpr static float max_x = 3.0f; 
  constexpr static float max_y = 3.0f; 
  static constexpr std::size_t gene_length = 10;

  using chromosome_array = std::array<int, gene_length>;

  Gene(){
   chromosome_x =  generate_random_chromosome();
   chromosome_y = generate_random_chromosome();
   set_phenotypes();
   calculate_fitness();
  }

  Gene(chromosome_array ch_x, chromosome_array ch_y): chromosome_x(ch_x), chromosome_y(ch_y){
    set_phenotypes();
    calculate_fitness();
  }

  [[nodiscard]] std::pair<chromosome_array, chromosome_array> get_chromosomes() const{
    return {chromosome_x, chromosome_y} ;
  }

   [[nodiscard]] std::pair<float, float> get_phenotype() const{
    return {phenotype_x, phenotype_y};
  }

  [[nodiscard]] float get_fitness() const{
    return fitness;
  }

private:

  [[nodiscard]] chromosome_array generate_random_chromosome() const;
  void set_phenotypes();  
  [[nodiscard]] float calculate_phenotype(float value, float min, float max) const;
  void calculate_fitness();
  [[nodiscard]] std::size_t chromosome_to_decimal(chromosome_array const& binary_array) const;

  chromosome_array chromosome_x{};
  chromosome_array chromosome_y{};
  float phenotype_x{};
  float phenotype_y{};
  float fitness{};

  
};