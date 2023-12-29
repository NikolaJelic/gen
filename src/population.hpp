#pragma  once

#include <array>
#include <utility>
#include <vector>
#include "gene.hpp"

class Population{
    public:
    //  algorithm variables
    static constexpr std::size_t max_generations = 100;
    static constexpr std::size_t population_size = 1000;
    static constexpr float mutation_probability = 0.1f;
    static constexpr float recombination_probability = 0.7f;

    Population();
    void run();
    [[nodiscard]] std::pair<Gene, Gene> recombine(std::pair<Gene, Gene>const& parents) const;
    [[nodiscard]] Gene mutate(Gene const& gene) const;

    private:
    // population tracking
    std::vector<Gene> population{};
    std::vector< std::vector<Gene>> population_history{};

    // statistics
    std::array<int, max_generations> mutation_count{};
    std::array<int, max_generations> recombination_count{};
    std::array<float, max_generations> mean_square_error{};
    std::array<float, max_generations> average_fitness{};
    std::array<Gene, max_generations> best_genes{};
};