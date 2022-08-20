#include <random>

#include "random.h"

std::mt19937& get_prng_engine() {
    thread_local static std::random_device rd{}; // initial seed for the pRNG
    thread_local static std::mt19937 prng{rd()};
    return prng;
}

int random_bit(double p) {
    std::bernoulli_distribution dist{p};
    return dist(get_prng_engine());
}
