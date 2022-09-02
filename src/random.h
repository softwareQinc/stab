#pragma once

#include <random>

namespace stab {
    std::mt19937 &get_prng_engine();

    int random_bit(double p = 0.5);

// template functions need to be defined in the header

    template<class T = int>
    T random_integer(T a = std::numeric_limits<T>::min(),
                     T b = std::numeric_limits<T>::max()) {
        std::uniform_int_distribution<T> dist{a, b};
        return dist(get_prng_engine());
    }

    template<class T = double>
    T random_real(T a = std::numeric_limits<T>::min(),
                  T b = std::numeric_limits<T>::max()) {
        std::uniform_real_distribution<T> dist{a, b};
        return dist(get_prng_engine());
    }
} // namespace stab
