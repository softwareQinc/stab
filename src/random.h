#ifndef STAB_RANDOM_H_
#define STAB_RANDOM_H_

#include <random>

namespace stab {
std::mt19937& get_prng_engine();

int random_bit(double p = 0.5);

// template functions need to be defined in the header

template <class T = int>
T random_integer(T a = std::numeric_limits<T>::min(),
                 T b = std::numeric_limits<T>::max()) {
    std::uniform_int_distribution<T> dist(a, b);
    if (a > b) {
        throw std::out_of_range("a > b");
    }
    if constexpr (!std::is_integral_v<T>) {
        throw std::runtime_error("Type must be integral");
    }
    return dist(get_prng_engine());
}

template <class T = double>
T random_real(T a = std::numeric_limits<T>::min(),
              T b = std::numeric_limits<T>::max()) {
    if (a >= b) {
        throw std::out_of_range("a >= b");
    }
    if constexpr (!std::is_floating_point_v<T>) {
        throw std::runtime_error("Type must be floating point");
    }
    std::uniform_real_distribution<T> dist(a, b);
    return dist(get_prng_engine());
}
} /* namespace stab */

#endif /* STAB_RANDOM_H_ */
