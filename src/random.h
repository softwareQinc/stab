#ifndef STAB_RANDOM_H_
#define STAB_RANDOM_H_

#include <random>

namespace stab {
std::mt19937& get_prng_engine();

int random_bit(double p = 0.5);

// template functions need to be defined in the header

template <class T = int, std::enable_if_t<std::is_integral_v<T>>* = nullptr>
T random_integer(T a = std::numeric_limits<T>::min(),
                 T b = std::numeric_limits<T>::max()) {
    std::uniform_int_distribution<T> dist(a, b);
    if (a > b) {
        throw std::out_of_range("a > b");
    }
    return dist(get_prng_engine());
}

template <class T = double,
          std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
T random_real(T a = std::numeric_limits<T>::min(),
              T b = std::numeric_limits<T>::max()) {
    if (a >= b) {
        throw std::out_of_range("a >= b");
    }
    std::uniform_real_distribution<T> dist(a, b);
    return dist(get_prng_engine());
}
} /* namespace stab */

#endif /* STAB_RANDOM_H_ */
