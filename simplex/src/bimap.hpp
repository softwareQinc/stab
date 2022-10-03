#pragma once

#include <iostream>
#include <map>
#include <optional>

/**
 * A bijection between two sets of unsigned integers.
 */
class Bimap {
public:
  Bimap() : fwd(), inv() {}

  std::optional<unsigned> fwd_at(unsigned i) const {
    auto j = fwd.find(i);
    if (j != fwd.end()) {
      return j->second;
    }
    return std::nullopt;
  }

  std::optional<unsigned> inv_at(unsigned j) const {
    auto i = inv.find(j);
    if (i != inv.end()) {
      return i->second;
    }
    return std::nullopt;
  }

  void fwd_erase(unsigned i) {
    auto j = fwd.find(i);
    if (j != fwd.end()) {
      unsigned j1 = j->second;
      fwd.erase(i);
      inv.erase(j1);
    }
  }

  void make_match(unsigned i, unsigned j) {
    std::optional<unsigned> j1 = fwd_at(i);
    std::optional<unsigned> i1 = inv_at(j);
    if (j1 && (*j1 == j)) return;
    insert_pair(i, j);
    if (i1) fwd.erase(*i1);
    if (j1) inv.erase(*j1);
  }

  // Swap fwd[i1] and fwd[i2]
  void swap_fwd(unsigned i1, unsigned i2) {
    std::optional<unsigned> j1 = fwd_at(i1);
    std::optional<unsigned> j2 = fwd_at(i2);
    if (j1) {
      fwd.erase(i1);
      inv.erase(*j1);
    }
    if (j2) {
      fwd.erase(i2);
      inv.erase(*j2);
    }
    if (j1) {
      insert_pair(i2, *j1);
    }
    if (j2) {
      insert_pair(i1, *j2);
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const Bimap& p);

private:
  void insert_pair(unsigned i, unsigned j) {
    fwd[i] = j;
    inv[j] = i;
  }

  std::map<unsigned, unsigned> fwd;
  std::map<unsigned, unsigned> inv;
};

inline std::ostream& operator<<(std::ostream& os, const Bimap& p) {
  for (auto ij : p.fwd) {
    unsigned i = ij.first;
    unsigned j = ij.second;
    os << i << ":" << j << " ";
  }
  return os;
}
