#pragma once

#include <vector>

enum optype {
  X,
  Y,
  Z,
  H,
  S,
  Sdg,
  CX,
  CZ,
  MeasX,
  MeasY,
  MeasZ,
  ResetX,
  ResetY,
  ResetZ
};

struct op {
  optype type;
  std::vector<unsigned> qubits;
};

struct instrs {
  unsigned n; // number of qubits
  std::vector<struct op> ops;
};

/**
 * Parse a Stim file.
 *
 * Only a subset of Stim syntax is supported.
 *
 * @param p path to file
 *
 * @return parsed list of instructions
 */
struct instrs parse_file(const char *p);
