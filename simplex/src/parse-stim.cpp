#include "parse-stim.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

struct opdatum {
  optype opt;
  std::vector<unsigned> args;
};

struct opdata {
  unsigned arity;
  std::vector<opdatum> expansion;
};

struct instrs parse_file(const char *p) {
  static const std::map<std::string, struct opdata> opmap = {
      {"I", {1, {}}},
      {"X", {1, {{optype::X, {0}}}}},
      {"Y", {1, {{optype::Y, {0}}}}},
      {"Z", {1, {{optype::Z, {0}}}}},
      {"C_XYZ", {1, {
        {optype::Sdg, {0}},
        {optype::H, {0}}}}},
      {"C_ZYX", {1, {
        {optype::H, {0}},
        {optype::S, {0}}}}},
      {"H", {1, {{optype::H, {0}}}}},
      {"H_XY", {1, {
        {optype::H, {0}},
        {optype::Z, {0}},
        {optype::H, {0}},
        {optype::S, {0}}}}},
      {"H_XZ", {1, {{optype::H, {0}}}}},
      {"H_YZ", {1, {
        {optype::H, {0}},
        {optype::S, {0}},
        {optype::H, {0}},
        {optype::Z, {0}}}}},
      {"S", {1, {{optype::S, {0}}}}},
      {"SQRT_X", {1, {
        {optype::H, {0}},
        {optype::S, {0}},
        {optype::H, {0}}}}},
      {"SQRT_X_DAG", {1, {
        {optype::S, {0}},
        {optype::H, {0}},
        {optype::S, {0}}}}},
      {"SQRT_Y", {1, {
        {optype::Z, {0}},
        {optype::H, {0}}}}},
      {"SQRT_Y_DAG", {1, {
        {optype::H, {0}},
        {optype::Z, {0}}}}},
      {"SQRT_Z", {1, {{optype::S, {0}}}}},
      {"SQRT_Z_DAG", {1, {{optype::Sdg, {0}}}}},
      {"S_DAG", {1, {{optype::Sdg, {0}}}}},
      {"CNOT", {2, {{optype::CX, {0, 1}}}}},
      {"CX", {2, {{optype::CX, {0, 1}}}}},
      {"CY", {2, {
        {optype::Sdg, {1}},
        {optype::CX, {0, 1}},
        {optype::S, {1}}}}},
      {"CZ", {2, {{optype::CZ, {0, 1}}}}},
      {"ISWAP", {2, {
        {optype::CX, {0, 1}},
        {optype::S, {1}},
        {optype::CX, {1, 0}},
        {optype::CX, {0, 1}}}}},
      {"ISWAP_DAG", {2, {
        {optype::CX, {0, 1}},
        {optype::Sdg, {1}},
        {optype::CX, {1, 0}},
        {optype::CX, {0, 1}}}}},
      {"SQRT_XX", {2, {
        {optype::CX, {0, 1}},
        {optype::H, {0}},
        {optype::S, {0}},
        {optype::H, {0}},
        {optype::CX, {0, 1}}}}},
      {"SQRT_XX_DAG", {2, {
        {optype::S, {0}},
        {optype::CX, {0, 1}},
        {optype::H, {0}},
        {optype::S, {0}},
        {optype::CX, {0, 1}}}}},
      {"SQRT_YY", {2, {
        {optype::S, {0}},
        {optype::CX, {1, 0}},
        {optype::Z, {0}},
        {optype::H, {1}},
        {optype::CX, {1, 0}},
        {optype::S, {0}}}}},
      {"SQRT_YY_DAG", {2, {
        {optype::CX, {0, 1}},
        {optype::S, {1}},
        {optype::H, {0}},
        {optype::S, {0}},
        {optype::H, {0}},
        {optype::CX, {1, 0}},
        {optype::CX, {0, 1}}}}},
      {"SQRT_ZZ", {2, {
        {optype::CX, {0, 1}},
        {optype::S, {1}},
        {optype::CX, {0, 1}}}}},
      {"SQRT_ZZ_DAG", {2, {
        {optype::H, {1}},
        {optype::CX, {0, 1}},
        {optype::H, {1}},
        {optype::Sdg, {0}},
        {optype::Sdg, {1}}}}},
      {"SWAP", {2, {
        {optype::CX, {0, 1}},
        {optype::CX, {1, 0}},
        {optype::CX, {0, 1}}}}},
      {"XCX", {2, {
        {optype::H, {0}},
        {optype::CX, {0, 1}},
        {optype::H, {0}}}}},
      {"XCY", {2, {
        {optype::CX, {0, 1}},
        {optype::H, {0}},
        {optype::S, {0}},
        {optype::CX, {0, 1}},
        {optype::H, {0}}}}},
      {"XCZ", {2, {{optype::CX, {1, 0}}}}},
      {"YCX", {2, {
        {optype::CX, {0, 1}},
        {optype::H, {1}},
        {optype::S, {1}},
        {optype::CX, {1, 0}},
        {optype::H, {1}}}}},
      {"YCY", {2, {
        {optype::H, {0}},
        {optype::S, {0}},
        {optype::H, {0}},
        {optype::CX, {0, 1}},
        {optype::H, {0}},
        {optype::CX, {1, 0}},
        {optype::S, {0}}}}},
      {"YCZ", {2, {
        {optype::Sdg, {0}},
        {optype::CX, {1, 0}},
        {optype::S, {0}}}}},
      {"ZCX", {2, {{optype::CX, {0, 1}}}}},
      {"ZCY", {2, {
        {optype::Sdg, {1}},
        {optype::CX, {0, 1}},
        {optype::S, {1}}}}},
      {"ZCZ", {2, {{optype::CZ, {0, 1}}}}},
      {"M", {1, {{optype::MeasZ, {0}}}}},
      {"MX", {1, {{optype::MeasX, {0}}}}},
      {"MY", {1, {{optype::MeasY, {0}}}}},
      {"MZ", {1, {{optype::MeasZ, {0}}}}},
      {"R", {1, {{optype::ResetZ, {0}}}}},
      {"RX", {1, {{optype::ResetX, {0}}}}},
      {"RY", {1, {{optype::ResetY, {0}}}}},
      {"RZ", {1, {{optype::ResetZ, {0}}}}},
      {"tick", {0, {}}}
  };

  std::ifstream file(p);
  unsigned max_n = 0;
  std::vector<struct op> ops;
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::vector<std::string> tokens(
      std::istream_iterator<std::string>(iss),
      std::istream_iterator<std::string>{});
    if (tokens.empty()) continue;
    const std::string &opname = tokens[0];
    struct opdata opda = opmap.at(opname);
    unsigned n_args = opda.arity;
    if (tokens.size() != 1 + n_args) {
      std::cerr << "Cannot parse line: " << line << std::endl;
      throw;
    }
    std::vector<unsigned> qubits(n_args);
    for (unsigned i = 0; i < n_args; i++) {
      unsigned k = std::stoul(tokens[1 + i]);
      if (k > max_n) max_n = k;
      qubits[i] = k;
    }
    for (const struct opdatum &opdm : opda.expansion) {
      std::vector<unsigned> qbs(n_args);
      for (unsigned i = 0; i < n_args; i++) {
        qbs[i] = qubits[opdm.args[i]];
      }
      ops.push_back({opdm.opt, qbs});
    }
  }
  return {max_n + 1, ops};
}
