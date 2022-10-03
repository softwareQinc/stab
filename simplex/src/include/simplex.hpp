#pragma once

#include <iostream>
#include <memory>
#include <optional>

/**
 * Clifford circuit simulator
 */
class Simplex {
public:
  /**
   * Construct a simulator initialized in the all-zero state.
   *
   * @param n number of qubits
   * @param seed seed for PRNG
   */
  Simplex(unsigned n, int seed = 0);

  /**
   * Construct a simulator initialized in the all-zero state and apply the
   * commands specified in a Stim-format file.
   *
   * https://github.com/quantumlib/Stim/blob/main/doc/file_format_stim_circuit.md
   *
   * Not all Stim instruction types are supported.
   *
   * @param p path to Stim file
   * @param seed seed for PRNG
   */
  Simplex(const char *p, int seed = 0);

  ~Simplex();
  Simplex(const Simplex& other);
  Simplex(Simplex&& other);
  Simplex& operator=(const Simplex& other);
  Simplex& operator=(Simplex&& other);

  friend std::ostream& operator<<(std::ostream& os, const Simplex& S);

  /**
   * Get the number of qubits
   *
   * @return number of qubits
   */
  unsigned n() const;

  /**
   * Apply an X gate
   *
   * @param j qubit index
   */
  void X(unsigned j);

  /**
   * Apply a Y gate
   *
   * @param j qubit index
   */
  void Y(unsigned j);

  /**
   * Apply a Z gate
   *
   * @param j qubit index
   */
  void Z(unsigned j);

  /**
   * Apply an H gate
   *
   * @param j qubit index
   */
  void H(unsigned j);

  /**
   * Apply an S gate
   *
   * @param j qubit index
   */
  void S(unsigned j);

  /**
   * Apply an inverse-S gate
   *
   * @param j qubit index
   */
  void Sdg(unsigned j);

  /**
   * Apply a CX gate
   *
   * @param j index of control qubit
   * @param k index of target qubit
   */
  void CX(unsigned j, unsigned k);

  /**
   * Apply a CZ gate
   *
   * @param j index of control qubit
   * @param k index of target qubit
   */
  void CZ(unsigned j, unsigned k);

  /**
   * Measure a qubit in the X basis.
   *
   * If the measurement is non-deterministic, a PRNG is used to determine it,
   * unless a coin is specified.
   *
   * @param j qubit index
   * @param coin value to use (instead of PRNG) if result is non-deterministic
   *
   * @return measurement result
   */
  int MeasX(unsigned j, std::optional<int> coin = std::nullopt);

  /**
   * Measure a qubit in the Y basis
   *
   * If the measurement is non-deterministic, a PRNG is used to determine it,
   * unless a coin is specified.
   *
   * @param j qubit index
   * @param coin value to use (instead of PRNG) if result is non-deterministic
   *
   * @return measurement result
   */
  int MeasY(unsigned j, std::optional<int> coin = std::nullopt);

  /**
   * Measure a qubit in the Z basis
   *
   * If the measurement is non-deterministic, a PRNG is used to determine it,
   * unless a coin is specified.
   *
   * @param j qubit index
   * @param coin value to use (instead of PRNG) if result is non-deterministic
   *
   * @return measurement result
   */
  int MeasZ(unsigned j, std::optional<int> coin = std::nullopt);

  /**
   * Reset a qubit in the X basis by measurement and conditional correction
   *
   * @param j qubit index
   */
  void ResetX(unsigned j);

  /**
   * Reset a qubit in the Y basis by measurement and conditional correction
   *
   * @param j qubit index
   */
  void ResetY(unsigned j);

  /**
   * Reset a qubit in the Z basis by measurement and conditional correction
   *
   * @param j qubit index
   */
  void ResetZ(unsigned j);

  /**
   * Global phase, in units of pi/4
   *
   * @return an integer in the range [0,8) representing the global phase
   */
  int phase() const;

  /**
   * Determine whether all measurements in the circuit are deterministic
   *
   * @retval true all measurements are deterministic
   * @retval false some measurements are non-deterministic
   */
  bool is_deterministic() const;

private:
  struct impl;
  std::unique_ptr<impl> pImpl;
};
