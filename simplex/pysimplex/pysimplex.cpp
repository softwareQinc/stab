#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <simplex.hpp>

#include <optional>
#include <ostream>
#include <sstream>

namespace py = pybind11;

PYBIND11_MODULE(_simplex, m) {
  py::class_<Simplex>(m, "Simplex",
    "Clifford circuit simulator")
    .def(py::init<unsigned>(),
        "Initialize a simulator with `n` qubits.",
        py::arg("n"))
    .def(py::init<char *, int>(),
        "Initialize a simulator from a Stim-format file."
        "\n\n",
        "Accepts a file path (as a string), and optionally an RNG seed.",
        py::arg("p"), py::arg("seed") = 0)
    .def("__repr__",
        [](const Simplex& S) {
            std::stringstream ss;
            ss << "Simplex(" << std::endl << S << ")";
            return ss.str();
        })
    .def("X",
        [](Simplex *S, unsigned j) { S->X(j); return S; },
        "Apply an X gate to qubit `j`.",
        py::arg("j"))
    .def("Y",
        [](Simplex *S, unsigned j) { S->Y(j); return S; },
        "Apply a Y gate to qubit `j`.",
        py::arg("j"))
    .def("Z",
        [](Simplex *S, unsigned j) { S->Z(j); return S; },
        "Apply a Z gate to qubit `j`.",
        py::arg("j"))
    .def("H",
        [](Simplex *S, unsigned j) { S->H(j); return S; },
        "Apply an H gate to qubit `j`.",
        py::arg("j"))
    .def("S",
        [](Simplex *S, unsigned j) { S->S(j); return S; },
        "Apply an S gate to qubit `j`.",
        py::arg("j"))
    .def("Sdg",
        [](Simplex *S, unsigned j) { S->Sdg(j); return S; },
        "Apply an inverse S gate to qubit `j`.",
        py::arg("j"))
    .def("CX",
        [](Simplex *S, unsigned j, unsigned k) { S->CX(j, k); return S; },
        "Apply a CX gate to qubits `j` and `k`.",
        py::arg("j"), py::arg("k"))
    .def("CZ",
        [](Simplex *S, unsigned j, unsigned k) { S->CZ(j, k); return S; },
        "Apply a CZ gate to qubits `j` and `k`.",
        py::arg("j"), py::arg("k"))
    .def("MeasX",
        [](Simplex& S, unsigned j, std::optional<int> coin) {
            return S.MeasX(j, coin);
        },
        "Measure qubit `j` in the X basis."
        "\n\n"
        "If the outcome is non-deterministic, then the outcome is the value of "
        "`coin` if specified (must be 0 or 1), otherwise random.",
        py::arg("j"), py::arg("coin") = std::nullopt)
    .def("MeasY",
        [](Simplex& S, unsigned j, std::optional<int> coin) {
            return S.MeasY(j, coin);
        },
        "Measure qubit `j` in the Y basis."
        "\n\n"
        "If the outcome is non-deterministic, then the outcome is the value of "
        "`coin` if specified (must be 0 or 1), otherwise random.",
        py::arg("j"), py::arg("coin") = std::nullopt)
    .def("MeasZ",
        [](Simplex& S, unsigned j, std::optional<int> coin) {
            return S.MeasZ(j, coin);
        },
        "Measure qubit `j` in the Z basis."
        "\n\n"
        "If the outcome is non-deterministic, then the outcome is the value of "
        "`coin` if specified (must be 0 or 1), otherwise random.",
        py::arg("j"), py::arg("coin") = std::nullopt)
    .def("ResetX",
        [](Simplex *S, unsigned j) { S->ResetX(j); },
        "Reset qubit `j` in the X basis by measurement and conditional "
        "correction.",
        py::arg("j"))
    .def("ResetY",
        [](Simplex *S, unsigned j) { S->ResetY(j); },
        "Reset qubit `j` in the Y basis by measurement and conditional "
        "correction.",
        py::arg("j"))
    .def("ResetZ",
        [](Simplex *S, unsigned j) { S->ResetZ(j); },
        "Reset qubit `j` in the Z basis by measurement and conditional "
        "correction.",
        py::arg("j"))
    .def_property_readonly("phase",
        &Simplex::phase,
        "Global phase, in units of pi/4 (an integer in the range [0,8))")
    .def("is_deterministic",
        &Simplex::is_deterministic,
        "Report whether all measurements are deterministic.");
}
