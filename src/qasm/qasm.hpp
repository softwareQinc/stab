#pragma once

#include "qasmtools_stab/ast/ast.hpp"
#include "qasmtools_stab/parser/parser.hpp"
#include "qasmtools_stab/utils/angle.hpp"

namespace stab::qasm_simulator {

namespace ast = qasmtools_stab::ast;
namespace parser = qasmtools_stab::parser;

class QASMSimulator final : public ast::Visitor {
    stab::AffineState psi;
    double
        temp_value; ///< stores intermediate values when computing expressions
    std::unordered_map<std::string, std::pair<int, int>> qregs;
    int num_qubits;

  public:
    explicit QASMSimulator(int nq) : psi(nq), temp_value(0), num_qubits(0) {}

    stab::AffineState run(ast::Program& prog) {
        prog.accept(*this);
        return psi;
    }

    // Variables
    void visit(ast::VarAccess&) override {}

    // Expressions
    // - Set temp_value to the value of the expression
    void visit(ast::BExpr& expr) override {
        expr.lexp().accept(*this);
        double lval = temp_value;
        expr.rexp().accept(*this);
        double rval = temp_value;

        switch (expr.op()) {
            case ast::BinaryOp::Plus:
                temp_value = lval + rval;
                break;
            case ast::BinaryOp::Minus:
                temp_value = lval - rval;
                break;
            case ast::BinaryOp::Times:
                temp_value = lval * rval;
                break;
            case ast::BinaryOp::Divide:
                temp_value = lval / rval;
                break;
            case ast::BinaryOp::Pow:
                temp_value = pow(lval, rval);
                break;
            default:
                temp_value = 0;
        }
    }

    void visit(ast::UExpr& expr) override {
        expr.subexp().accept(*this);
        double val = temp_value;

        switch (expr.op()) {
            case ast::UnaryOp::Neg:
                temp_value = -val;
                break;
            case ast::UnaryOp::Sin:
                temp_value = std::sin(val);
                break;
            case ast::UnaryOp::Cos:
                temp_value = std::cos(val);
                break;
            case ast::UnaryOp::Tan:
                temp_value = std::tan(val);
                break;
            case ast::UnaryOp::Ln:
                temp_value = std::log(val);
                break;
            case ast::UnaryOp::Sqrt:
                temp_value = std::sqrt(val);
                break;
            case ast::UnaryOp::Exp:
                temp_value = std::exp(val);
                break;
            default:
                temp_value = 0;
        }
    }

    void visit(ast::PiExpr&) override { temp_value = qasmtools_stab::utils::pi; }

    void visit(ast::IntExpr& expr) override {
        temp_value = static_cast<double>(expr.value());
    }

    void visit(ast::RealExpr& expr) override { temp_value = expr.value(); }

    void visit(ast::VarExpr& expr) override {
        throw std::logic_error("Not supported");
    }

    // Statements
    void visit(ast::MeasureStmt& stmt) override {
        std::vector<int> id = get_ids(stmt.q_arg());
        for (int i : id) {
            psi.MeasureZ(i);
        }
    }

    void visit(ast::ResetStmt& stmt) override {
        std::vector<int> id = get_ids(stmt.arg());
        for (int i : id) {
            psi.Reset(i);
        }
    }

    void visit(ast::IfStmt& stmt) override {
        throw std::logic_error("If statements are not supported");
    }

    // Gates
    void visit(ast::UGate& gate) override {
        throw std::logic_error("Custom U-gates are not supported");
    }

    void visit(ast::CNOTGate& gate) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::BarrierGate&) override {}

    void visit(ast::DeclaredGate& dgate) override { // TODO: Optimize?
        // First, get qubit ids.
        std::vector<int> id1 = get_ids(dgate.qarg(0));
        std::vector<int> id2; // Only needed for two-qubit gates:
        if (dgate.name() == "cx" || dgate.name() == "cz" ||
            dgate.name() == "swap") {
            id2 = get_ids(dgate.qarg(1));
        }
        int s1 = id1.size();
        int s2 = id2.size();

        // Now apply the gates.
        // One-qubit gates:
        if (dgate.name() == "h") {
            for (int i : id1) {
                psi.H(i);
            }
        } else if (dgate.name() == "s") {
            for (int i : id1) {
                psi.S(i);
            }
        } else if (dgate.name() == "sdg") {
            for (int i : id1) {
                psi.SDG(i);
            }
        } else if (dgate.name() == "x") {
            for (int i : id1) {
                psi.X(i);
            }
        } else if (dgate.name() == "y") {
            for (int i : id1) {
                psi.Y(i);
            }
        } else if (dgate.name() == "z") {
            for (int i : id1) {
                psi.Z(i);
            }
        } /*Now two-qubit gates*/ else if (dgate.name() == "cx") {
            if (s1 == 1 && s2 == 1) {
                psi.CX(id1[0], id2[0]);
            } else if (s1 == 1) {
                for (int i : id2) {
                    psi.CX(id1[0], i);
                }
            } else if (s2 == 1) {
                for (int i : id1) {
                    psi.CX(i, id2[0]);
                }
            } else {
                for (size_t i = 0; i < s1; ++i) {
                    psi.CX(id1[i], id2[i]);
                }
            }
        } else if (dgate.name() == "cz") {
            if (s1 == 1 && s2 == 1) {
                psi.CZ(id1[0], id2[0]);
            } else if (s1 == 1) {
                for (size_t i = 0; i < s2; ++i) {
                    psi.CZ(id1[0], id2[i]);
                }
            } else if (s2 == 1) {
                for (size_t i = 0; i < s1; ++i) {
                    psi.CZ(id1[i], id2[0]);
                }
            } else {
                for (size_t i = 0; i < s1; ++i) {
                    psi.CZ(id1[i], id2[i]);
                }
            }
        } else if (dgate.name() == "swap") {
            for (size_t i = 0; i < id1.size(); ++i) {
                psi.SWAP(id1[i], id2[i]);
            }
        } /*Otherwise it's some other unsupported gate*/ else {
            throw std::logic_error("Gate must be one of {x, y, z, h, s, sdg, cx, cz, swap}.");
        }
        return;
    }

    // Declarations
    void visit(ast::GateDecl& decl) override {
        //throw std::logic_error("Not supported");
    }

    void visit(ast::OracleDecl& decl) override {
        throw std::logic_error("Oracles are not supported");
    }

    void visit(ast::AncillaDecl& decl) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::RegisterDecl& decl) override {
        if (decl.is_quantum()) {
            qregs[decl.id()] = {num_qubits, decl.size()};
            num_qubits += decl.size();
        }
    }

    // Program
    void visit(ast::Program& prog) override {
        prog.foreach_stmt([this](auto& stmt) { stmt.accept(*this); });
    }

  private:
    std::vector<int> get_ids(ast::VarAccess& va) {
        std::vector<int> ids;
        if (va.offset()) {
            ids.push_back(qregs[va.var()].first + *va.offset());
        } else {
            for (int i = 0; i < qregs[va.var()].second; ++i) {
                ids.push_back(qregs[va.var()].first + i);
            }
        }
        return ids;
    }
};

inline void simulate(std::istream& stream) {
    ast::ptr<ast::Program> program = parser::parse_stream(stream);

    QASMSimulator qs(program->qubits());
    auto psi = qs.run(*program);
    std::cout << psi << std::endl;
}

inline void simulate_file(const std::string& fname) {
    ast::ptr<ast::Program> program = parser::parse_file(fname);

    QASMSimulator qs(program->qubits());
    auto psi = qs.run(*program);
    std::cout << psi << std::endl;
}

} /* namespace stab::qasm_simulator */
