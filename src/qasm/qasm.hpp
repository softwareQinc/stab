#ifndef QASM_QASM_HPP_
#define QASM_QASM_HPP_

#include "qasmtools/ast/ast.hpp"
#include "qasmtools/parser/parser.hpp"
#include "qasmtools/utils/angle.hpp"

namespace stab::qasm {

namespace ast = qasmtools::ast;
namespace parser = qasmtools::parser;

class QASMSimulator final : public ast::Visitor {
    stab::AffineState psi;
    double
        temp_value; ///< stores intermediate values when computing expressions

  public:
    explicit QASMSimulator(int num_qubits) : psi(num_qubits), temp_value(0) {}

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

    void visit(ast::PiExpr&) override { temp_value = qasmtools::utils::pi; }

    void visit(ast::IntExpr& expr) override {
        temp_value = static_cast<double>(expr.value());
    }

    void visit(ast::RealExpr& expr) override { temp_value = expr.value(); }

    void visit(ast::VarExpr& expr) override {
        throw std::logic_error("Not supported");
    }

    // Statements
    void visit(ast::MeasureStmt& stmt) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::ResetStmt& stmt) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::IfStmt& stmt) override {
        throw std::logic_error("Not supported");
    }

    // Gates
    void visit(ast::UGate& gate) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::CNOTGate& gate) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::BarrierGate&) override {}

    void visit(ast::DeclaredGate& dgate) override {
        throw std::logic_error("Not supported");
    }

    // Declarations
    void visit(ast::GateDecl& decl) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::OracleDecl& decl) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::RegisterDecl& decl) override {
        throw std::logic_error("Not supported");
    }

    void visit(ast::AncillaDecl& decl) override {
        throw std::logic_error("Not supported");
    }

    // Program
    void visit(ast::Program& prog) override {
        prog.foreach_stmt([this](auto& stmt) { stmt.accept(*this); });
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

} /* namespace stab::qasm */

#endif /* QASM_QASM_HPP_ */
