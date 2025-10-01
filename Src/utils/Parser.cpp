#include "Parser.hpp"
#include "exprtk.hpp"

class FunctionParser {
private:
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;
    
    symbol_table_t symbol_table_;
    expression_t expression_;
    parser_t parser_;
    
    double x_, y_, u_;
    
public:
    FunctionParser() {
        symbol_table_.add_variable("x", x_);
        symbol_table_.add_variable("y", y_);
        symbol_table_.add_variable("u", u_);
        symbol_table_.add_constants();
        
        expression_.register_symbol_table(symbol_table_);
    }
    
    bool compile(const std::string& expression_str) {
        return parser_.compile(expression_str, expression_);
    }
    
    double evaluate(double x, double y, double u = 0.0) {
        x_ = x;
        y_ = y;
        u_ = u;
        return expression_.value();
    }
};

std::function<double(double, double)> Parser::parse2D(const std::string& expr) {
    auto parser = std::make_shared<FunctionParser>();
    if (!parser->compile(expr)) {
        throw std::runtime_error("Failed to parse expression: " + expr);
    }
    
    return [parser](double x, double y) {
        return parser->evaluate(x, y);
    };
}

std::function<double(double)> Parser::parse1D(const std::string& expr) {
    auto parser = std::make_shared<FunctionParser>();
    if (!parser->compile(expr)) {
        throw std::runtime_error("Failed to parse expression: " + expr);
    }
    
    return [parser](double u) {
        return parser->evaluate(0, 0, u);
    };
}