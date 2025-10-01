#pragma once
#include <string>
#include <functional>
#include <memory>

class Parser {
public:
    static std::function<double(double, double)> parse2D(const std::string& expr);
    static std::function<double(double)> parse1D(const std::string& expr);
    
private:
    class FunctionParserImpl;
};