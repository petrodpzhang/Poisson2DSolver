#pragma once
#include <functional>
#include <vector>

class Integration {
public:
    static double gaussQuadrature1D(const std::function<double(double)>& func, 
                                   double a, double b, int order = 5);
    
    static double gaussQuadrature2D(const std::function<double(double, double)>& func,
                                   double a, double b, double c, double d, 
                                   int order = 3);
    
    static std::pair<std::vector<double>, std::vector<double>> 
    getGaussPointsWeights(int order);
    
private:
    static void validateOrder(int order);
};