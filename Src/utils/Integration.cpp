#include "utils/Integration.hpp"
#include <boost/math/quadrature/gauss.hpp>
#include <stdexcept>

double Integration::gaussQuadrature1D(const std::function<double(double)>& func, 
                                     double a, double b, int order) {
    
    validateOrder(order);
    
    auto rule = boost::math::quadrature::gauss<double, 0>();
    switch (order) {
        case 1: rule = boost::math::quadrature::gauss<double, 1>(); break;
        case 2: rule = boost::math::quadrature::gauss<double, 2>(); break;
        case 3: rule = boost::math::quadrature::gauss<double, 3>(); break;
        case 4: rule = boost::math::quadrature::gauss<double, 4>(); break;
        case 5: rule = boost::math::quadrature::gauss<double, 5>(); break;
        default: rule = boost::math::quadrature::gauss<double, 5>(); break;
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < rule.size(); ++i) {
        double x = 0.5 * ((b - a) * rule.nodes()[i] + (a + b));
        sum += rule.weights()[i] * func(x);
    }
    
    return 0.5 * (b - a) * sum;
}

double Integration::gaussQuadrature2D(const std::function<double(double, double)>& func,
                                     double a, double b, double c, double d, 
                                     int order) {
    
    validateOrder(order);
    
    auto rule = boost::math::quadrature::gauss<double, 0>();
    switch (order) {
        case 1: rule = boost::math::quadrature::gauss<double, 1>(); break;
        case 2: rule = boost::math::quadrature::gauss<double, 2>(); break;
        case 3: rule = boost::math::quadrature::gauss<double, 3>(); break;
        default: rule = boost::math::quadrature::gauss<double, 3>(); break;
    }
    
    double sum = 0.0;
    for (size_t i = 0; i < rule.size(); ++i) {
        double xi = 0.5 * ((b - a) * rule.nodes()[i] + (a + b));
        for (size_t j = 0; j < rule.size(); ++j) {
            double eta = 0.5 * ((d - c) * rule.nodes()[j] + (c + d));
            sum += rule.weights()[i] * rule.weights()[j] * func(xi, eta);
        }
    }
    
    return 0.25 * (b - a) * (d - c) * sum;
}

std::pair<std::vector<double>, std::vector<double>> 
Integration::getGaussPointsWeights(int order) {
    
    validateOrder(order);
    
    std::vector<double> points, weights;
    
    auto rule = boost::math::quadrature::gauss<double, 0>();
    switch (order) {
        case 1: rule = boost::math::quadrature::gauss<double, 1>(); break;
        case 2: rule = boost::math::quadrature::gauss<double, 2>(); break;
        case 3: rule = boost::math::quadrature::gauss<double, 3>(); break;
        case 4: rule = boost::math::quadrature::gauss<double, 4>(); break;
        case 5: rule = boost::math::quadrature::gauss<double, 5>(); break;
        default: rule = boost::math::quadrature::gauss<double, 5>(); break;
    }
    
    for (size_t i = 0; i < rule.size(); ++i) {
        points.push_back(rule.nodes()[i]);
        weights.push_back(rule.weights()[i]);
    }
    
    return {points, weights};
}

void Integration::validateOrder(int order) {
    if (order < 1 || order > 10) {
        throw std::invalid_argument("Gauss quadrature order must be between 1 and 10");
    }
}