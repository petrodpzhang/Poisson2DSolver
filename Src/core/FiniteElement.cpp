#include "core/FiniteElement.hpp"
#include "elements/TriangleElement.hpp"
#include "elements/QuadElement.hpp"
#include <stdexcept>

namespace FiniteElementFactory {
    std::unique_ptr<Element> createElement(const std::string& type, 
                                          const std::vector<int>& node_indices) {
        if (type == "triangle") {
            return std::make_unique<TriangleElement>(node_indices);
        } else if (type == "quad") {
            return std::make_unique<QuadElement>(node_indices);
        } else {
            throw std::invalid_argument("Unknown element type: " + type);
        }
    }
}