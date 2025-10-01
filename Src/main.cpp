#include "core/Config.hpp"
#include "core/Solver.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file>" << std::endl;
        return 1;
    }
    
    try {
        // read config
        Config config = Config::fromJSON(argv[1]);
        
        // create solver
        NonlinearPoissonSolver solver(config);
        
        // solve
        if (solver.solve()) {
            std::cout << "Solution converged successfully!" << std::endl;
            solver.outputResults(config.output_path);
        } else {
            std::cout << "Solution failed to converge!" << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}