// main.cpp
#include "VariantCaller.hpp"
#include <iostream>

int main(int argc, char **argv) {
    VariantCaller vc;
    if (!vc.parse_arguments(argc, argv)) {
        vc.print_usage();
        return 1;
    }

    if (vc.help) {
        vc.print_usage();
        return 1;
    }

    if (!vc.run()) {
        std::cerr << "Error running variant caller" << std::endl;
        return 1;
    }

    return 0;
}