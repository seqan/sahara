#include <ivsigma/ivsigma.h>
#include <iostream>

int main() {
    auto ranks = ivs::convert_char_to_rank<ivs::dna4>("acgt");
    for (uint8_t c : ranks) {
        std::cout << (int)c << "\n";
    }
}
