#include <ivsigma/ivsigma.h>
#include <iostream>

int main()
{
    std::vector<uint8_t> input = ivs::convert_char_to_rank<ivs::dna4>("GCGACGTAC");
    for (auto enc : ivs::winnowing_minimizer<ivs::dna4>(input, /*._k=*/ 3, /*._window=*/ 2)) {
        std::cout << enc << ' ';
    }
    std::cout << '\n';
}
