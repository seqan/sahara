#include <ivsigma/ivsigma.h>
#include <iostream>

int main()
{
    //                                 G  G  C  T  A
    auto input = std::vector<uint8_t>{ 2, 2, 1, 3, 0};
    { // Version 1
        auto output = std::string{};
        output.resize(input.size());
        ivs::convert_rank_to_char<ivs::dna4>(input, output);
        std::cout << output << '\n';
    }
    { // Version 2
        std::string output = ivs::convert_rank_to_char<ivs::dna4>(input);
        std::cout << output << '\n';
    }
    { // Version 3
        auto output_view = input | ivs::view_rank_to_char<ivs::dna4>;
        auto output = std::string(output_view.begin(), output_view.end());
        std::cout << output << '\n';
    }
}
