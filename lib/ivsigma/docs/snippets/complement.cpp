#include <ivsigma/ivsigma.h>
#include <iostream>

int main()
{
    { // complement_rank
        //                                 A  C  G  T
        auto input = std::vector<uint8_t>{ 0, 1, 2, 3};
        { // Version 1
            auto output = std::vector<uint8_t>{};
            output.resize(input.size());
            ivs::complement_rank<ivs::dna4>(input, output);
            for (auto r : output) {
                std::cout << (int)r;
            }
            std::cout << '\n';
        }

        { // Version 2
            std::vector<uint8_t> output = ivs::complement_rank<ivs::dna4>(input);
            for (auto r : output) {
                std::cout << (int)r;
            }
            std::cout << '\n';
        }
        { // Version 3
            auto output_view = ivs::complement_rank<ivs::dna4>(input);
            auto output = std::vector<uint8_t>(output_view.begin(), output_view.end());
            for (auto r : output) {
                std::cout << (int)r;
            }
            std::cout << '\n';
        }
    }

    { // complement_char
        auto input = std::string{"AaCcGgTtUu"};
        { // Version 4
            auto output = std::string{};
            output.resize(input.size());
            ivs::complement_char<ivs::dna4>(input, output);
            std::cout << output << '\n';
        }
        { // Version 5
            std::string output = ivs::complement_char<ivs::dna4>(input);
            std::cout << output << '\n';
        }
        { // Version 6
            auto output_view = input | ivs::view_complement_char<ivs::dna4>;
            auto output = std::string(output_view.begin(), output_view.end());
            std::cout << output << '\n';
        }
    }
}
