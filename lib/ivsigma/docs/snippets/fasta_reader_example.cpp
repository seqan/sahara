#include <fmt/format.h>
#include <fmt/ranges.h>
#include <ivsigma/ivsigma.h>

int main() {
    auto input  = std::string{"ACGnACGt"};
    auto output = ivs::convert_char_to_rank<ivs::dna5>(input);
    fmt::print("{} => {}\n", input, output);
}
