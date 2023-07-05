// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <iostream>
#include <ivsigma/ivsigma.h>
#include <string>

template <ivs::alphabet_c Alphabet>
static void check_normalize(std::string input, std::string expected) {
    {
        auto output = ivs::normalize_char<Alphabet>(input);
        assert(expected == output);
    }
    {
        auto output = std::string{};
        output.resize(input.size());
        ivs::normalize_char<Alphabet>(input, output);
        assert(expected == output);
    }
    {
        auto output_view = input | ivs::view_normalize_char<Alphabet>;
        assert(std::ranges::equal(expected, output_view));
    }

    // Invalid symbols are converted to '\0'
    {
        auto invalid_input = std::string{};
        for (int c{std::numeric_limits<char>::min()}; c < std::numeric_limits<char>::max(); ++c) {
            if (std::string::npos == input.find(c)) {
                invalid_input += c;
            }
        }
        auto expected = std::string{};
        expected.resize(invalid_input.size(), '\0');
        auto output = ivs::normalize_char<Alphabet>(invalid_input);
        assert(expected == output);
    }
};

template <ivs::alphabet_c Alphabet>
static void check_char_to_rank(std::string input, std::vector<uint8_t> expected) {
    {
        auto output = ivs::convert_char_to_rank<Alphabet>(input);
        assert(expected == output);
    }
    {
        auto output = std::vector<uint8_t>{};
        output.resize(input.size());
        ivs::convert_char_to_rank<Alphabet>(input, output);
        assert(expected == output);
    }
    {
        auto output_view = input | ivs::view_char_to_rank<Alphabet>;
        assert(std::ranges::equal(expected, output_view));
    }

    // Invalid symbols are converted to 255
    {
        auto invalid_input = std::string{};
        for (int c{std::numeric_limits<char>::min()}; c < std::numeric_limits<char>::max(); ++c) {
            if (std::string::npos == input.find(c)) {
                invalid_input += c;
            }
        }
        auto expected = std::vector<uint8_t>{};
        expected.resize(invalid_input.size(), 255);
        auto output = ivs::convert_char_to_rank<Alphabet>(invalid_input);
        assert(expected == output);
    }
};

template <ivs::alphabet_c Alphabet>
static void check_rank_to_char(std::vector<uint8_t> input, std::string expected) {
    {
        auto output = ivs::convert_rank_to_char<Alphabet>(input);
        assert(expected == output);
    }
    {
        auto output = std::string{};
        output.resize(input.size());
        ivs::convert_rank_to_char<Alphabet>(input, output);
        assert(expected == output);
    }
    {
        auto output_view = input | ivs::view_rank_to_char<Alphabet>;
        assert(std::ranges::equal(expected, output_view));
    }

    // Invalid ranks are converted to '\0'
    {
        auto invalid_input = std::vector<uint8_t>{};
        for (int r{0}; r < 256; ++r) {
            if (input.end() == std::ranges::find(input, r)) {
                invalid_input.push_back(r);;
            }
        }
        auto expected = std::string{};
        expected.resize(invalid_input.size(), '\0');
        auto output = ivs::convert_rank_to_char<Alphabet>(invalid_input);
        assert(expected == output);
    }
};

template <ivs::alphabet_c Alphabet>
static void check_complement_char(std::string input, std::string expected) {
    {
        auto output = ivs::complement_char<Alphabet>(input);
        assert(expected == output);
    }
    {
        auto output = std::string{};
        output.resize(input.size());
        ivs::complement_char<Alphabet>(input, output);
        assert(expected == output);
    }
    {
        auto output_view = input | ivs::view_complement_char<Alphabet>;
        assert(std::ranges::equal(expected, output_view));
    }

    // Invalid symbols are converted to '\0'
    {
        auto invalid_input = std::string{};
        for (int c{std::numeric_limits<char>::min()}; c < std::numeric_limits<char>::max(); ++c) {
            if (std::string::npos == input.find(c)) {
                invalid_input += c;
            }
        }
        auto expected = std::string{};
        expected.resize(invalid_input.size(), '\0');
        auto output = ivs::complement_char<Alphabet>(invalid_input);
        assert(expected == output);
    }
};

template <ivs::alphabet_c Alphabet>
static void check_complement_rank(std::vector<uint8_t> input, std::vector<uint8_t> expected) {
    {
        auto output = ivs::complement_rank<Alphabet>(input);
        assert(expected == output);
    }
    {
        auto output = std::vector<uint8_t>{};
        output.resize(input.size());
        ivs::complement_rank<Alphabet>(input, output);
        assert(expected == output);
    }
    {
        auto output_view = input | ivs::view_complement_rank<Alphabet>;
        assert(std::ranges::equal(expected, output_view));
    }

    // Invalid ranks are converted to 255
    {
        auto invalid_input = std::vector<uint8_t>{};
        for (int r{0}; r < 256; ++r) {
            if (input.end() == std::ranges::find(input, r)) {
                invalid_input.push_back(r);;
            }
        }
        auto expected = std::vector<uint8_t>{};
        expected.resize(invalid_input.size(), 255);
        auto output = ivs::complement_rank<Alphabet>(invalid_input);
        assert(expected == output);
    }
};

template <ivs::alphabet_c Alphabet>
static void check_reverse_complement_char(std::string input, std::string expected) {
    {
        auto output = ivs::reverse_complement_char<Alphabet>(input);
        assert(expected == output);
    }
    {
        auto output = std::string{};
        output.resize(input.size());
        ivs::reverse_complement_char<Alphabet>(input, output);
        assert(expected == output);
    }
    {
        auto output_view = input | ivs::view_reverse_complement_char<Alphabet>;
        assert(std::ranges::equal(expected, output_view));
    }

    // Invalid symbols are converted to '\0'
    {
        auto invalid_input = std::string{};
        for (int c{std::numeric_limits<char>::min()}; c < std::numeric_limits<char>::max(); ++c) {
            if (std::string::npos == input.find(c)) {
                invalid_input += c;
            }
        }
        auto expected = std::string{};
        expected.resize(invalid_input.size(), '\0');
        auto output = ivs::reverse_complement_char<Alphabet>(invalid_input);
        assert(expected == output);
    }
};

template <ivs::alphabet_c Alphabet>
static void check_reverse_complement_rank(std::vector<uint8_t> input, std::vector<uint8_t> expected) {
    {
        auto output = ivs::reverse_complement_rank<Alphabet>(input);
        assert(expected == output);
    }
    {
        auto output = std::vector<uint8_t>{};
        output.resize(input.size());
        ivs::reverse_complement_rank<Alphabet>(input, output);
        assert(expected == output);
    }
    {
        auto output_view = input | ivs::view_reverse_complement_rank<Alphabet>;
        assert(std::ranges::equal(expected, output_view));
    }

    // Invalid ranks are converted to 255
    {
        auto invalid_input = std::vector<uint8_t>{};
        for (int r{0}; r < 256; ++r) {
            if (input.end() == std::ranges::find(input, r)) {
                invalid_input.push_back(r);;
            }
        }
        auto expected = std::vector<uint8_t>{};
        expected.resize(invalid_input.size(), 255);
        auto output = ivs::reverse_complement_rank<Alphabet>(invalid_input);
        assert(expected == output);
    }
};



void test_nucliotides() {
    check_normalize<ivs::dna4>("AaCcGgTtUu", "AACCGGTTTT");
    check_normalize<ivs::dna5>("AaCcGgTtUuNn", "AACCGGTTTTNN");
    check_normalize<ivs::rna4>("AaCcGgUuTt", "AACCGGUUUU");
    check_normalize<ivs::rna5>("AaCcGgUuTtNn", "AACCGGUUUUNN");
    check_normalize<ivs::iupac>("AaCcGgTtNnRrYySsWwKkMmBbDdHhVv-.Uu", "AACCGGTTNNRRYYSSWWKKMMBBDDHHVV--TT");
    check_normalize<ivs::dna3bs>("AaCcGgTtUu", "AATTGGTTTT");

    check_char_to_rank<ivs::dna4>("AaCcGgTtUu", {0, 0, 1, 1, 2, 2, 3, 3, 3, 3});
    check_char_to_rank<ivs::dna5>("AaCcGgTtUuNn", {0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4});
    check_char_to_rank<ivs::rna4>("AaCcGgUuTt", {0, 0, 1, 1, 2, 2, 3, 3, 3, 3});
    check_char_to_rank<ivs::rna5>("AaCcGgUuTtNn", {0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4});
    check_char_to_rank<ivs::iupac>("AaCcGgTtNnRrYySsWwKkMmBbDdHhVv-.Uu", {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 3, 3});
    check_char_to_rank<ivs::dna3bs>("AaGgTtCcUu", {0, 0, 1, 1, 2, 2, 2, 2, 2, 2});


    check_rank_to_char<ivs::dna4>( {0, 1, 2, 3}, "ACGT");
    check_rank_to_char<ivs::dna5>( {0, 1, 2, 3, 4}, "ACGTN");
    check_rank_to_char<ivs::rna4>( {0, 1, 2, 3}, "ACGU");
    check_rank_to_char<ivs::rna5>( {0, 1, 2, 3, 4}, "ACGUN");
    check_rank_to_char<ivs::iupac>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, "ACGTNRYSWKMBDHV-");
    check_rank_to_char<ivs::dna3bs>( {0, 1, 2}, "AGT");

    check_complement_char<ivs::dna4>("AaCcGgTtUu", "TTGGCCAAAA");
    check_complement_char<ivs::dna5>("AaCcGgTtUuNn", "TTGGCCAAAANN");
    check_complement_char<ivs::rna4>("AaCcGgUuTt", "UUGGCCAAAA");
    check_complement_char<ivs::rna5>("AaCcGgUuTtNn", "UUGGCCAAAANN");
    check_complement_char<ivs::iupac>("AaCcGgTtNnRrYySsWwKkMmBbDdHhVv-.Uu", "TTGGCCAANNYYRRSSWWMMKKVVHHDDBB--AA");

    check_complement_rank<ivs::dna4>({0, 1, 2, 3}, {3, 2, 1, 0});
    check_complement_rank<ivs::dna5>({0, 1, 2, 3, 4}, {3, 2, 1, 0, 4});
    check_complement_rank<ivs::rna4>({0, 1, 2, 3}, {3, 2, 1, 0});
    check_complement_rank<ivs::rna5>({0, 1, 2, 3, 4}, {3, 2, 1, 0, 4});
    check_complement_rank<ivs::iupac>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, {3, 2, 1, 0, 4, 6, 5, 7, 8, 10, 9, 14, 13, 12, 11, 15});

    check_reverse_complement_char<ivs::dna4>("AaCcGgTtUu", "AAAACCGGTT");
    check_reverse_complement_char<ivs::dna5>("AaCcGgTtUuNn", "NNAAAACCGGTT");
    check_reverse_complement_char<ivs::rna4>("AaCcGgUuTt", "AAAACCGGUU");
    check_reverse_complement_char<ivs::rna5>("AaCcGgUuTtNn", "NNAAAACCGGUU");
    check_reverse_complement_char<ivs::iupac>("AaCcGgTtNnRrYySsWwKkMmBbDdHhVv-.Uu", "AA--BBDDHHVVKKMMWWSSRRYYNNAACCGGTT");

    check_reverse_complement_rank<ivs::dna4>({0, 1, 2, 3}, {0, 1, 2, 3});
    check_reverse_complement_rank<ivs::dna5>({0, 1, 2, 3, 4}, {4, 0, 1, 2, 3});
    check_reverse_complement_rank<ivs::rna4>({0, 1, 2, 3}, {0, 1, 2, 3});
    check_reverse_complement_rank<ivs::rna5>({0, 1, 2, 3, 4}, {4, 0, 1, 2, 3});
    check_reverse_complement_rank<ivs::iupac>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, {15, 11, 12, 13, 14, 9, 10, 8, 7, 5, 6, 4, 0, 1, 2, 3});
}

void test_aminoacids() {

    check_normalize<ivs::aa27>("AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYyUuBbZzJjOoXx*",
                               "AACCDDEEFFGGHHIIKKLLMMNNPPQQRRSSTTVVWWYYUUBBZZJJOOXX*");
    check_normalize<ivs::aa20>("AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYyUuBbZzJjOoXx*",
                               "AACCDDEEFFGGHHIIKKLLMMNNPPQQRRSSTTVVWWYYCCDDEELLLLSSW");
    check_normalize<ivs::aa10li>("AaSsTtXxBbDdEeQqZzCcUuFfWwYy*GgHhNnIiVvJjLlMmKkOoRrPp",
                                 "AAAAAAAABBBBBBBBBBCCCCFFFFFFFGGHHHHIIIIJJJJJJKKKKKKPP");
    check_normalize<ivs::aa10murphy>("AaBbDdEeNnQqZzCcUuFfWwYy*GgHhIiJjLlMmVvKkOoRrPpSsTtXx",
                                     "AABBBBBBBBBBBBCCCCFFFFFFFGGHHIIIIIIIIIIKKKKKKPPSSSSSS");


    check_char_to_rank<ivs::aa27>("AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYyUuBbZzJjOoXx*",
                                  {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,
                                   10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17,
                                   18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25,
                                   26});
    check_char_to_rank<ivs::aa20>("AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYyUuBbZzJjOoXx*",
                                  {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,
                                   10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17,
                                   18, 18, 19, 19,  1,  1,  2,  2,  3,  3,  9,  9,  9,  9, 15, 15,
                                   18});
    check_char_to_rank<ivs::aa10li>("AaSsTtXxBbDdEeQqZzCcUuFfWwYy*GgHhNnIiVvJjLlMmKkOoRrPp",
                                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                                     3, 3, 3, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7,
                                     7, 8, 8, 8, 8, 8, 8, 9, 9});
    check_char_to_rank<ivs::aa10murphy>("AaBbDdEeNnQqZzCcUuFfWwYy*GgHhIiJjLlMmVvKkOoRrPpSsTtXx",
                                        {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3,
                                         3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7,
                                         7, 7, 7, 8, 8, 9, 9, 9, 9, 9, 9});

    check_rank_to_char<ivs::aa27>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                                   18, 19, 20, 21, 22, 23, 24, 25, 26},
                                   "ACDEFGHIKLMNPQRSTVWYUBZJOX*");
    check_rank_to_char<ivs::aa20>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                                   18, 19},
                                   "ACDEFGHIKLMNPQRSTVWY");
    check_rank_to_char<ivs::aa10li>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                                    "ABCFGHIJKP");
    check_rank_to_char<ivs::aa10murphy>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                                        "ABCFGHIKPS");

}
void test_qualities() {
    //!TODO
}

void test_compact_encoding() {
    auto v = std::vector<uint8_t>{3, 2, 0, 1, 4, 3, 2};

    {
        // attention: dna5 has a reverse complement, so the canonical kmers are being taken
        auto result = std::vector<size_t>{};
        for (auto h : ivs::compact_encoding<ivs::dna5>{v, /*.k=*/ 3}) {
            result.push_back(h);
        }

        assert(result.size() == 5);
        assert(result[0] ==  80);
        assert(result[1] ==  51);
        assert(result[2] ==   9);
        assert(result[3] ==  22);
        assert(result[4] ==  29);
    }

    {
        auto result = std::vector<size_t>{};
        for (auto h : ivs::compact_encoding<ivs::aa10li>{v, /*.k=*/ 3}) {
            result.push_back(h);
        }

        assert(result.size() == 5);
        assert(result[0] == 320);
        assert(result[1] == 201);
        assert(result[2] ==  14);
        assert(result[3] == 143);
        assert(result[4] == 432);
    }
}

void test_winnowing_minimizer() {
    auto v = std::vector<uint8_t>{3, 2, 0, 1, 4, 3, 2};

    {
        // attention: dna5 has a reverse complement, so the canonical kmers are being taken
        auto result = std::vector<size_t>{};
        for (auto h : ivs::winnowing_minimizer<ivs::dna5>{v, /*.k=*/ 3, /*.window=*/ 1}) {
            result.push_back(h);
        }

        assert(result.size() == 5);
        assert(result[0] ==  80);
        assert(result[1] ==  51);
        assert(result[2] ==   9);
        assert(result[3] ==  22);
        assert(result[4] ==  29);
    }

    {
        auto result = std::vector<size_t>{};
        for (auto h : ivs::winnowing_minimizer<ivs::aa10li>{v, /*.k=*/ 3, /*.window=*/ 1}) {
            result.push_back(h);
        }

        assert(result.size() == 5);
        assert(result[0] == 320);
        assert(result[1] == 201);
        assert(result[2] ==  14);
        assert(result[3] == 143);
        assert(result[4] == 432);
    }

    {
        // attention: dna5 has a reverse complement, so the canonical kmers are being taken
        auto result = std::vector<size_t>{};
        for (auto h : ivs::winnowing_minimizer<ivs::dna5>{v, /*.k=*/ 3, /*.window=*/ 2}) {
            result.push_back(h);
        }

        assert(result.size() == 3);
        assert(result[0] ==  51);
        assert(result[1] ==   9);
        assert(result[2] ==  22);
    }

    {
        auto result = std::vector<size_t>{};
        for (auto h : ivs::winnowing_minimizer<ivs::aa10li>{v, /*.k=*/ 3, /*.window=*/ 2}) {
            result.push_back(h);
        }

        assert(result.size() == 3);
        assert(result[0] == 201);
        assert(result[1] ==  14);
        assert(result[2] == 143);
    }

    // Test with many duplicates
    {
        auto v = std::vector<uint8_t>{3, 3, 3, 3, 3, 3};

        auto result = std::vector<size_t>{};
        for (auto h : ivs::winnowing_minimizer<ivs::dna4>{v, /*.k=*/ 3, /*.window=*/ 2}) {
            result.push_back(h);
        }

        assert(result.size() == 3);
        assert(result[0] == 0);
        assert(result[1] == 0);
        assert(result[2] == 0);
    }
    // Test with many duplicates, but skipping them
    {
        auto v = std::vector<uint8_t>{3, 3, 3, 3, 3, 3};

        auto result = std::vector<size_t>{};
        for (auto h : ivs::winnowing_minimizer<ivs::dna4, /*DuplicatesAllowed=*/ false>{v, /*.k=*/ 3, /*.window=*/ 2}) {
            result.push_back(h);
        }

        assert(result.size() == 1);
        assert(result[0] == 0);
    }

}

int main() {
    test_nucliotides();
    test_aminoacids();
    test_qualities();
    test_compact_encoding();
    test_winnowing_minimizer();
    using namespace std::literals;
    assert(!ivs::verify_char("ACGT"s));
    assert(ivs::verify_char("ACG\0T"s).value() == 3);
    assert(!ivs::verify_rank(std::vector<uint8_t>{0, 1, 2, 3, 128, 254}));
    assert(ivs::verify_rank(std::vector<uint8_t>{0, 1, 2, 3, 255, 128, 254}).value() == 4);

    // using convert function to convert string to ranks
    {
        auto input  = std::string{"ACGnACGt"};
        if (ivs::dna4::char_to_rank(input[3]) == 4) {
            input[3] = 'A';
        }
        auto output = ivs::convert_char_to_rank<ivs::dna5>(input);
        auto s = fmt::format("{} => {}", input, output);
        assert(s == "ACGnACGt => [0, 1, 2, 4, 0, 1, 2, 3]");
    }
    // using view to convert string to ranks
    {
        auto input  = std::string{"ACGTACGT"};
        auto output = input | ivs::view_char_to_rank<ivs::dna4>;
        auto s = fmt::format("{} => {}", input, output);
        assert(s == "ACGTACGT => [0, 1, 2, 3, 0, 1, 2, 3]");
    }
}
