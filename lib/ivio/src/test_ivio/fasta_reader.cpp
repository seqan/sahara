// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include <catch2/catch.hpp>
#include <filesystem>
#include <fstream>
#include <ivio/ivio.h>

TEST_CASE("reading fasta files", "[fasta][reader]") {
    auto tmp = std::filesystem::temp_directory_path() / "ivio_test";
    std::filesystem::create_directory(tmp);

    auto expected = std::vector<ivio::fasta::record> {
        ivio::fasta::record{.id = "I am genomic data",    .seq = "ACGTACGTACGTACGT"},
        ivio::fasta::record{.id = "Some other data",      .seq = "GNNNGNJKL"},
        ivio::fasta::record{.id = "Weird multiline data", .seq = "AGgagatAGagagAGAJKNNNN"},
    };

    auto test_data = std::string {
        ">I am genomic data\n"
        "ACGTACGTACGTACGT\n"
        ">Some other data\n"
        "GNNNGNJKL\n"
        ">Weird multiline data\n"
        "AGgagat\n"
        "AGagagAGA\n"
        "JKNNNN\n"
    };

    SECTION("Creating Test file") {
        auto ofs = std::ofstream{tmp / "file.fa", std::ios::binary};
        ofs << test_data;
    }

    SECTION("Read from std::filesystem::path") {
        auto reader = ivio::fasta::reader{{tmp / "file.fa"}};
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(expected == vec);
    }

    SECTION("Read from std::ifstream") {
        auto fs = std::ifstream{tmp / "file.fa"};
        auto reader = ivio::fasta::reader{{fs}};
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(expected == vec);
    }

    SECTION("Read from std::stringstream") {
        auto ss = std::stringstream{test_data};
        auto reader = ivio::fasta::reader{{ss}};
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(expected == vec);
    }

    SECTION("cleanup - deleting temp folder") {
        std::filesystem::remove_all(tmp);
    }
}
