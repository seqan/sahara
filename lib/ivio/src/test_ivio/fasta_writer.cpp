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

static auto read_file(std::filesystem::path p) -> std::string {
    auto fs = std::ifstream{p};
    auto buffer = std::stringstream{};
    buffer << fs.rdbuf();
    return buffer.str();
}

TEST_CASE("writing fasta files", "[fasta][writer]") {
    auto tmp = std::filesystem::temp_directory_path() / "ivio_test";
    std::filesystem::create_directory(tmp);

    auto test_data = std::vector<ivio::fasta::record> {
        ivio::fasta::record{.id = "I am genomic data",    .seq = "ACGTACGTACGTACGT"},
        ivio::fasta::record{.id = "Some other data",      .seq = "GNNNGNJKL"},
        ivio::fasta::record{.id = "Weird multiline data", .seq = "AGgagatAGagagAGAJKNNNN"},
    };

    auto expected = std::string {
        ">I am genomic data\n"
        "ACGTACGTACGTACGT\n"
        ">Some other data\n"
        "GNNNGNJKL\n"
        ">Weird multiline data\n"
        "AGgagatAGagagAGAJKNNNN\n"
    };

    SECTION("Write to std::filesystem::path") {
        auto writer = ivio::fasta::writer{{tmp / "file.fa"}};
        for (auto r : test_data) {
            writer.write(r);
        }
        writer.close();
        CHECK(read_file(tmp / "file.fa") == expected);
    }

    SECTION("Write to std::ofstream") {
        auto fs = std::ofstream{tmp / "file.fa"};
        auto writer = ivio::fasta::writer{{fs}};
        for (auto r : test_data) {
            writer.write(r);
        }
        writer.close();
        fs.close();
        CHECK(read_file(tmp / "file.fa") == expected);
    }

    SECTION("Write to std::stringstream") {
        auto ss = std::stringstream{};
        auto writer = ivio::fasta::writer{{ss}};
        for (auto r : test_data) {
            writer.write(r);
        }
        writer.close();
        CHECK(ss.str() == expected);
    }

    SECTION("cleanup - deleting temp folder") {
        std::filesystem::remove_all(tmp);
    }
}
