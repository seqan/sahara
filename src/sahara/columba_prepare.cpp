// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <clice/clice.h>
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <libsais64.h>

namespace {
void app();
auto cli = clice::Argument {
    .args   = "columba_prepare",
    .desc   = "takes a fasta file and prepares it for columba",
    .cb     = app,
};

auto cliInput = clice::Argument {
    .parent = &cli,
    .args   = {"-i", "--input"},
    .desc   = "path to a fasta file",
    .value  = std::filesystem::path{},
    .tags   = {"required"},
};

auto cliOutput = clice::Argument {
    .parent = &cli,
    .args   = {"-o", "--output"},
    .desc   = "base path (without extensions)",
    .value  = std::string{},
    .tags   = {"required"},
};

char randomPick() {
    switch(rand()% 4) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
    }
    throw std::runtime_error("should never happen");
}

auto loadFastaAsSingleText(std::filesystem::path input) -> std::vector<uint8_t> {
    // read all queries into one giant text (Columba can not handle mutlistrings
    auto text = std::vector<uint8_t>{};
    for (auto record : ivio::fasta::reader {{input}}) {
        text.reserve(text.size() + record.seq.size()+1);
        for (auto c : record.seq) {
            c = ivs::dna4::normalize_char(c);
            if (!ivs::verify_char(c)) {
                c = randomPick();
            }
            text.emplace_back(c);
        }
    }
    text.emplace_back('$');
    return text;
}

auto createSA(std::span<uint8_t const> text) -> std::vector<int64_t> {
    auto sa = std::vector<int64_t>{};
    sa.resize(text.size());
    int e = libsais64(text.data(), sa.data(), sa.size(), 0, nullptr);
    if (e != 0) {
        throw std::runtime_error{"error creating suffix array with libsais64 (" + std::to_string(e) + ")"};
    }
    return sa;
}

void writeText(std::filesystem::path output, std::span<uint8_t const> text) {
    auto ofs = std::ofstream{output, std::ios::binary};
    ofs.write((char const*)text.data(), text.size());
}


void writeSA(std::filesystem::path output, std::span<int64_t const> sa) {
    auto ofs = std::ofstream{output};

    if (!sa.empty()) {
        ofs << sa[0];
        for (size_t i{1}; i < sa.size(); ++i) {
            ofs << ' ' << sa[i];
        }
    }
}

void app() {
    fmt::print("reading string T from fasta file...\n");
    auto text = loadFastaAsSingleText(*cliInput);
    {
        fmt::print("saving text T to disk...\n");
        writeText(*cliOutput + ".txt", text);
        fmt::print("-> {}\n", *cliOutput + ".txt");

        fmt::print("constructing Suffix Array for T...\n");
        auto sa   = createSA(text);

        fmt::print("saving Suffix Array disk...\n");
        writeSA(*cliOutput + ".sa", sa);
        fmt::print("-> {}\n", *cliOutput + ".sa");
    }
    {
        fmt::print("reversing text T...\n");
        std::ranges::reverse(text);

        fmt::print("saving reversed text T to disk...\n");
        writeText(*cliOutput + ".rev.txt", text);
        fmt::print("-> {}\n", *cliOutput + ".rev.txt");

        fmt::print("constructing Suffix Array for reverse T...\n");
        auto sa   = createSA(text);

        fmt::print("saving Suffix Array (reversed T) disk...\n");
        writeSA(*cliOutput + ".rev.sa", sa);
        fmt::print("-> {}\n", *cliOutput + ".rev.sa");
    }
}
}
