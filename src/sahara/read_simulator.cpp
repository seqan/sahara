// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "error_fmt.h"

#include <clice/clice.h>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <random>

namespace {
void app();
auto cli = clice::Argument{ .args   = "read_simulator",
                            .desc   = "simulates reads of a certain length",
                            .cb     = app,
};

auto cliInput = clice::Argument{ .parent = &cli,
                                 .args   = {"-i", "--input"},
                                 .desc   = "path to a fasta file",
                                 .value  = std::filesystem::path{},
};

auto cliOutput = clice::Argument{ .parent = &cli,
                                  .args   = {"-o", "--output"},
                                  .desc   = "path to the output fasta file",
                                  .value  = std::filesystem::path{},
                                  .tags   = {"required"},
};

auto cliFastaLineLength = clice::Argument{ .parent = &cli,
                                           .args   = {"--fasta_line_length"},
                                           .desc   = "How long should each fasta line be (0: infinite)",
                                           .value  = size_t{80},
};

auto cliReadLength = clice::Argument{ .parent = &cli,
                                      .args   = {"-l", "--read_length"},
                                      .desc   = "length of the simulated reads",
                                      .value  = size_t{150},
};

auto cliNumberOfReads = clice::Argument{ .parent = &cli,
                                         .args   = {"-n", "--number_of_reads"},
                                         .desc   = "number of reads to simulate",
                                         .value  = size_t{1000},
};
auto cliErrorSubstitutions = clice::Argument{ .parent = &cli,
                                              .args   = {"--substitution_errors"},
                                              .desc   = "number of substitution errors per read",
                                              .value  = size_t{0},
};
auto cliErrorInsertions    = clice::Argument{ .parent = &cli,
                                              .args   = {"--insertion_errors"},
                                              .desc   = "number of insert errors per read",
                                              .value  = size_t{0},
};
auto cliErrorDeletions     = clice::Argument{ .parent = &cli,
                                              .args   = {"--deletion_errors"},
                                              .desc   = "number of deletion errors per read",
                                              .value  = size_t{0},
};
auto cliErrorRandom        = clice::Argument{ .parent = &cli,
                                              .args   = {"-e", "--errors"},
                                              .desc   = "number of errors (randomly chosen S, I or D)",
                                              .value  = size_t{0},
};

auto cliSeed               = clice::Argument{ .parent = &cli,
                                              .args   = {"--seed"},
                                              .desc   = "seed to initialize the random generator",
                                              .value  = (unsigned int){0},
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

auto loadFasta(std::filesystem::path input) -> std::vector<std::string> {
    auto sequences = std::vector<std::string>{};
    for (auto record : ivio::fasta::reader {{input}}) {
        auto seq = std::string{};
        seq.reserve(record.seq.size());
        for (auto c : record.seq) {
            c = ivs::dna4::normalize_char(c);
            if (!ivs::verify_char(c)) {
                c = randomPick();
            }
            seq.push_back(c);
        }
        sequences.emplace_back(std::move(seq));
    }
    return sequences;
}

std::mt19937_64 generator;

// Creates a new edit transcript of a certain length
struct Transcript {
    std::string transcript;
    size_t matches{};

    Transcript(size_t len, size_t sub=0, size_t ins=0, size_t del=0)
        : transcript(len, 'M')
        , matches{len}
    {
        addErrors(sub, ins, del);
    }

    // replaces a match with a substitution
    void addSubstitution() {
        if (matches == 0) throw std::runtime_error{"no more matches for this transcript possible"};

        auto pos = std::uniform_int_distribution<size_t>{0, transcript.size()-1}(generator);
        while (transcript[pos] != 'M') pos = std::uniform_int_distribution<size_t>{0, transcript.size()-1}(generator);
        transcript[pos] = 'S';
        matches -= 1;
    }
    // adds an insertion
    void addInsertion() {
        if (matches == 0) throw std::runtime_error{"no more matches for this transcript possible"};

        auto pos = std::uniform_int_distribution<size_t>{0, transcript.size()-1}(generator);
        while (transcript[pos] != 'M') pos = std::uniform_int_distribution<size_t>{0, transcript.size()-1}(generator);
        transcript[pos] = 'I';
        matches -= 1;
    }
    // adds an deletion
    void addDeletion() {
        auto pos = std::uniform_int_distribution<size_t>{0, transcript.size()}(generator);
        transcript.insert(transcript.begin() + pos, 'D');
    }
    void addErrors(size_t substitutions, size_t insertions, size_t deletions) {
        for (size_t i{0}; i < substitutions; ++i) addSubstitution();
        for (size_t i{0}; i < insertions; ++i)    addInsertion();
        for (size_t i{0}; i < deletions; ++i)     addDeletion();
    }
    auto lengthOfRef() const {
        size_t a = transcript.size();
        for (auto t : transcript) {
            if (t == 'I') {
                a -= 1;
            }
        }
        return a;
    }
};


struct ReadGenerator {
    std::vector<std::string> const& sequences;
    size_t readLength;
    size_t totalLength = [&]() {
        size_t l{};
        for (auto s : sequences) {
            l += s.size();
        }
        return l;
    }();
    std::mt19937_64 generator;
    std::uniform_int_distribution<size_t> uniform_pos{0, totalLength-1};

    auto generate(size_t len) -> std::tuple<size_t, size_t, std::string_view> {
        while (true) {
            // Simulating a single read
            auto pos = uniform_pos(generator);

            size_t seqId = 0;
            // pick correct sequence
            for (std::string_view seq : sequences) {
                if (pos + len > seq.size()) {
                    break;
                }
                if (pos < seq.size()) {
                    return {seqId, pos, seq.substr(pos, len)};
                }

                seqId += 1;
                pos = pos + readLength - seq.size() - 1;
            }
        }
    }

    auto applyTranscript(std::string_view v, std::span<char const> transcript) {
        std::string res;
        size_t p{0};
        auto uniform_char02 = std::uniform_int_distribution<size_t>{0, 2};
        auto uniform_char03 = std::uniform_int_distribution<size_t>{0, 3};

        auto substitute_char = [&](char c) {
            auto r = uniform_char02(generator);
            return ivs::dna4::rank_to_char((ivs::dna4::char_to_rank(c) + r + 1) % 4);
        };
        auto generate_char = [&]() {
            auto r = uniform_char03(generator);
            return ivs::dna4::rank_to_char(r);
        };

        for (auto t : transcript) {
            switch (t) {
            case 'M':
                res.push_back(v[p]);
                ++p;
                break;
            case 'S':
                res.push_back(substitute_char(v[p]));
                ++p;
                break;
            case 'I':
                res.push_back(generate_char());
                break;
            case 'D':
                ++p;
                break;
            default:
                throw error_fmt{"Invalid transcript \"{}\"", t};
            }
        }
        return res;
    }
};


void app() {
    srand(*cliSeed);
    if (cliInput) {
        auto sequences = loadFasta(*cliInput);
        fmt::print("loaded fasta file - start simulating\n");

        auto readGenerator = ReadGenerator { .sequences  = sequences,
                                             .readLength = *cliReadLength };


        auto writer = ivio::fasta::writer {{ .output = *cliOutput,
                                             .length = *cliFastaLineLength==0?std::numeric_limits<int>::max():*cliFastaLineLength,
                                           }};
        for (size_t i{0}; i < *cliNumberOfReads; ++i) {
            auto error_sub = *cliErrorSubstitutions;
            auto error_ins = *cliErrorInsertions;
            auto error_del = *cliErrorDeletions;
            for (size_t i{0}; i < *cliErrorRandom; ++i) {
                switch(rand()%3) {
                    case 0: error_sub += 1; break;
                    case 1: error_ins += 1; break;
                    case 2: error_del += 1; break;
                }
            }
            auto transcript = Transcript{*cliReadLength, error_sub, error_ins, error_del};
            auto [seqId, pos, read] = readGenerator.generate(transcript.lengthOfRef());

            auto faultyRead = readGenerator.applyTranscript(read, transcript.transcript);
            writer.write({
                .id  = fmt::format("simulated-{} (seqid:{}, pos:{}, trans:{})", i, seqId, pos, transcript.transcript),
                .seq = faultyRead,
            });
        }
    } else {
        fmt::print("no fasta file - start pure random simulating\n");
        auto writer = ivio::fasta::writer {{*cliOutput}};
        for (size_t i{0}; i < *cliNumberOfReads; ++i) {
            auto seq = std::string{};
            seq.reserve(*cliReadLength);
            for (size_t i{0}; i < *cliReadLength; ++i) {
                seq.push_back(randomPick());
            }
            writer.write({
                .id  = fmt::format("simulated-{}", i),
                .seq = seq,
            });
        }

    }
}
}
