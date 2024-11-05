// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <clice/clice.h>
#include <fmt/format.h>

namespace {
auto cliHelp = clice::Argument { .args     = {"-h", "--help"},
                                 .desc     = "prints the help page",
                                 .cb       = []{ fmt::print("{}", clice::generateHelp()); exit(0); },
};
}

void slix_script_main(std::string script);

int main(int argc, char** argv) {
    try {
        if (auto failed = clice::parse(argc, argv); failed) {
            fmt::print(stderr, "parsing failed {}\n", *failed);
            return 1;
        }
        if (auto ptr = std::getenv("CLICE_COMPLETION"); ptr) {
            return 0;
        }
    } catch (std::exception const& e) {
        fmt::print(stderr, "error {}\n", e.what());
    }
}
