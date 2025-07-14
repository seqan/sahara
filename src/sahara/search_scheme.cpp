// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "error_fmt.h"
#include "tikz.h"

#include <clice/clice.h>
#include <fmindex-collection/search/all.h>

namespace {
void app();
auto cli = clice::Argument {
    .args   = "search_scheme",
    .desc   = "generates and info about search schemes",
    .cb     = app,
};
auto cliListGenerator = clice::Argument {
    .parent = &cli,
    .args   = {"list-generators"},
    .desc   = "show a list of generators"
};
auto cliGenerator = clice::Argument {
    .parent = &cli,
    .args   = {"-g", "--generator"},
    .desc   = "which generator to use?",
    .value  = std::string{"pigeon"},
};
auto cliQueryLength = clice::Argument {
    .parent = &cli,
    .args   = {"-l", "--length"},
    .desc   = "the assumed query length, when applying node count",
    .value  = int64_t{150},
};
auto cliReferenceLength = clice::Argument {
    .parent = &cli,
    .args   = {"--ref-length"},
    .desc   = "the assumed length of the reference text",
    .value  = int64_t{1'000'000'000},
};
auto cliMinAllowedErrors = clice::Argument {
    .parent = &cli,
    .args   = {"--min-error"},
    .desc   = "minimum errors that have to appear, such that the search scheme accepts it",
    .value  = int64_t{0},
};
auto cliMaxAllowedErrors = clice::Argument {
    .parent = &cli,
    .args   = {"-k", "--max-error"},
    .desc   = "maximum errors that can appear",
    .value  = int64_t{2},
};
auto cliAlphabetSize = clice::Argument {
    .parent = &cli,
    .args   = {"--sigma"},
    .desc   = "Size of the alphabet, e.g.: '4' for ACGT or  '5' for ACGTN",
    .value  = int64_t{4},
};
auto cliAll = clice::Argument {
    .parent = &cli,
    .args   = {"-a", "--all"},
    .desc   = "print information table about all generators",
};
auto cliYaml = clice::Argument {
    .parent = &cli,
    .args   = {"-y", "--yaml"},
    .desc   = "print in a yaml compatible format",
};
auto cliColumba = clice::Argument {
    .parent = &cli,
    .args   = {"--columba"},
    .desc   = "generates columba compatible files",
    .value  = std::filesystem::path{}
};
auto cliTikz = clice::Argument {
    .parent = &cli,
    .args   = {"--tikz"},
    .desc   = "generate a tikz diagram",
    .value  = std::string{}
};
auto cliExpansionMode = clice::Argument {
    .parent = &cli,
    .args   = {"--expansion_mode"},
    .desc   = "mode to use for generation: uniform, bottomup, topdown",
    .value  = std::string{"uniform"}
};

auto generateCounts(fmc::search_scheme::Scheme const& ss) -> std::vector<size_t> {
    if (ss.size() == 0) return {};
    if (*cliExpansionMode == "uniform") {
        return fmc::search_scheme::expandCount(ss[0].pi.size(), *cliQueryLength);
    } else if (*cliExpansionMode == "bottomup") {
        return fmc::search_scheme::optimizeByWNC(ss, *cliQueryLength, *cliAlphabetSize, *cliReferenceLength);
    } else if (*cliExpansionMode == "topdown") {
        return fmc::search_scheme::optimizeByWNC(ss, *cliQueryLength, *cliAlphabetSize, *cliReferenceLength);
    }
    throw std::runtime_error{"invalid parameter for expansion mode"};
}


void printSingleScheme() {
    // pick the correct generator
    auto iter = fmc::search_scheme::generator::all.find(*cliGenerator);
    if (iter == fmc::search_scheme::generator::all.end()) {
        throw error_fmt{"can not find generator \"{}\"", *cliGenerator};
    }
    auto const& e = iter->second;

    // generate search schemes
    auto sss = e.generator(*cliMinAllowedErrors, *cliMaxAllowedErrors, *cliAlphabetSize, *cliReferenceLength);

    // expand search schemes, so they match the length of the query
    auto ss = expand(sss, *cliQueryLength);

    // expand search schemes, but "dynamically" trying to minimize weighted node count
    auto dss = expandByWNC</*Edit=*/true>(sss, *cliQueryLength, *cliAlphabetSize, *cliReferenceLength);

    // expand search schemes, but "dynamically" trying to minimize weighted node count
    auto dss_td = expandByWNCTopDown</*Edit=*/true>(sss, *cliQueryLength, *cliAlphabetSize, *cliReferenceLength, 1);


    auto parts = sss[0].pi.size();

    fmt::print("# Search Scheme Information\n");
    fmt::print("name:                       {}\n", e.name);
    fmt::print("description:                {}\n", e.description);
    fmt::print("alphabet size:              {}\n", *cliAlphabetSize);
    fmt::print("min errors:                 {}\n", *cliMinAllowedErrors);
    fmt::print("max errors:                 {}\n", *cliMaxAllowedErrors);
    fmt::print("reference length:           {}\n", *cliReferenceLength);
    fmt::print("number of parts:            {}\n", parts);
    fmt::print("number of searches:         {}\n", ss.size());
    fmt::print("valid:                      {}\n", isValid(sss));
    fmt::print("complete:                   {}\n", isComplete(sss, *cliMinAllowedErrors, *cliMaxAllowedErrors));
    fmt::print("non-redundant:              {}\n", isNonRedundant(sss,*cliMinAllowedErrors, *cliMaxAllowedErrors));
    fmt::print("node count (ham):           {}\n", nodeCount</*Edit=*/false>(ss, *cliAlphabetSize));
    fmt::print("weighted node count (ham):  {}\n", weightedNodeCount</*Edit=*/false>(ss, *cliAlphabetSize, *cliReferenceLength));
    fmt::print("dynamic wnc (ham):          {}\n", weightedNodeCount</*Edit=*/false>(dss, *cliAlphabetSize, *cliReferenceLength));
    fmt::print("dynamic wnc td (ham):       {}\n", weightedNodeCount</*Edit=*/false>(dss_td, *cliAlphabetSize, *cliReferenceLength));
    fmt::print("node count (edit):          {}\n", nodeCount</*Edit=*/true>(ss, *cliAlphabetSize));
    fmt::print("weighted node count (edit): {}\n", weightedNodeCount</*Edit=*/true>(ss, *cliAlphabetSize, *cliReferenceLength));
    fmt::print("dynamic wnc (edit):         {}\n", weightedNodeCount</*Edit=*/true>(dss, *cliAlphabetSize, *cliReferenceLength));
    fmt::print("dynamic wnc td (edit):      {}\n", weightedNodeCount</*Edit=*/true>(dss_td, *cliAlphabetSize, *cliReferenceLength));


    fmt::print("searches:  {:^{}}  {:^{}}  {:^{}}\n", "pi", parts*3, "L", parts*3, "U", parts*3);
    for (auto const& s : sss) {
        fmt::print("           {{{}}}, {{{}}}, {{{}}}\n", fmt::join(s.pi, ", "), fmt::join(s.l, ", "), fmt::join(s.u, ", "));
    }

    fmt::print("expanded:\n");
    for (auto const& s : ss) {
        fmt::print("           {{{}}}, {{{}}}, {{{}}}\n", fmt::join(s.pi, ", "), fmt::join(s.l, ", "), fmt::join(s.u, ", "));
    }

    fmt::print("limited for hamming distance:\n");
    auto hss = limitToHamming(ss);
    for (auto const& s : hss) {
        fmt::print("           {{{}}}, {{{}}}, {{{}}}\n", fmt::join(s.pi, ", "), fmt::join(s.l, ", "), fmt::join(s.u, ", "));
    }

}

void printTikz(std::string const& path_prefix) {
    // pick the correct generator
    auto iter = fmc::search_scheme::generator::all.find(*cliGenerator);
    if (iter == fmc::search_scheme::generator::all.end()) {
        throw error_fmt{"can not find generator \"{}\"", *cliGenerator};
    }
    auto const& e = iter->second;

    // generate search schemes
    auto sss = e.generator(*cliMinAllowedErrors, *cliMaxAllowedErrors, *cliAlphabetSize, *cliReferenceLength);

    auto counts = generateCounts(sss);
    for (size_t i{0}; i < sss.size(); ++i) {
        auto filename = fmt::format("{}-{:02}.tikz", path_prefix, i);
        auto ofs = std::ofstream(filename);
        fmt::print(ofs, "{}\n", generateTIKZ(sss[i], counts, false, 4, true));
    }
}


void printTable() {
    fmt::print("# Search Scheme Information\n");
    fmt::print("alphabet size:       {}\n", *cliAlphabetSize);
    fmt::print("min errors:          {}\n", *cliMinAllowedErrors);
    fmt::print("max errors:          {}\n", *cliMaxAllowedErrors);
    fmt::print("reference length:    {}\n", *cliReferenceLength);

    fmt::print("{:^15} | {:^6} {:^8} {:^6} {:^8} {:^10} | {:^32} | {:^25} | {:^25} | {:^25}\n", "name", "parts", "searches", "valid", "complete", "non-red", "node count ham/edit", "weighted nnc ham/edit", "dyn exp bu", "dyn exp td");
    auto order = std::vector<std::string>{"backtracking", "optimum", "01*0", "01*0_opt", "pigeon", "pigeon_opt", "suffix", "h2-k1", "h2-k2", "h2-k3", "kianfar", "kucherov-k1", "kucherov-k2", "lam", "hato", "pex-td", "pex-td-l", "pex-bu", "pex-bu-l"};
    for (auto const& [key, e] : fmc::search_scheme::generator::all) {
        if (std::find(order.begin(), order.end(), key) == order.end()) {
            order.push_back(key);
            fmt::print("WARNING: missing {} in order list\n", key);
        }
    }
    for (auto const& o : order) {
        if (auto iter = fmc::search_scheme::generator::all.find(o); iter == fmc::search_scheme::generator::all.end()) {
            fmt::print("Warning: generator {} doesn't exists\n", o);
            continue;
        }

        auto sigma = *cliAlphabetSize;
        auto N     = *cliReferenceLength;

        auto const& e = fmc::search_scheme::generator::all[o];
        // generate search schemes
        auto sss = e.generator(*cliMinAllowedErrors, *cliMaxAllowedErrors, sigma, N);

        // expand search schemes, so they match the length of the query
        auto counts = generateCounts(sss);
        auto ss = expand(sss, counts);

        // expand by node count
//        auto dss_ham  = expandByNC</*Edit=*/false>(sss, *cliQueryLength, sigma);
//        auto dss_edit = expandByNC</*Edit=*/true> (sss, *cliQueryLength, sigma);

        // expand by weighted node count
        auto dess_ham  = expandByWNC</*Edit=*/false>(sss, *cliQueryLength, sigma, N);
        auto dess_edit = expandByWNC</*Edit=*/true> (sss, *cliQueryLength, sigma, N);

        // expand by weighted node count top down
        auto dess_ham_td  = expandByWNCTopDown</*Edit=*/false>(sss, *cliQueryLength, sigma, N, 1);
        auto dess_edit_td = expandByWNCTopDown</*Edit=*/true> (sss, *cliQueryLength, sigma, N, 1);

        auto parts = ss.size()>0?sss[0].pi.size():0;

        struct {
            long double countHam;
            long double countEdit;
        } stat_ss, stat_ss_w, stat_dess, stat_dess_td;

        auto valid         = isValid(sss);
        auto complete      = isComplete(sss, *cliMinAllowedErrors, *cliMaxAllowedErrors);
        auto non_redundant = isNonRedundant(sss, *cliMinAllowedErrors, *cliMaxAllowedErrors);

        stat_ss      = { .countHam  = nodeCount</*Edit=*/false>(ss, sigma),
                         .countEdit = nodeCount</*Edit=*/true>(ss, sigma)};
        stat_ss_w    = { .countHam  = weightedNodeCount</*Edit=*/false>(ss, sigma, N),
                         .countEdit = weightedNodeCount</*Edit=*/true>(ss, sigma, N)};
        stat_dess    = { .countHam  = weightedNodeCount</*Edit=*/false>(dess_ham, sigma, N),
                      .countEdit = weightedNodeCount</*Edit=*/true>(dess_edit, sigma, N)};
        stat_dess_td = { .countHam  = weightedNodeCount</*Edit=*/false>(dess_ham_td, sigma, N),
                         .countEdit = weightedNodeCount</*Edit=*/true>(dess_edit_td, sigma, N)};

        fmt::print("{:>15} | {:>6} {:>8} {:^6} {:^8} {:^10} | {:>15.0f} {:>15.0f}  | {:>12.2f} {:>12.2f} | {:>12.2f} {:>12.2f} | {:>12.2f} {:>12.2f}\n", e.name, parts, sss.size(), valid, complete, non_redundant, stat_ss.countHam, stat_ss.countEdit, stat_ss_w.countHam, stat_ss_w.countEdit, stat_dess.countHam, stat_dess.countEdit, stat_dess_td.countHam, stat_dess_td.countEdit);
    }
}

void printColumba() {
    std::filesystem::create_directories(*cliColumba);
    for (auto const& [key, e] : fmc::search_scheme::generator::all) {
        std::filesystem::create_directories(*cliColumba / key);

        // print name
        {
            auto ofs = std::ofstream{*cliColumba / key / "name.txt"};
            ofs << key;
        }
        for (auto k{*cliMinAllowedErrors}; k <= *cliMaxAllowedErrors; ++k) {
            // generate search schemes
            auto sss = e.generator(*cliMinAllowedErrors, k, *cliAlphabetSize, *cliReferenceLength);

            if (sss.empty()) continue; // no search scheme exists

            std::filesystem::create_directories(*cliColumba / key / std::to_string(k));

            auto ofs = std::ofstream{*cliColumba / key / std::to_string(k) / "searches.txt"};
            for (auto const& s : sss) {
                fmt::print(ofs, "{{{}}} {{{}}} {{{}}}\n", fmt::join(s.pi, ","), fmt::join(s.l, ","), fmt::join(s.u, ","));
            }
        }
    }
}

void printYaml() {
    fmt::print("# Search Scheme Information\n");
    fmt::print("alphabet size:       {}\n", *cliAlphabetSize);
    fmt::print("min errors:          {}\n", *cliMinAllowedErrors);
    fmt::print("max errors:          {}\n", *cliMaxAllowedErrors);
    fmt::print("reference length:    {}\n", *cliReferenceLength);
    fmt::print("---\n");
    for (auto k{*cliMinAllowedErrors}; k <= *cliMaxAllowedErrors; ++k) {
        for (auto const& [key, e] : fmc::search_scheme::generator::all) {
            // generate search schemes
            auto sss = e.generator(*cliMinAllowedErrors, k, *cliAlphabetSize, *cliReferenceLength);

            // expand search schemes, so they match the length of the query
            auto counts = generateCounts(sss);
            auto ss = expand(sss, counts);

            auto name          = e.name;
            auto parts         = ss.size()>0?sss[0].pi.size():0;
            auto searches      = ss.size();
            auto valid         = isValid(sss);
            auto complete      = isComplete(sss, *cliMinAllowedErrors, k);
            auto nodeCount     = fmc::search_scheme::nodeCount</*Edit=*/false>(ss, *cliAlphabetSize);
            auto weightedCount = weightedNodeCount</*Edit=*/false>(ss, *cliAlphabetSize, *cliReferenceLength);
            fmt::print("- name: \"{}\"\n", name);
            fmt::print("  parts: {}\n", parts);
            fmt::print("  counts: [{}]\n", fmt::join(counts, ", "));
            fmt::print("  searchCt: {}\n", searches);
            fmt::print("  valid: {}\n", valid);
            fmt::print("  complete: {}\n", complete);
            fmt::print("  nodeCount: {}\n", nodeCount);
            fmt::print("  weightedNodeCount: {:.2f}\n", weightedCount);
            fmt::print("  searches:\n");
            for (auto const& s : sss) {
                fmt::print("  - pi: [{}]\n", fmt::join(s.pi, ", "));
                fmt::print("    l: [{}]\n", fmt::join(s.l, ", "));
                fmt::print("    u: [{}]\n", fmt::join(s.u, ", "));
            }
        }
    }
}

void app() {
    if (cliListGenerator) {
        for (auto const& [key, e] : fmc::search_scheme::generator::all) {
            fmt::print("{:>15} - {}\n", e.name, e.description);
        }
        return;
    }

    if (cliAll && cliColumba) {
        printColumba();
    } else if (cliAll && cliYaml) {
        printYaml();
    } else if (cliAll) {
        printTable();
    } else if (cliTikz) {
        printTikz(*cliTikz);
    } else {
        printSingleScheme();
    }

}
}
