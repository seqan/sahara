// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cassert>
#include <fmindex-collection/search_scheme/all.h>
#include <set>
#include <unordered_set>
#include <fmt/ranges.h>

template <typename CB>
void allErrorConfig(fmc::search_scheme::Search s, CB cb, std::vector<int>& errorConf, size_t minError, size_t maxStep) {
    auto l = errorConf.size();
    if (l == s.pi.size()) return;

    errorConf.push_back(0);
    for (auto i{std::max(minError, s.l[l])}; i <= s.u[l]; ++i) {
        if (i-minError > maxStep) continue;
        errorConf.back() = i-minError;
        cb(errorConf);
        allErrorConfig(s, cb, errorConf, i, maxStep);
    }
    errorConf.pop_back();
}

template <typename CB>
void allErrorConfig(fmc::search_scheme::Search s, CB cb, int maxStep = std::numeric_limits<int>::max()) {
    auto errorConf = std::vector<int>{};
    allErrorConfig(s, cb, errorConf, 0, maxStep);
}

auto generateTIKZ(fmc::search_scheme::Search _s, std::vector<size_t> const& counts, bool displayAlphabet, double fontSize, bool zeroIndex) -> std::string {
    assert(isValid(_s));
    auto os = expand(_s, counts);
    assert(os);
    for (auto& v : os->pi) {
        v += 1;
    }

    auto s = limitToHamming(*os);
    for (auto& v : _s.pi) {
        v += 1;
    }


    std::string out;
    out += R"(
\begin{tikzpicture}[scale=1.]
\tikzstyle{node}=[fill=white, shape=circle, draw, minimum size=0.25cm,scale=2.]
\tikzstyle{edge}=[left,scale=1.]
\tikzstyle{medge}=[scale=1.]
\tikzstyle{redge}=[right,scale=1.]
\tikzstyle{bedge}=[below,scale=1.]
)";

    out += fmt::format("\\node[node] (n)       at (0, 0) {{}};\n");

    int leafs{0};
    auto knownErrors = std::set<std::vector<int>>{};
    int maxLevel = s.pi.size();
    allErrorConfig(s, [&](auto error) {
        knownErrors.insert(error);
        auto level = static_cast<int>(error.size());
        if (error.back() == 1) {
            auto sibError = error;
            sibError.back() = 0;
            //if (knownErrors.count(sibError) == 0) {
                leafs += 1;
            //}
        }

        auto name = fmt::format("(n{})", fmt::join(error, ""));
        auto cord = fmt::format("({:2}, {:2})", leafs, -level*2);
        out += fmt::format("\\node[node] {} at {} {{}};\n", name, cord);

/*        if (level == maxLevel) {
            leafs+= 1;
        }*/
    }, 1);

    allErrorConfig(s, [&](auto error) {
        auto level = static_cast<int>(error.size());
        auto name1 = fmt::format("(n{})", fmt::join(begin(error), --end(error), ""));
        auto name2 = fmt::format("(n{})", fmt::join(begin(error), end(error), ""));

        if (error.back() == 0) { // no error added
            char c = displayAlphabet?'M':' ';
            out += fmt::format("\\draw {} to node[edge] {{{}}} {};\n", name1, c, name2);

        } else { // error added
            char c = displayAlphabet?'S':' ';

            if (level < maxLevel) {
                out += fmt::format("\\draw[dashed] {} to node[bedge] {{{}}} {};\n", name1, c, name2);
            } else {
                out += fmt::format("\\draw[dashed] {} to node[redge] {{{}}} {};\n", name1, c, name2);
            }
        }
    }, 1);

    int  accum  = 0;
    out += fmt::format("\\node[] (sl0) at (-1, 0) {{}};\n");
    for (size_t i{1}; i < counts.size(); ++i) {
        accum += counts[_s.pi[i-1]-1];
        auto namel = fmt::format("(sl{})", i);
        auto namer = fmt::format("(sr{})", i);

        auto cordl = fmt::format("({:2}, {:2})",      -1, -accum*2);
        auto cordr = fmt::format("({:2}, {:2})",   leafs, -accum*2);
        out += fmt::format("\\node[] {} at {} {{}};\n", namel, cordl);
        out += fmt::format("\\node[] {} at {} {{}};\n", namer, cordr);
        out += fmt::format("\\draw [dashed] {} -- {};\n", namel, namer);
    }
    accum += counts.back();
    out += fmt::format("\\node[] (sl{}) at (-1, {:2}) {{}};\n", counts.size(), -accum*2);

    for (size_t i{0}; i < counts.size(); ++i) {
        out += fmt::format("\\path [] (sl{}) -- node [midway,left,scale={}] {{P{}}} (sl{});\n", i, fontSize, _s.pi[i] - ((zeroIndex)?1:0), i+1);
    }
    out += R"(
\end{tikzpicture})";


    return out;

}
