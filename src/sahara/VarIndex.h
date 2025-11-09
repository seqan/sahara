// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <fmindex-collection/fmindex-collection.h>

template <size_t Sigma>
struct VarIndex {
    std::string type;
    using Vs = std::variant<
        fmc::BiFMIndex<Sigma, fmc::string::InterleavedBitvector16>,
        fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_64_64k>,
        fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_512_64k>,
        typename fmc::BiFMIndex<Sigma, fmc::string::InterleavedBitvector16>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_64_64k>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_512_64k>::NoDelim
    >;
    Vs vs;

    template <typename Archive>
    void save(Archive& ar) const {
        ar(type);
        std::visit([&](auto const& v) {
            ar(v);
        }, vs);
    }
    template <typename... Args>
    void emplace(std::string _type, Args&&... args) {
        type = _type;
        if (type == "ibv16")             vs = std::variant_alternative_t<0, Vs>{std::forward<Args>(args)...};
        else if (type == "fbv64_64")     vs = std::variant_alternative_t<1, Vs>{std::forward<Args>(args)...};
        else if (type == "fbv512_64")    vs = std::variant_alternative_t<2, Vs>{std::forward<Args>(args)...};
        else if (type == "ibv16-nd")     vs = std::variant_alternative_t<3, Vs>{std::forward<Args>(args)...};
        else if (type == "fbv64_64-nd")  vs = std::variant_alternative_t<4, Vs>{std::forward<Args>(args)...};
        else if (type == "fbv512_64-nd") vs = std::variant_alternative_t<5, Vs>{std::forward<Args>(args)...};
        else throw std::runtime_error{"unknown index type"};
    }
    template <typename Archive>
    void load(Archive& ar) {
        ar(type);
        emplace(type);
        std::visit([&](auto& v) {
            ar(v);
        }, vs);
    }
};

template <>
struct VarIndex<2> {
    static constexpr size_t Sigma = 2;
    std::string type;
    using Vs = std::variant<
        typename fmc::BiFMIndex<Sigma, fmc::string::WrappedBitvectorImpl<2, fmc::bitvector::Bitvector2L<64, 65536>>::RmSigma>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::WrappedBitvectorImpl<2, fmc::bitvector::Bitvector2L<512, 65536>>::RmSigma>::NoDelim
    >;
    Vs vs;

    template <typename Archive>
    void save(Archive& ar) const {
        ar(type);
        std::visit([&](auto const& v) {
            ar(v);
        }, vs);
    }
    template <typename... Args>
    void emplace(std::string _type, Args&&... args) {
        type = _type;
        if (type == "fbv64_64-nd")  vs = std::variant_alternative_t<0, Vs>{std::forward<Args>(args)...};
        else if (type == "fbv512_64-nd") vs = std::variant_alternative_t<1, Vs>{std::forward<Args>(args)...};
        else throw std::runtime_error{"unknown index type"};
    }
    template <typename Archive>
    void load(Archive& ar) {
        ar(type);
        emplace(type);
        std::visit([&](auto& v) {
            ar(v);
        }, vs);
    }
};

