// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/fmindex/BiFMIndexKStepCursor.h>
#include <fmindex-collection/string/PairedFlattenedBitvectors2L_b.h>

template <typename Index, typename Alphabet>
auto _emplaceRev() {
    return Index{};
}

template <typename ADEntry>
//using SparseArray = fmc::suffixarray::SparseArray<ADEntry, fmc::bitvector::OptSparseRBBitvector<>>;
using SparseArray = fmc::suffixarray::SparseArray<ADEntry>;

template <typename Index, typename Alphabet>
auto _emplaceRev(fmc::Sequences auto const& _input, size_t samplingRate, size_t threadNbr) {
    bool includeReversedInput = true;
    auto [totalSize, inputText, inputSizes] = fmc::createSequences(_input, /*._addReversed=*/includeReversedInput, /*._useDelimiters=*/Index::Delim_v);
    // Adjust reversed text and compute complement
    for (size_t i{inputText.size()/2}; i < inputText.size(); ++i) {
        auto& r = inputText[i];
        r = Alphabet::complement_rank(r);
    }
    size_t refId{0};
    size_t pos{0};
    size_t seqOffset{};

    //!TODO what about empty strings?
    while (inputSizes[refId] == pos) {
        refId += 1;
        assert(refId < inputSizes.size());
    }
    auto const startRefId = refId;

    using ADEntry = std::tuple<uint32_t, uint32_t, bool>;

    auto annotatedSequence = SparseArray<ADEntry> {
        std::views::iota(size_t{0}, totalSize) | std::views::transform([&](size_t phase) -> std::optional<ADEntry> {
            if (phase == 0) {
                refId = startRefId;
                pos = 0;
            }
            if (!includeReversedInput || phase*2 < totalSize) { // going forward
                assert(refId < inputSizes.size());
                assert(pos < inputSizes[refId]);

                auto ret = std::optional<ADEntry>{std::nullopt};

                if (pos % samplingRate == 0) {
                    ret = std::make_tuple(refId+seqOffset, pos, false);
                }

                ++pos;
                if (inputSizes[refId] == pos) {
                    refId += 1;
                    pos = 0;
                }
                return ret;
            } else { // going backwards
                assert(refId < inputSizes.size());
                assert(pos < inputSizes[refId]);

                auto ret = std::optional<ADEntry>{std::nullopt};

                if (pos % samplingRate == 0) {
                    auto _refId = _input.size() + inputSizes.size() - refId-1+seqOffset;
                    size_t extra = Index::Delim_v?1:0;
                    auto _pos   = (inputSizes[refId] - pos + inputSizes[refId] - 1 - extra) % inputSizes[refId];
                    ret = std::make_tuple(_refId, _pos, true);
                }

                ++pos;
                if (inputSizes[refId] == pos) {
                    refId += 1;
                    pos = 0;
                }
                return ret;
            }
        })
    };
    return Index{inputText, annotatedSequence, threadNbr, /*includeReversedInput=*/false};
}

template <typename Alphabet, size_t Sigma=Alphabet::size()>
struct VarIndex {
    size_t sigma{Sigma};
    size_t samplingRate;
    std::string type;
    using Vs = std::variant<
        fmc::BiFMIndex<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t>>>,
        fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>,
        fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>,
        fmc::BiFMIndex<Sigma, fmc::string::PairedFlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>,
        fmc::BiFMIndex<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>,
        typename fmc::BiFMIndex<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::PairedFlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim,
        fmc::BiFMIndex<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t, bool>>>,
        fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t, bool>>>,
        fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t, bool>>>,
        typename fmc::BiFMIndex<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t, bool>>>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t, bool>>>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t, bool>>>::NoDelim,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<2>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<3>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::InterleavedBitvector16, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::FlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_64_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<4>,
        typename fmc::BiFMIndexKStep<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim::template SetKStep<4>
    >;
    Vs vs;

    template <typename Archive>
    void save(Archive& ar) const {
        ar(sigma);
        ar(samplingRate);
        ar(type);
        std::visit([&](auto const& v) {
            ar(v);
        }, vs);
    }

    template <size_t I, typename ...Args>
    void _emplace(Args&&... args) {
        using Index = std::variant_alternative_t<I, Vs>;
        vs = Index{std::forward<Args>(args)...};
    }


    template <size_t I, typename ...Args>
    void _emplaceRev(Args&&... args) {
        using Index = std::variant_alternative_t<I, Vs>;
        vs = ::_emplaceRev<Index, Alphabet>(std::forward<Args>(args)...);
    }

    template <typename... Args>
    void emplace(std::string _type, Args&&... args) {
        type = _type;
        if (type == "ibv16"             || type=="ibv16_1step")           _emplace< 0>(std::forward<Args>(args)...);
        else if (type == "fbv64_64"     || type == "fbv64_64_1step")      _emplace< 1>(std::forward<Args>(args)...);
        else if (type == "fbv512_64"    || type == "fbv512_64_1step")     _emplace< 2>(std::forward<Args>(args)...);
        else if (type == "pfbv64_64"    || type == "pfbv64_64_1step")     _emplace< 3>(std::forward<Args>(args)...);
        else if (type == "pfbv512_64"   || type == "pfbv512_64_1step")    _emplace< 4>(std::forward<Args>(args)...);
        else if (type == "ibv16-nd"     || type == "ibv16_1step-nd"     ) _emplace< 5>(std::forward<Args>(args)...);
        else if (type == "fbv64_64-nd"  || type == "fbv64_64_1step-nd"  ) _emplace< 6>(std::forward<Args>(args)...);
        else if (type == "fbv512_64-nd" || type == "fbv512_64_1step-nd" ) _emplace< 7>(std::forward<Args>(args)...);
        else if (type == "pfbv64_64-nd" || type == "pfbv64_64_1step-nd" ) _emplace< 8>(std::forward<Args>(args)...);
        else if (type == "pfbv512_64-nd"|| type == "pfbv512_64_1step-nd") _emplace< 9>(std::forward<Args>(args)...);
        else if (type == "ibv16-rev")           _emplaceRev<10>(std::forward<Args>(args)...);
        else if (type == "fbv64_64-rev")        _emplaceRev<11>(std::forward<Args>(args)...);
        else if (type == "fbv512_64-rev")       _emplaceRev<12>(std::forward<Args>(args)...);
        else if (type == "ibv16-nd-rev")        _emplaceRev<13>(std::forward<Args>(args)...);
        else if (type == "fbv64_64-nd-rev")     _emplaceRev<14>(std::forward<Args>(args)...);
        else if (type == "fbv512_64-nd-rev")    _emplaceRev<15>(std::forward<Args>(args)...);
        else if (type == "ibv16_2step")         _emplace<16>(std::forward<Args>(args)...);
        else if (type == "fbv64_64_2step")      _emplace<17>(std::forward<Args>(args)...);
        else if (type == "fbv512_64_2step")     _emplace<18>(std::forward<Args>(args)...);
        else if (type == "pfbv64_64_2step")     _emplace<19>(std::forward<Args>(args)...);
        else if (type == "pfbv512_64_2step")    _emplace<20>(std::forward<Args>(args)...);
        else if (type == "ibv16_3step")         _emplace<21>(std::forward<Args>(args)...);
        else if (type == "fbv64_64_3step")      _emplace<22>(std::forward<Args>(args)...);
        else if (type == "fbv512_64_3step")     _emplace<23>(std::forward<Args>(args)...);
        else if (type == "fbv64_64_3step")      _emplace<24>(std::forward<Args>(args)...);
        else if (type == "fbv512_64_3step")     _emplace<25>(std::forward<Args>(args)...);
        else if (type == "ibv16_4step")         _emplace<26>(std::forward<Args>(args)...);
        else if (type == "fbv64_64_4step")      _emplace<27>(std::forward<Args>(args)...);
        else if (type == "pfbv512_64_4step")    _emplace<28>(std::forward<Args>(args)...);
        else if (type == "pfbv64_64_4step")     _emplace<29>(std::forward<Args>(args)...);
        else if (type == "fbv512_64_4step")     _emplace<30>(std::forward<Args>(args)...);
        else if (type == "ibv16_2step-nd")      _emplace<31>(std::forward<Args>(args)...);
        else if (type == "fbv64_64_2step-nd")   _emplace<32>(std::forward<Args>(args)...);
        else if (type == "fbv512_64_2step-nd")  _emplace<33>(std::forward<Args>(args)...);
        else if (type == "pfbv64_64_2step-nd")  _emplace<34>(std::forward<Args>(args)...);
        else if (type == "pfbv512_64_2step-nd") _emplace<35>(std::forward<Args>(args)...);
        else if (type == "ibv16_3step-nd")      _emplace<36>(std::forward<Args>(args)...);
        else if (type == "fbv64_64_3step-nd")   _emplace<37>(std::forward<Args>(args)...);
        else if (type == "fbv512_64_3step-nd")  _emplace<38>(std::forward<Args>(args)...);
        else if (type == "pfbv64_64_3step-nd")  _emplace<39>(std::forward<Args>(args)...);
        else if (type == "pfbv512_64_3step-nd") _emplace<40>(std::forward<Args>(args)...);
        else if (type == "ibv16_4step-nd")      _emplace<41>(std::forward<Args>(args)...);
        else if (type == "fbv64_64_4step-nd")   _emplace<42>(std::forward<Args>(args)...);
        else if (type == "fbv512_64_4step-nd")  _emplace<43>(std::forward<Args>(args)...);
        else if (type == "pfbv64_64_4step-nd")  _emplace<44>(std::forward<Args>(args)...);
        else if (type == "pfbv512_64_4step-nd") _emplace<45>(std::forward<Args>(args)...);
        else throw std::runtime_error{"unknown index type: " + type};
    }
    template <typename Archive>
    void load(Archive& ar) {
        ar(sigma);
        ar(samplingRate);
        ar(type);
        emplace(type);
        std::visit([&](auto& v) {
            ar(v);
        }, vs);
    }
    template <typename Archive>
    void saveSize(Archive& ar) const {
        ar(sigma);
        ar(samplingRate);
        ar(type);
        std::visit([&](auto& v) {
            ar(v);
        }, vs);
    }
};

template <typename Alphabet>
struct VarIndex<Alphabet, 2> {
    static constexpr size_t Sigma = 2;
    size_t sigma{Sigma};
    size_t samplingRate;
    std::string type;
    using Vs = std::variant<
        typename fmc::BiFMIndex<Sigma, fmc::string::WrappedBitvectorImpl<2, fmc::bitvector::Bitvector2L<64, 65536>>::RmSigma, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::WrappedBitvectorImpl<2, fmc::bitvector::Bitvector2L<512, 65536>>::RmSigma, SparseArray<std::tuple<uint32_t, uint32_t>>>::NoDelim,
        typename fmc::BiFMIndex<Sigma, fmc::string::WrappedBitvectorImpl<2, fmc::bitvector::Bitvector2L<64, 65536>>::RmSigma, SparseArray<std::tuple<uint32_t, uint32_t, bool>>>::NoDelim::ReuseRev,
        typename fmc::BiFMIndex<Sigma, fmc::string::WrappedBitvectorImpl<2, fmc::bitvector::Bitvector2L<512, 65536>>::RmSigma, SparseArray<std::tuple<uint32_t, uint32_t, bool>>>::NoDelim::ReuseRev
    >;
    Vs vs;

    template <typename Archive>
    void save(Archive& ar) const {
        ar(sigma);
        ar(samplingRate);
        ar(type);
        std::visit([&](auto const& v) {
            ar(v);
        }, vs);
    }

    template <size_t I, typename ...Args>
    void _emplace(Args&&... args) {
        using Index = std::variant_alternative_t<I, Vs>;
        vs = Index{std::forward<Args>(args)...};
    }

    template <size_t I, typename ...Args>
    void _emplaceRev(Args&&... args) {
        using Index = std::variant_alternative_t<I, Vs>;
        vs = ::_emplaceRev<Index, Alphabet>(std::forward<Args>(args)...);
    }

    template <typename... Args>
    void emplace(std::string _type, Args&&... args) {
        type = _type;
        if (type == "fbv64_64-nd")           _emplace<0>(std::forward<Args>(args)...);
        else if (type == "fbv512_64-nd")     _emplace<1>(std::forward<Args>(args)...);
        else if (type == "fbv64_64-nd-rev")  _emplaceRev<2>(std::forward<Args>(args)...);
        else if (type == "fbv512_64-nd-rev") _emplaceRev<3>(std::forward<Args>(args)...);
        else throw std::runtime_error{"unknown index type: " + type};
    }
    template <typename Archive>
    void load(Archive& ar) {
        ar(sigma);
        ar(samplingRate);
        ar(type);
        emplace(type);
        std::visit([&](auto& v) {
            ar(v);
        }, vs);
    }
    template <typename Archive>
    void saveSize(Archive& ar) const {
        ar(sigma);
        ar(samplingRate);
        ar(type);
        std::visit([&](auto& v) {
            ar(v);
        }, vs);
    }
};

