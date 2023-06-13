// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "simple_alphabet.h"

namespace ivs {

//                                                             rank, symbol, compl, alts...
struct dna4 : alphabet_with_compl<rank_char_mapping_with_compl<   0,    'A',   'T', 'a'>,
                                  rank_char_mapping_with_compl<   1,    'C',   'G', 'c'>,
                                  rank_char_mapping_with_compl<   2,    'G',   'C', 'g'>,
                                  rank_char_mapping_with_compl<   3,    'T',   'A', 't'>>
{};
static_assert(alphabet_c<dna4>, "Unit test: is supposed to model an alphabet");

struct dna5 : alphabet_with_compl<rank_char_mapping_with_compl<   0,    'A',   'T', 'a'>,
                                  rank_char_mapping_with_compl<   1,    'C',   'G', 'c'>,
                                  rank_char_mapping_with_compl<   2,    'G',   'C', 'g'>,
                                  rank_char_mapping_with_compl<   3,    'T',   'A', 't'>,
                                  rank_char_mapping_with_compl<   4,    'N',   'N', 'n'>>
{};
static_assert(alphabet_c<dna5>, "Unit test: is supposed to model an alphabet");

struct iupac : alphabet_with_compl<rank_char_mapping_with_compl<   0,   'A',   'T', 'a'>,
                                   rank_char_mapping_with_compl<   1,   'C',   'G', 'c'>,
                                   rank_char_mapping_with_compl<   2,   'G',   'C', 'g'>,
                                   rank_char_mapping_with_compl<   3,   'T',   'A', 't', 'U', 'u'>,
                                   rank_char_mapping_with_compl<   4,   'N',   'N', 'n'>,
                                   rank_char_mapping_with_compl<   5,   'R',   'Y', 'r'>,
                                   rank_char_mapping_with_compl<   6,   'Y',   'R', 'y'>,
                                   rank_char_mapping_with_compl<   7,   'S',   'S', 's'>,
                                   rank_char_mapping_with_compl<   8,   'W',   'W', 'w'>,
                                   rank_char_mapping_with_compl<   9,   'K',   'M', 'k'>,
                                   rank_char_mapping_with_compl<  10,   'M',   'K', 'm'>,
                                   rank_char_mapping_with_compl<  11,   'B',   'V', 'b'>,
                                   rank_char_mapping_with_compl<  12,   'D',   'H', 'd'>,
                                   rank_char_mapping_with_compl<  13,   'H',   'D', 'h'>,
                                   rank_char_mapping_with_compl<  14,   'V',   'B', 'v'>,
                                   rank_char_mapping_with_compl<  15,   '-',   '-', '.'>>
{};
static_assert(alphabet_c<iupac>, "Unit test: is supposed to model an alphabet");

struct rna4 : alphabet_with_compl<rank_char_mapping_with_compl<   0,    'A',   'U', 'a'>,
                                  rank_char_mapping_with_compl<   1,    'C',   'G', 'c'>,
                                  rank_char_mapping_with_compl<   2,    'G',   'C', 'g'>,
                                  rank_char_mapping_with_compl<   3,    'U',   'A', 'u'>>
{};
static_assert(alphabet_c<rna4>, "Unit test: is supposed to model an alphabet");

struct rna5 : alphabet_with_compl<rank_char_mapping_with_compl<   0,    'A',   'U', 'a'>,
                                  rank_char_mapping_with_compl<   1,    'C',   'G', 'c'>,
                                  rank_char_mapping_with_compl<   2,    'G',   'C', 'g'>,
                                  rank_char_mapping_with_compl<   3,    'U',   'A', 'u'>,
                                  rank_char_mapping_with_compl<   4,    'N',   'N', 'n'>>
{};
static_assert(alphabet_c<rna5>, "Unit test: is supposed to model an alphabet");

//                                                rank, symbol,   alts...
struct dna3bs : simple_alphabet<rank_char_mapping<   0,    'A',   'a'>,
                                rank_char_mapping<   1,    'G',   'g'>,
                                rank_char_mapping<   2,    'T',   't', 'C', 'c'>>
{};
static_assert(alphabet_c<dna3bs>, "Unit test: is supposed to model an alphabet");

}
