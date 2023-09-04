#pragma once
#include <ivsigma/ivsigma.h>

struct dr_dna4 : ivs::simple_alphabet<    ivs::rank_char_mapping<   0,    '$'>,
                                          ivs::rank_char_mapping<   1,    'A',   'a', 'T', 't', 'U', 'u'>,
                                          ivs::rank_char_mapping<   2,    'C',   'c', 'G', 'g'>>
{};
static_assert(ivs::alphabet_c<dr_dna4>, "Unit test: is supposed to model an alphabet");

struct dr_dna5 : ivs::simple_alphabet<    ivs::rank_char_mapping<   0,    '$'>,
                                          ivs::rank_char_mapping<   1,    'A',   'a', 'T', 't', 'U', 'u'>,
                                          ivs::rank_char_mapping<   2,    'C',   'c', 'G', 'g'>,
                                          ivs::rank_char_mapping<   3,    'N',   'n'>>
{};
static_assert(ivs::alphabet_c<dr_dna5>, "Unit test: is supposed to model an alphabet");
