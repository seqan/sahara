// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <xxhash.h>

/*struct xxhash {
    XXH64_state_t* state;

    xxhash() {
        state = XXH64_createState();
    }
    ~xxhash() {
    }
    auto hash(uint64_t v) const -> uint64_t {
        XXH64_hash const seed = 0;
        XXH64_reset(state, seed);
        XXH64_update(state, &v, sizeof(v));
        XXH64_hash_t const hash = XXH64_digest(state);
        return hash;
    }
};*/
inline auto hash(uint64_t v) -> uint64_t {
    return XXH64(&v, sizeof(v), uint64_t{0});
}
