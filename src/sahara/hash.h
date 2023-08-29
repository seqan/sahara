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
