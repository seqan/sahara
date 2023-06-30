# Utility functions

---
## Convert chars to ranks
1. `#!cpp void ivs::convert_char_to_rank<Alphabet>(std::span<char const> in, std::span<uint8_t> out)`
2. `#!cpp auto ivs::convert_char_to_rank<Alphabet>(std::span<char const> in) -> std::vector<uint8_t>`
3. `#!cpp auto ivs::view_char_to_rank<Alphabet> = /*unspecified*/`

Conversion from **char space** to **rank space**. Version 1 and 3 will never throw. Version 2 might throw inside of `std::vector`.
Characters invalid to `Alphabet` will be converted to `255` and can be checked via `verify_rank`.

### Example
```cpp
{% include-markdown "snippets/char_to_rank.cpp" %}
```
**Output:**
```bash
{% include-markdown "snippets/char_to_rank.cpp.out" %}
```

---
## Convert ranks to chars
1. `#!cpp void ivs::convert_rank_to_char<Alphabet>(std::span<char const> in, std::span<uint8_t> out)`
2. `#!cpp auto ivs::convert_rank_to_char<Alphabet>(std::span<char const> in) -> std::string`
3. `#!cpp auto ivs::view_rank_to_charchar<Alphabet> = /*unspecified*/`

Conversion from **rank space** to **char space**.
Version 1 and 3 will never throw. Version 2 might throw inside of `std::string`.
Ranks invalid to `Alphabet` will be converted to `\0` and can be checked via `verify_char`.

### Example
```cpp
{% include-markdown "snippets/rank_to_char.cpp" %}
```
**Output:**
```bash
{% include-markdown "snippets/rank_to_char.cpp.out" %}
```

---
## Normalize chars
1. `#!cpp void ivs::normalize_char<Alphabet>(std::span<char const> in, std::span<char> out)`
2. `#!cpp auto ivs::normalize_char<Alphabet>(std::span<char const> in) -> std::string`
3. `#!cpp auto ivs::view_normalize_char<Alphabet> = /*unspecified*/`

Normalizes string according to `Alphabet`. For example using `dna4` and providing the string `AaCcGgTt` will convert it
to `AACCGGTT`. Additionally some alphabets have some special rules. For `dna4` the letter `U` will also be normalized to `T`.
Version 1 and 3 will never throw. Version 2 might throw.
Characters invalid to `Alphabet` will be converted to `255` and can be checked via `verify_rank`.

### Example
```cpp
{% include-markdown "snippets/normalize_char.cpp" %}
```
**Output:**
```bash
{% include-markdown "snippets/normalize_char.cpp.out" %}
```

---
## Complement
1. `#!cpp void ivs::complement_rank<Alphabet>(std::span<uint8_t const> in, std::span<uint8_t> out)`
2. `#!cpp auto ivs::complement_rank<Alphabet>(std::span<uint8_t const> in) -> std::vector<uint8_t>`
3. `#!cpp auto ivs::view_complement_rank<Alphabet> = /*unspecified*/`
4. `#!cpp void ivs::complement_char<Alphabet>(std::span<char const> in, std::span<char> out)`
5. `#!cpp auto ivs::complement_char<Alphabet>(std::span<char const> in) -> std::string`
6. `#!cpp auto ivs::view_complement_char<Alphabet> = /*unspecified*/`

Computes the complement according to `Alphabet`. It is required that `Alphabet` has the concept `alphabet_with_complement_c`.
Version 1, 3, 4 and 6 never throw. Version 2 and 5 might throw inside of `std::vector` or `std::string`.
Invalid ranks in `*_rank` functions will be converted to `0`.
in the `*_char` version, invalid letters will be converted to `\0`.

### Example
```cpp
{% include-markdown "snippets/complement.cpp" %}
```
**Output:**
```bash
{% include-markdown "snippets/complement.cpp.out" %}
```

---
## Reverse complement
1. `#!cpp void ivs::reverse_complement_rank<Alphabet>(std::span<uint8_t const> in, std::span<uint8_t> out)`
2. `#!cpp auto ivs::reverse_complement_rank<Alphabet>(std::span<uint8_t const> in) -> std::vector<uint8_t>`
3. `#!cpp auto ivs::view_reverse_complement_rank<Alphabet> = /*unspecified*/`
4. `#!cpp void ivs::reverse_complement_char<Alphabet>(std::span<char const> in, std::span<char> out)`
5. `#!cpp auto ivs::reverse_complement_char<Alphabet>(std::span<char const> in) -> std::string`
6. `#!cpp auto ivs::view_reverse_complement_char<Alphabet> = /*unspecified*/`

Computes the reverse complement according to `Alphabet`. It is required that `Alphabet` has the concept `alphabet_with_complement_c`.
Version 1, 3, 4 and 6 never throw. Version 2 and 5 might throw inside of `std::vector` or `std::string`.
Invalid ranks in `*_rank` functions will be converted to `0`.
In the `*_char` version, invalid letters will be converted to `\0`.

### Example
```cpp
{% include-markdown "snippets/reverse_complement.cpp" %}
```
**Output:**
```bash
{% include-markdown "snippets/reverse_complement.cpp.out" %}
```

---
## Verification
1. `#!cpp std::optional<size_t> verify_char(std::span<char const> in)`
2. `#!cpp std::optional<size_t> verify_rank(std::span<uint8_t const> in)`

Verifies if characters and ranks are valid. Invalid ranks are `255` and invalid characters `\0`. If an invalid
char/rank is found, this function reports its position. These function never throw.

### Example
```cpp
{% include-markdown "snippets/reverse_complement.cpp" %}
```
**Output:**
```bash
{% include-markdown "snippets/reverse_complement.cpp.out" %}
```

---
