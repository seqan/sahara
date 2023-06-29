# Utility functions
## Convert chars to ranks
1. `#!cpp void ivs::convert_char_to_rank<Alphabet>(std::span<char const> in, std::span<uint8_t> out)`
2. `#!cpp auto ivs::convert_char_to_rank<Alphabet>(std::span<char const> in) -> std::vector<uint8_t>`
3. `#!cpp auto ivs::view_char_to_rank<Alphabet> = /*unspecified*/`

Conversion from **char space** to **rank space**.
All versions are guaranteed to not throw.

### Example
```cpp
{%include-markdown "snippets/char_to_rank.cpp"%}
```
### Output
```bash
{% include-markdown "snippets/char_to_rank.cpp.out" %}
```

## Convert ranks to chars
1. `#!cpp void ivs::char_to_rank<Alphabet>(std::span<char const> in, std::span<uint8_t> out)`
2. `#!cpp auto ivs::char_to_rank<Alphabet>(std::span<char const> in) -> std::vector<uint8_t>`
3. `#!cpp auto ivs::view_char_to_rank<Alphabet> = /*unspecified*/`

## Normalize chars
1. `#!cpp void ivs::normalize_char<Alphabet>(std::span<char const> in, std::span<char> out)`
2. `#!cpp auto ivs::normalize_char<Alphabet>(std::span<char const> in) -> std::string`
3. `#!cpp auto ivs::view_normalize_char<Alphabet> = /*unspecified*/`

## Complement
1. `#!cpp void ivs::complement_rank<Alphabet>(std::span<uint8_t const> in, std::span<uint8_t> out)`
2. `#!cpp auto ivs::complement_rank<Alphabet>(std::span<uint8_t const> in) -> std::vector<uint8_t>`
3. `#!cpp auto ivs::view_complement_rank<Alphabet> = /*unspecified*/`
4. `#!cpp void ivs::complement_char<Alphabet>(std::span<char const> in, std::span<char> out)`
5. `#!cpp auto ivs::complement_char<Alphabet>(std::span<char const> in) -> std::string`
6. `#!cpp auto ivs::view_complement_char<Alphabet> = /*unspecified*/`

## Reverse complement
1. `#!cpp void ivs::reverse_complement_rank<Alphabet>(std::span<uint8_t const> in, std::span<uint8_t> out)`
2. `#!cpp auto ivs::reverse_complement_rank<Alphabet>(std::span<uint8_t const> in) -> std::vector<uint8_t>`
3. `#!cpp auto ivs::view_reverse_complement_rank<Alphabet> = /*unspecified*/`
4. `#!cpp void ivs::reverse_complement_char<Alphabet>(std::span<char const> in, std::span<char> out)`
5. `#!cpp auto ivs::reverse_complement_char<Alphabet>(std::span<char const> in) -> std::string`
6. `#!cpp auto ivs::view_reverse_complement_char<Alphabet> = /*unspecified*/`

## Verification
1. `#!cpp bool verify_char(std::span<char const> in)`
2. `#!cpp bool verify_rank(std::span<uint8_t const> in)`
