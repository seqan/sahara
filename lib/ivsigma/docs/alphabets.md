# Alphabets
The following alphabets are available

| Nucleotide based | Amino acid based | Quality based         |
|------------------|------------------|-----------------------|
| `ivs::dna4`      | `ivs::aa27`      | `ivs::pthred42`       |
| `ivs::dna5`      | `ivs::aa20`      | `ivs::pthred63`       |
| `ivs::rna4`      | `ivs::aa10li`    | `ivs::pthred68solexa` |
| `ivs::rna5`      | `ivs::aa10murphy`| `ivs::pthred94`       |
| `ivs::iupac`     |                  |                       |
| `ivs::dna3bs`    |                  |                       |

## Nucleotide based
- `#!cpp uint8_t char_to_rank(char)`
- `#!cpp char rank_to_char(uint8_t)`
- `#!cpp char normalize_char(char)`
- `#!cpp char complement_char(char)`
- `#!cpp uint8_t complement_rank(uint8_t)`

All functions are guarenteed to not throw.

## Amino acid based
- `#!cpp uint8_t char_to_rank(char)`
- `#!cpp char rank_to_char(uint8_t)`
- `#!cpp char normalize_char(char)`

All functions are guarenteed to not throw.

## Quality based
- `#!cpp uint8_t char_to_rank(char)`
- `#!cpp char rank_to_char(uint8_t)`
- `#!cpp char normalize_char(char)`

All functions are guarenteed to not throw.
