# K-mers/Minimizers
K-mers (also known as q-grams) are substrings of length `k`.
These allow for a convenient handling of long string, by splitting it into overlapping
substrings and perform analysis on these substrings.


Shorter k-mers can even be placed into a 64bit integer, allowing for quick comparision with other k-mers.

## K-mers
```
    template <alphabet_c Alphabet>
    struct compact_encoding;
```

This structure creates a view over a rank based values and encodes k-mers
according to some `Alphabet`. If `Alphabet` allow for `complement` this encoding will take the canonical kmer.
It allows for typical for-range loop syntax.

### Example
```cpp
{% include-markdown "snippets/compact_encoding.cpp" %}
```
**Output:**
```bash
{% include-markdown "snippets/compact_encoding.cpp.out" %}
```

---
## Winnowing minimizers
Winnowing minimizers select a smaller subset from all k-mers. These are representatives of all the k-mers of a sequence.
Dealing with a smaller set of k-mers allows faster analysis of sequences. These work similar to the `compact_encoding` view.
Additionally it takes a window, from which it picks the smallest value.

### Example
```cpp
{% include-markdown "snippets/winnowing_minimizers.cpp" %}
```
**Output:**
```bash
{% include-markdown "snippets/winnowing_minimizers.cpp.out" %}
```


