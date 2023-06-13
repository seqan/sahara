# VCF

Support reading and writing [vcf files](https://samtools.github.io/hts-specs/VCFv4.3.pdf)

## Record and Record-View
IVio provides the data structs `ivio::vcf::record` and `ivio::vcf::record_view`.
A `record_view` has following layout:
``` c++
struct record_view {
    std::string_view            chrom;
    int32_t                     pos;
    std::string_view            id;
    std::string_view            ref;
    std::string_view            alts;
    std::optional<float>        qual;
    std::string_view            filters;
    std::string_view            infos;
    std::string_view            formats;
    std::string_view            samples;
};
```
The type `record` looks like this:
```c++
struct record {
    std::string                 chrom;
    int32_t                     pos;
    std::string                 id;
    std::string                 ref;
    std::string                 alts;
    std::optional<float>        qual;
    std::string                 filters;
    std::string                 infos;
    std::string                 formats;
    std::string                 samples;
};
```

## Reading
**Example**
```cpp
{% include-markdown "snippets/read_vcf.cpp" %}
```
**Output**
```sh
{% include-markdown "snippets/read_vcf.cpp.out" %}
```

## Writing
**Example**
```cpp
{% include-markdown "snippets/write_vcf.cpp" %}
```
**Output**
```cpp
{% include-markdown "snippets/write_vcf.cpp.out" %}
```
