# FASTA

Support reading and writing [fasta files](https://blast.ncbi.nlm.nih.gov/doc/blast-topics/).
All classes/functions are available in the `ivio::fasta` namespace.

[TOC]

## Record and RecordView (`record`, `record_view`) { #markdown data-toc-label="Record and RecordView" }
IVio provides the datastructures `ivio::fasta::record` and `ivio::fasta::record_view`. These are the typical input/output arguments of the readers and writers.
These datastructures convertible to each other.

The `record_view` layout
``` c++
struct record_view {
    std::string_view id;
    std::string_view seq;
};
```
The `record` layout
```c++
struct record {
    std::string id;
    std::string seq;
};
```

## Reader (`reader`, `reader::config`) { #markdown data-toc-label="Reader" }
The `ivio::fasta::reader` accepts a `ivio::fasta::reader::config` object for initialization. It initializes an object
that enables you to iterate over the records of a file. It is is for-range compatible, meaning it fulfils c++ concepts of  [range](https://en.cppreference.com/w/cpp/ranges/range) and [LegacyInputIterator](https://en.cppreference.com/w/cpp/named_req/InputIterator).
When looping over a reader it returns `record_view` objects that are only valid until next record is being requested from the reader.
To get ownership of the data, it is required to create an object of type `record`.

To configure the `reader` one must fill the `ivio::fasta::reader::config` struct which declares the input file and if gzip compression is expected.
```c++
struct ivio::fasta::reader::config {
    // Source: file or stream
    std::variant<std::filesystem::path, std::reference_wrapper<std::istream>> input;

    // This is only relevant if a stream is being used
    bool compressed{};
};
```

Overview of the memberfunctions of `ivio::fasta::reader`
```c++
struct ivio::fasta::reader {
    /*...*/
    auto next() -> std::optional<record_view>; // reads the next record
    void close() const; // closes the read file
};
```

## Writer (`writer`, `writer::config`) { #markdown data-toc-label="Writer" }
The `ivio::fasta::writer` provides a single function `write` which takes a `ivio::fasta::record_view` as input.
The class is initialized with a `ivio::fasta::writer::config` object which has the options:
```
struct ivio::fasta::writer::config {
    // Source: file or stream
    std::variant<std::filesystem::path, std::reference_wrapper<std::ostream>> output;

    // This is only relevant if a stream is being used
    bool compressed{};

    size_t length{80}; // Break after 80 characters
};
```
## Examples
### Example - Reading record by record
In this example a file is being opened and print to command line
```c++
{% include-markdown "snippets/fasta_example_01.cpp" %}
```
**Output**
```sh
{% include-markdown "snippets/fasta_example_01.cpp.out" %}
```
### Example - Copying a file
In this example the data is read from the standard input and writen to standard output
```c++
{% include-markdown "snippets/fasta_example_02.cpp" %}
```

### Example - Reading complete file
Load complete fasta file into memory:
```c++
{% include-markdown "snippets/fasta_example_03.cpp" %}
```
