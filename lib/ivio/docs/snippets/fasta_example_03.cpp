#include <iostream>
#include <ivio/ivio.h>

int main(int, char**) {
    auto reader = ivio::fasta::reader{{.input = std::cin,
                                       .compressed = false, // false by default, if true a gzip file is expected
    }};

    // creates a vector of type `std::vector<ivio::fasta::record>` even though reader returns `record_view`
    auto vec = std::vector{begin(reader), end(reader)};
}

