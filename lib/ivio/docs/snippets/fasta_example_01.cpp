#include <iostream>
#include <ivio/ivio.h>

int main(int argc, char** argv) {
    auto inputFile = std::filesystem::path{argv[1]};
    auto reader = ivio::fasta::reader{{.input = inputFile,
                                       .compressed = false, // false by default, if true a gzip file is expected
    }};
    for (auto record_view : reader) {
        std::cout << "id: " << record_view.id << "\n";
        std::cout << "seq: " << record_view.seq << "\n";
    }
}

