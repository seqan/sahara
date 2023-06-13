#include <iostream>
#include <ivio/ivio.h>

int main(int, char**) {
    auto reader = ivio::fasta::reader{{.input = std::cin}};
    auto writer = ivio::fasta::writer{{.input = std::cout}};

    for (auto record_view : reader) {
        writer.write(record_view);
    }
}

