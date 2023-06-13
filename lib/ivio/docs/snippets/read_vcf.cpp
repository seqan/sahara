#include <ivio/ivio.h>
#include <iostream>

int main(int argc, char** argv) {
    auto file = std::filesystem::path{argv[1]};
    for (auto view : ivio::vcf::reader{{file}}) {
        std::cout << view.ref << "\n";
    }
}
