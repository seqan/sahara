# IVio
**IVio**, pronounced **for you** aka **io4**.
A c++20 library for **bioinformatics**. Providing functions for reading and writing different file formats.

It supports many file like **fasta**, **fastq**, **faidx**, **vcf**, **bcf**, **sam**, **bam**.

See [documentation](https://iv-project.github.io/IVio/) for detailed information.

## Usage CMake
To use this package with cmake you must add this repository via `find_package` and add a `target_link_libraries`.
e.g:
```cmake
find_package(ivio PATHS lib/IVio)
add_executable(your_exec file1.cpp)
target_link_libraries(your_exec PUBLIC ivio::ivio)
```

## Benchmarks

For detailed benchmarks see [IVbenchmark](https://github.com/iv-project/IVbenchmark#results).

## Sample snippet
This demonstrates how reading a file and printing it to the command line looks like:
```c++
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
```
