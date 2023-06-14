# clice
A CLI argument parser for c++20. It is the spiritual successor of sargparse. It shares many of the goals
and design principles.


```c++
#include <clice/clice.h>
#include <iostream>

auto cliAdd     = clice::Argument{ .arg    = "add",
                                   .desc   = "adds some stuff",
                                 };
auto cliVerbose = clice::Argument{ .parent = &cliAdd,
                                   .arg    = "--verbose",
                                   .desc   = "detailed description of what is happening",
                                  };
int main(int argc, char** argv) {
    if (auto failed = clice::parse(argc, argv); failed) {
        std::cerr << "parsing failed: " << *failed << "\n";
        return 1;
    }
    if (auto ptr = std::getenv("CLICE_COMPLETION"); ptr) {
        return 0;
    }


    std::cout << cliAdd << "\n";
    std::cout << "  " << cliVerbose << "\n";
    return 0;
}
```

# Bash/Zsh completion
Just run `eval "$(CLICE_GENERATE_COMPLETION=$$ ./clice-demo)"` and enjoy
tab-completion when running `./clice-demo` programs.
