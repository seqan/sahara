#include <clice/clice.h>

#include <iostream>

namespace {
auto cliHelp    = clice::Argument{ .args   = "--help",
                                   .desc   = "prints the help page",
                                   .cb     = []{ std::cout << clice::generateHelp(); exit(0); },
                                 };

auto cliAdd     = clice::Argument{ .args   = "add",
                                   .symlink = true, // allows access via symlink to clice-demo-add
                                   .desc   = "adds some stuff",
                                 };
auto cliVerbose = clice::Argument{ .parent = &cliAdd,
                                   .args   = "--verbose",
                                   .desc   = "detailed description of what is happening",
                                  };
auto cliNbr     = clice::Argument{ .args   = "--nbr",
                                   .desc   = "setting some nbr",
                                   .value  = 5,
                                 };
auto cliInts    = clice::Argument{ .args   = "--ints",
                                   .desc   = "a list of numbers",
                                   .value  = std::vector<int>{},
                                 };
auto cliAuto    = clice::Argument{ .args   = "--auto",
                                   .desc   = "value depending on cliNbr +1",
                                   .value  = []() { return *cliNbr+1; },
                                   .cb     = []() {
                                        std::cout << "done\n";
                                   }
                                 };
auto cliAuto2   = clice::Argument{ .args   = "--auto2",
                                   .desc   = "same as --auto, but with lower priority",
                                   .value  = []() { return *cliNbr+1; },
                                   .cb     = []() {
                                        std::cout << "auto2 done, after auto\n";
                                   },
                                   .cb_priority = 200,
                                 };

auto cliAuto3   = clice::Argument{ .args   = "--auto3",
                                   .desc   = "same as --auto, but with higher priority",
                                   .value  = []() { return *cliNbr+1; },
                                   .cb     = []() {
                                        std::cout << "auto3 done, before auto\n";
                                   },
                                   .cb_priority = 50,
                                 };

auto cliBasicTypes = clice::Argument{ .args = "basic",
                                      .desc = "basic values like bool, ints and doubles",
                                    };
auto cliBTFlag     = clice::Argument{ .parent = &cliBasicTypes,
                                      .args   = "--flag",
                                      .desc   = "flag value",
                                    };
auto cliBTBool     = clice::Argument{ .parent = &cliBasicTypes,
                                      .args   = "--bool",
                                      .desc   = "bool value",
                                      .value  = false,
                                    };
auto cliBTUInt8   = clice::Argument{ .parent = &cliBasicTypes,
                                     .args   = "--uint8",
                                     .desc   = "uint8 value",
                                     .value  = uint8_t{},
                                   };
auto cliBTInt8    = clice::Argument{ .parent = &cliBasicTypes,
                                     .args   = "--int8",
                                     .desc   = "int8 value",
                                     .value  = int8_t{},
                                   };
auto cliBTUInt16   = clice::Argument{ .parent = &cliBasicTypes,
                                      .args   = "--uint16",
                                      .desc   = "uint16 value",
                                      .value  = uint16_t{},
                                   };
auto cliBTInt16    = clice::Argument{ .parent = &cliBasicTypes,
                                      .args   = "--int16",
                                      .desc   = "int16 value",
                                      .value  = int16_t{},
                                   };
auto cliBTUInt32   = clice::Argument{ .parent = &cliBasicTypes,
                                      .args   = "--uint32",
                                      .desc   = "uint32 value",
                                      .value  = uint32_t{},
                                   };
auto cliBTInt32    = clice::Argument{ .parent = &cliBasicTypes,
                                      .args   = "--int32",
                                      .desc   = "int32 value",
                                      .value  = int32_t{},
                                   };
auto cliBTUInt64   = clice::Argument{ .parent = &cliBasicTypes,
                                      .args   = "--uint64",
                                      .desc   = "uint64 value",
                                      .value  = uint64_t{},
                                   };
auto cliBTInt64    = clice::Argument{ .parent = &cliBasicTypes,
                                      .args   = "--int64",
                                      .desc   = "int64 value",
                                      .value  = int64_t{},
                                    };
auto cliBTChar    = clice::Argument{ .parent = &cliBasicTypes,
                                     .args   = "--char",
                                     .desc   = "char value",
                                     .value  = char{'A'},
                                    };
auto cliBTString  = clice::Argument{ .parent = &cliBasicTypes,
                                     .args   = "--string",
                                     .desc   = "string value",
                                     .value  = std::string{},
                                    };
auto cliBTMappedBool = clice::Argument{ .parent = &cliBasicTypes,
                                        .args   = "--mapped_bool",
                                        .desc   = "takes \"good\" and \"bad\" as input",
                                        .value  = bool{},
                                        .mapping = {{{"good", true}, {"bad", false}}},
                                      };

enum class MyEnumType {Foo, Bar};
auto cliBTEnum = clice::Argument{ .parent = &cliBasicTypes,
                                  .args   = "--enum",
                                  .desc   = "takes \"foo\" and \"bar\" as input",
                                  .value  = MyEnumType::Foo,
                                  .mapping = {{{"foo", MyEnumType::Foo}, {"bar", MyEnumType::Bar}}},
                                 };

auto cliBTInputPath = clice::Argument{ .parent = &cliBasicTypes,
                                       .args   = "--input",
                                       .desc   = "input path",
                                       .value  = std::filesystem::path{},
                                     };
auto cliBTOutputPath = clice::Argument{ .parent = &cliBasicTypes,
                                        .args   = "--output",
                                        .desc   = "output path",
                                        .value  = std::filesystem::path{},
                                      };
auto cliBTVectorInt = clice::Argument{ .parent = &cliBasicTypes,
                                       .args   = "--vector_int",
                                       .desc   = "vector of ints",
                                       .value  = std::vector<int>{},
                                     };

auto cliBTOption    = clice::Argument{ .parent = &cliBasicTypes,
                                       .args   = {"--option1", "--option2"},
                                       .desc   = "option can be set in two different ways",
                                       .value  = int{},
                                     };
}

auto cliSingleDash = clice::Argument{ .args = "single_dash",
                                      .desc   = "single dash grouping",
};
auto cliSD_a       = clice::Argument{ .parent = &cliSingleDash,
                                      .args   = "-a",
                                      .desc   = "some a",
};
auto cliSD_b       = clice::Argument{ .parent = &cliSingleDash,
                                      .args   = "-b",
                                      .desc   = "some b",
};


int main(int argc, char** argv) {
    try {
        if (auto failed = clice::parse(argc, argv, /*.allowDashCombi=*/true); failed) {
            std::cerr << "parsing failed: " << *failed << "\n";
            return 1;
        }

        std::cout << cliAdd << "\n";
        std::cout << "  " << cliVerbose << "\n";
        std::cout << cliHelp << "\n";
        std::cout << *cliNbr << "\n";
        std::cout << cliInts << "\n";
        for (auto i : *cliInts) {
            std::cout << " - " << i << "\n";
        }
        std::cout << "auto: " << *cliAuto << " " << *cliAuto2 << " " << *cliAuto3 << "\n";

        std::cout << "\n\nbasic " << cliBasicTypes << "\n";
        std::cout << "  --flag " << cliBTFlag << " " << "\n";
        std::cout << "  --bool " << cliBTBool << " " << *cliBTBool << "\n";
        std::cout << "  --uint8 " << cliBTUInt8 << " " << (int)*cliBTUInt8 << "\n";
        std::cout << "  --int8 " << cliBTInt8 << " " << (int)*cliBTInt8 << "\n";
        std::cout << "  --uint16 " << cliBTUInt16 << " " << *cliBTUInt16 << "\n";
        std::cout << "  --int16 " << cliBTInt16 << " " << *cliBTInt16 << "\n";
        std::cout << "  --uint32 " << cliBTUInt32 << " " << *cliBTUInt32 << "\n";
        std::cout << "  --int32 " << cliBTInt32 << " " << *cliBTInt32 << "\n";
        std::cout << "  --uint64 " << cliBTUInt64 << " " << *cliBTUInt64 << "\n";
        std::cout << "  --int64 " << cliBTInt64 << " " << *cliBTInt64 << "\n";
        std::cout << "  --char " << cliBTChar << " " << *cliBTChar << "\n";
        std::cout << "  --string " << cliBTString << " " << *cliBTString << " has len " << cliBTString->size() << "\n";
        std::cout << "  --mapped_bool " << cliBTMappedBool << " " << *cliBTMappedBool << "\n";
        std::cout << "  --enum " << cliBTEnum << " " << (*cliBTEnum==MyEnumType::Foo?"foo":"bar") << "\n";
        std::cout << "  --input " << cliBTInputPath << " " << *cliBTInputPath << "\n";
        std::cout << "  --output " << cliBTOutputPath << " " << *cliBTOutputPath << "\n";
        std::cout << "  --vector_int " << cliBTVectorInt << " " << cliBTVectorInt->size() << "\n";
        std::cout << "  --option1, --option2 " << cliBTOption << " " << *cliBTOption << "\n";
        std::cout << "\n\nsingle dash " << cliSingleDash << "\n";
        std::cout << "  -a " << cliSD_a << "\n";
        std::cout << "  -b " << cliSD_b << "\n";
        for (auto x : *cliBTVectorInt) {
            std::cout << "    " << x << "\n";
        }
    } catch (std::exception const& e) {
        std::cerr << "error: " << e.what() << "\n";
    }
}
