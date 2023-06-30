#include <ivsigma/ivsigma.h>
#include <iostream>

int main()
{
    auto checkString = [](std::string_view s) {
        auto r = ivs::verify_char(s);
        if (r) {
            std::cout << "error, unexpected char " << (int)s[*r] << " at position " << *r << '\n';
        } else {
            std::cout << "no error in '" << s << "'\n";
        }
    };

    using namespace std::literals;
    checkString("ACGT"s);   // doesn't fail
    checkString("AC\0GT"s); // will fail

    auto checkRank = [](std::span<uint8_t const> v) {
        auto r = ivs::verify_rank(v);
        if (r) {
            std::cout << "error, unexpected rank " << (int)v[*r] << " at position " << *r << '\n';
        } else {
            std::cout << "no error\n";
        }
    };

    checkRank(std::vector<uint8_t>{0, 1, 2, 3});      // doesn't fail
    checkRank(std::vector<uint8_t>{0, 1, 2, 3, 255}); // will fail

}
