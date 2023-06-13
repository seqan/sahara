#pragma once

#include <filesystem>
#include <ostream>

namespace io2 {

/** Wrapper to allow path and stream inputs
 */
template <typename FileOut, auto format = nullptr>
struct Output {
    std::optional<FileOut> fileOut;

    Output() {}

    Output(char const* _path) {
        fileOut.emplace(_path);
    }

    Output(std::string const& _path)
        : Output(_path.c_str())
    {}

    Output(std::filesystem::path const& _path)
        : Output(_path.c_str())
    {}

    template <typename format_t>
    Output(std::ostream& ostr, format_t _format) {
        convert_format(_format, [&](auto _format) {
            if (!open(*fileOut, ostr, _format)) {
                throw std::runtime_error("couldn't open ostream");
            }
        });
    }

    Output(std::ostream& ostr) {
        if (!open(*fileOut, ostr, format)) {
            throw std::runtime_error("couldn't open ostream");
        }
    }

    auto operator=(std::filesystem::path const& _path) -> Output& {
        fileOut.emplace(_path.c_str());
        return *this;
    }



};

}
