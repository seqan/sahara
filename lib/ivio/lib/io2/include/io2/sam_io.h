#pragma once

#include "sam_io_reader.h"
#include "sam_io_writer.h"

#include <vector>

namespace io2::sam_io {

class FileExtensions {
    std::vector<std::string> validExtensions {
        ".sam",
        ".bam",
    };

public:
    auto list() const {
        return validExtensions;
    }
    auto contains(std::filesystem::path const& _p) const {
        auto p = _p.string();
        for (auto ext : validExtensions) {
            if (p.size() < ext.size()) continue;
            if (p.substr(p.size() - ext.size()) == ext) return true;
        }
        return false;
    }
};

inline static FileExtensions fileExtensions{};

auto validExtension(std::filesystem::path const& p) {
    return fileExtensions.contains(p);
}

}
