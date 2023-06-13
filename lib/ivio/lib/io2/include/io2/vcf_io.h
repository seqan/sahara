#pragma once

#include "vcf_io_reader.h"
#include "vcf_io_writer.h"

#include <vector>

namespace io2::vcf_io {

class FileExtensions {
    std::vector<std::string> validExtensions {
        ".vcf",
        ".bcf",
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
