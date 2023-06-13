#pragma once

#include "seqan2_seqout_writer.h"
#include "fastq_reader.h"

namespace io2::fastq_io {

template <typename AlphabetS3 = seqan3::dna5,
          typename QualitiesS3 = seqan3::phred42>
using writer = seqan2_seqout_io::writer<AlphabetS3, QualitiesS3, ExtensionAndFormat>;

}
