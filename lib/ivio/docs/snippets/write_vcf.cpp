#include <ivio/ivio.h>

int main(int argc, char** argv) {
    if (argc != 2) {
        return 1;
    }
    auto file = std::filesystem::path{argv[1]};
    auto header = ivio::vcf::header { // fileformat is being set automatically
        .table = {
            {"fileDate", "20090805"},
            {"reference", "file:///somefile.fasta"},
        },
        .genotypes = {"NA00001", "NA00002", "NA00003" },
    };

    auto writer = ivio::vcf::writer{{file, header}};

    auto record = ivio::vcf::record {
        .chrom = "Name",
        .pos   = 1,
        .id    = "some id",
        .ref   = "A",
        .alt   = "G",
        .qual  = 29,
        .filters = "PASS",
        .infos   = "NS=3;DP=14;AF=0.5;DB;H2",
        .formats = "GT:GQ:DP:HQ",
        .sampmles = "0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.",
    }
    writer.write(record);
}
