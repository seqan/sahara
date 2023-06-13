// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include <catch2/catch.hpp>
#include <filesystem>
#include <fstream>
#include <ivio/ivio.h>

TEST_CASE("reading sam files", "[sam][reader]") {
    auto tmp = std::filesystem::temp_directory_path() / "ivio_test";
    std::filesystem::create_directory(tmp);

    auto expected = std::vector<ivio::sam::record> {
        { .qname = "D00360:64:HBAP3ADXX:2:1115:16968:20035", .flag = 2147, .rname = "chrM", .pos = 1, .mapq = 6,  .cigar = "61H40M",  .rnext = "=", .pnext = 8953,  .tlen = 9053,  .seq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCT",                                                                                       .qual = "IGGDGGEHEA=BDFFFFCECEECE>CCCBCCCCC=80>B<",                                                                                      .tags = "NM:i:0	MD:Z:40	AS:i:40	XS:i:0	SA:Z:chrM,16511,+,61M40S,6,1;"},
        { .qname = "D00360:72:C6UP4ANXX:1:1115:17362:52449", .flag = 161,  .rname = "chrM", .pos = 1, .mapq = 60, .cigar = "14S111M", .rnext = "=", .pnext = 10470, .tlen = 10594, .seq = "AAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCA",  .qual = "CBBBBGGGGGGGGGD>CGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG//<<GGEFCFGGGGCEEG>FGGFFFFABGE:1:EAGGGGEE/@DGGGG..EGGBED>CAA..C>C@GBD..989C.", .tags = "NM:i:1	MD:Z:72G38	AS:i:106	XS:i:70"},
        { .qname = "D00360:72:C6UP4ANXX:1:1304:9139:2419",   .flag = 161,  .rname = "chrM", .pos = 1, .mapq = 32, .cigar = "26S99M",  .rnext = "=", .pnext = 9579,  .tlen = 9703,  .seq = "TTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCT",  .qual = "BCBBBGGGGGGGE>@GGEGGGGGGGGGGGGGGGGGGGGCFEGGEGGGGGGGGGGGGGGG:FGGGGBGGGGGGGGGGDGG>>EFCDF@GGDGGGGGGGC?EGGGBGGDGGGGGGGGGGGGDGGGGG", .tags = "NM:i:1	MD:Z:72G26	AS:i:94	XS:i:81"},
        { .qname = "D00360:71:C61DEANXX:1:2303:17550:46600", .flag = 161,  .rname = "chrM", .pos = 1, .mapq = 60, .cigar = "13S112M", .rnext = "=", .pnext = 9323,  .tlen = 9420,  .seq = "AGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCAC",  .qual = "3?AABGEGGB/EFGGGE1B<FGGFGEGGEGEFGGGGGGGEGBGFBGG<//>CCFGGG>G>F1FGGGGGGEGGGGGFGG01=ED;6CGBBEGGDGGGD.DDB@C@CGCCGGGGCD@GEGGG@C8DG", .tags = "NM:i:1	MD:Z:72G39	AS:i:107	XS:i:69"},
        { .qname = "D00360:71:C61DEANXX:1:2304:1384:69650",  .flag = 73,   .rname = "chrM", .pos = 1, .mapq = 60, .cigar = "8S107M",  .rnext = "=", .pnext = 1,     .tlen = 0,     .seq = "TCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGGATGCACGCGATAGCATTGCGAGACGCTGGAGCCGG",            .qual = "AABA=/CGBDFGGEG>0;=F>FBFDGGGG@FFGGEGGGGGGBCCGGG>FFF>C>1>1F1:?1CF>1EDFGGGGD1:99?/9C.://89.9C>8/7C@@D8>..6.6.888.:CC6",           .tags = "NM:i:2	MD:Z:71T0G34	AS:i:97	XS:i:0"},
        { .qname = "D00360:71:C61DEANXX:1:2304:1384:69650",  .flag = 133,  .rname = "chrM", .pos = 1, .mapq = 0,  .cigar = "*",       .rnext = "=", .pnext = 1,     .tlen = 0,     .seq = "TTATTTTTCCTTTGAGTGTTTCAATCGGTAGTACTACTCTATTGTCTTCGTTATGGTGACACATAGCTTGTTGAGTTTTGTTTAATTGGTGAAGATTGTATGAGTTTAAGATTAGTTTACCTTAT",  .qual = "3=A=B1@111111110==11=111@10//>/11:1>1EB1EG11E1111011>011111>1>E1:1>111111111111:11?FGE:0000>0=0:00=00:;0000=000=0000000000?@0", .tags = "AS:i:0	XS:i:0"},
        { .qname = "D00360:64:HBAP3ADXX:1:1103:13477:91844", .flag = 2129, .rname = "chrM", .pos = 1, .mapq = 26, .cigar = "64H37M",  .rnext = "=", .pnext = 8699,  .tlen = 8663,  .seq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGA",                                                                                          .qual = "D:::1AEB<IFCF=BA:ED@E>CA:DDC2+=:DD;=?",                                                                                         .tags = "NM:i:0	MD:Z:37	AS:i:37	XS:i:0	SA:Z:chrM,16508,-,64M37S,26,1;"},
        { .qname = "D00360:64:HBAP3ADXX:1:1111:9203:55170",  .flag = 2193, .rname = "chrM", .pos = 1, .mapq = 60, .cigar = "70H31M",  .rnext = "=", .pnext = 6094,  .tlen = 6064,  .seq = "GATCACAGGTCTATCACCCTATTAACCACTC",                                                                                                .qual = "EEIFEFFEC<A<:E?EEEDDDDD7DA44?8=",                                                                                               .tags = "NM:i:0	MD:Z:31	AS:i:31	XS:i:0	SA:Z:chrM,16502,-,70M31S,60,2;"},
        { .qname = "D00360:64:HBAP3ADXX:1:1116:10699:77936", .flag = 147,  .rname = "chrM", .pos = 1, .mapq = 60, .cigar = "49M",     .rnext = "=", .pnext = 6526,  .tlen = 6478,  .seq = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCA",                                                                              .qual = "IJIGIJJJJJJJJIHGJJJJJJJIHHHCIIIHJJJJHHHHHFFFFFCCC",                                                                             .tags = "NM:i:0	MD:Z:49	AS:i:49	XS:i:0"},
    };

    auto expected_header = std::vector<std::string> {
        std::string{"@HD	VN:1.3	SO:coordinate"},
        std::string{"@SQ	SN:chrM	LN:16571"},
        std::string{"@SQ	SN:chr1	LN:249250621"},
        std::string{"@SQ	SN:chr2	LN:243199373"},
        std::string{"@SQ	SN:chr3	LN:198022430"},
        std::string{"@SQ	SN:chr4	LN:191154276"},
        std::string{"@SQ	SN:chr5	LN:180915260"},
        std::string{"@SQ	SN:chr6	LN:171115067"},
        std::string{"@SQ	SN:chr7	LN:159138663"},
        std::string{"@SQ	SN:chr8	LN:146364022"},
        std::string{"@SQ	SN:chr9	LN:141213431"},
        std::string{"@SQ	SN:chr10	LN:135534747"},
        std::string{"@SQ	SN:chr11	LN:135006516"},
        std::string{"@SQ	SN:chr12	LN:133851895"},
        std::string{"@SQ	SN:chr13	LN:115169878"},
        std::string{"@SQ	SN:chr14	LN:107349540"},
        std::string{"@SQ	SN:chr15	LN:102531392"},
        std::string{"@SQ	SN:chr16	LN:90354753"},
        std::string{"@SQ	SN:chr17	LN:81195210"},
        std::string{"@SQ	SN:chr18	LN:78077248"},
        std::string{"@SQ	SN:chr19	LN:59128983"},
        std::string{"@SQ	SN:chr20	LN:63025520"},
        std::string{"@SQ	SN:chr21	LN:48129895"},
        std::string{"@SQ	SN:chr22	LN:51304566"},
        std::string{"@SQ	SN:chrX	LN:155270560"},
        std::string{"@SQ	SN:chrY	LN:59373566"},
        std::string{"@PG	ID:bwa	PN:bwa	VN:0.7.7-r441	CL:bwa mem -t 2 ref.fa left.fq.gz right.fq.gz"},
        std::string{"@PG	ID:samtools	PN:samtools	PP:bwa	VN:1.15.1	CL:samtools view -s 0.01 -b HG002.mate_pair.sorted.bam"},
        std::string{"@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.15.1	CL:samtools view -h sampled.bam"},
    };


    auto test_data = std::string {
        "@HD	VN:1.3	SO:coordinate\n"
        "@SQ	SN:chrM	LN:16571\n"
        "@SQ	SN:chr1	LN:249250621\n"
        "@SQ	SN:chr2	LN:243199373\n"
        "@SQ	SN:chr3	LN:198022430\n"
        "@SQ	SN:chr4	LN:191154276\n"
        "@SQ	SN:chr5	LN:180915260\n"
        "@SQ	SN:chr6	LN:171115067\n"
        "@SQ	SN:chr7	LN:159138663\n"
        "@SQ	SN:chr8	LN:146364022\n"
        "@SQ	SN:chr9	LN:141213431\n"
        "@SQ	SN:chr10	LN:135534747\n"
        "@SQ	SN:chr11	LN:135006516\n"
        "@SQ	SN:chr12	LN:133851895\n"
        "@SQ	SN:chr13	LN:115169878\n"
        "@SQ	SN:chr14	LN:107349540\n"
        "@SQ	SN:chr15	LN:102531392\n"
        "@SQ	SN:chr16	LN:90354753\n"
        "@SQ	SN:chr17	LN:81195210\n"
        "@SQ	SN:chr18	LN:78077248\n"
        "@SQ	SN:chr19	LN:59128983\n"
        "@SQ	SN:chr20	LN:63025520\n"
        "@SQ	SN:chr21	LN:48129895\n"
        "@SQ	SN:chr22	LN:51304566\n"
        "@SQ	SN:chrX	LN:155270560\n"
        "@SQ	SN:chrY	LN:59373566\n"
        "@PG	ID:bwa	PN:bwa	VN:0.7.7-r441	CL:bwa mem -t 2 ref.fa left.fq.gz right.fq.gz\n"
        "@PG	ID:samtools	PN:samtools	PP:bwa	VN:1.15.1	CL:samtools view -s 0.01 -b HG002.mate_pair.sorted.bam\n"
        "@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.15.1	CL:samtools view -h sampled.bam\n"
        "D00360:64:HBAP3ADXX:2:1115:16968:20035	2147	chrM	1	6	61H40M	=	8953	9053	GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCT	IGGDGGEHEA=BDFFFFCECEECE>CCCBCCCCC=80>B<	NM:i:0	MD:Z:40	AS:i:40	XS:i:0	SA:Z:chrM,16511,+,61M40S,6,1;\n"
        "D00360:72:C6UP4ANXX:1:1115:17362:52449	161	chrM	1	60	14S111M	=	10470	10594	AAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCA	CBBBBGGGGGGGGGD>CGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG//<<GGEFCFGGGGCEEG>FGGFFFFABGE:1:EAGGGGEE/@DGGGG..EGGBED>CAA..C>C@GBD..989C.	NM:i:1	MD:Z:72G38	AS:i:106	XS:i:70\n"
        "D00360:72:C6UP4ANXX:1:1304:9139:2419	161	chrM	1	32	26S99M	=	9579	9703	TTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCT	BCBBBGGGGGGGE>@GGEGGGGGGGGGGGGGGGGGGGGCFEGGEGGGGGGGGGGGGGGG:FGGGGBGGGGGGGGGGDGG>>EFCDF@GGDGGGGGGGC?EGGGBGGDGGGGGGGGGGGGDGGGGG	NM:i:1	MD:Z:72G26	AS:i:94	XS:i:81\n"
        "D00360:71:C61DEANXX:1:2303:17550:46600	161	chrM	1	60	13S112M	=	9323	9420	AGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCAC	3?AABGEGGB/EFGGGE1B<FGGFGEGGEGEFGGGGGGGEGBGFBGG<//>CCFGGG>G>F1FGGGGGGEGGGGGFGG01=ED;6CGBBEGGDGGGD.DDB@C@CGCCGGGGCD@GEGGG@C8DG	NM:i:1	MD:Z:72G39	AS:i:107	XS:i:69\n"
        "D00360:71:C61DEANXX:1:2304:1384:69650	73	chrM	1	60	8S107M	=	1	0	TCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGGATGCACGCGATAGCATTGCGAGACGCTGGAGCCGG	AABA=/CGBDFGGEG>0;=F>FBFDGGGG@FFGGEGGGGGGBCCGGG>FFF>C>1>1F1:?1CF>1EDFGGGGD1:99?/9C.://89.9C>8/7C@@D8>..6.6.888.:CC6	NM:i:2	MD:Z:71T0G34	AS:i:97	XS:i:0\n"
        "D00360:71:C61DEANXX:1:2304:1384:69650	133	chrM	1	0	*	=	1	0	TTATTTTTCCTTTGAGTGTTTCAATCGGTAGTACTACTCTATTGTCTTCGTTATGGTGACACATAGCTTGTTGAGTTTTGTTTAATTGGTGAAGATTGTATGAGTTTAAGATTAGTTTACCTTAT	3=A=B1@111111110==11=111@10//>/11:1>1EB1EG11E1111011>011111>1>E1:1>111111111111:11?FGE:0000>0=0:00=00:;0000=000=0000000000?@0	AS:i:0	XS:i:0\n"
        "D00360:64:HBAP3ADXX:1:1103:13477:91844	2129	chrM	1	26	64H37M	=	8699	8663	GATCACAGGTCTATCACCCTATTAACCACTCACGGGA	D:::1AEB<IFCF=BA:ED@E>CA:DDC2+=:DD;=?	NM:i:0	MD:Z:37	AS:i:37	XS:i:0	SA:Z:chrM,16508,-,64M37S,26,1;\n"
        "D00360:64:HBAP3ADXX:1:1111:9203:55170	2193	chrM	1	60	70H31M	=	6094	6064	GATCACAGGTCTATCACCCTATTAACCACTC	EEIFEFFEC<A<:E?EEEDDDDD7DA44?8=	NM:i:0	MD:Z:31	AS:i:31	XS:i:0	SA:Z:chrM,16502,-,70M31S,60,2;\n"
        "D00360:64:HBAP3ADXX:1:1116:10699:77936	147	chrM	1	60	49M	=	6526	6478	GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCA	IJIGIJJJJJJJJIHGJJJJJJJIHHHCIIIHJJJJHHHHHFFFFFCCC	NM:i:0	MD:Z:49	AS:i:49	XS:i:0\n"
    };

    SECTION("Creating Test file") {
        auto ofs = std::ofstream{tmp / "file.sam", std::ios::binary};
        ofs << test_data;
    }

    SECTION("Read from std::filesystem::path") {
        auto reader = ivio::sam::reader{{tmp / "file.sam"}};
        CHECK(reader.header == expected_header);
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(vec == expected);
    }

    SECTION("Read from std::ifstream") {
        auto ifs = std::ifstream{tmp / "file.sam"};
        auto reader = ivio::sam::reader{{ifs}};
        CHECK(reader.header == expected_header);
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(vec == expected);
    }

    SECTION("Read from std::stringstream") {
        auto ss = std::stringstream{test_data};
        auto reader = ivio::sam::reader{{ss}};
        CHECK(reader.header == expected_header);
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(vec == expected);
    }

    SECTION("cleanup - deleting temp folder") {
        std::filesystem::remove_all(tmp);
    }
}

TEST_CASE("reading invalid sam files", "[sam][reader][invalid]") {
    auto tmp = std::filesystem::temp_directory_path() / "ivio_test";
    std::filesystem::create_directory(tmp);

    auto test_data = std::string {
        "@HD	VN:1.3	SO:coordinate\n"
        "@SQ	SN:chrM	LN:16571\n"
        "@SQ	SN:chr1	LN:249250621\n"
        "@SQ	SN:chr2	LN:243199373\n"
        "@SQ	SN:chr3	LN:198022430\n"
        "@SQ	SN:chr4	LN:191154276\n"
        "@SQ	SN:chr5	LN:180915260\n"
        "@SQ	SN:chr6	LN:171115067\n"
        "@SQ	SN:chr7	LN:159138663\n"
        "@SQ	SN:chr8	LN:146364022\n"
        "@SQ	SN:chr9	LN:141213431\n"
        "@SQ	SN:chr10	LN:135534747\n"
        "@SQ	SN:chr11	LN:135006516\n"
        "@SQ	SN:chr12	LN:133851895\n"
        "@SQ	SN:chr13	LN:115169878\n"
        "@SQ	SN:chr14	LN:107349540\n"
        "@SQ	SN:chr15	LN:102531392\n"
        "@SQ	SN:chr16	LN:90354753\n"
        "@SQ	SN:chr17	LN:81195210\n"
        "@SQ	SN:chr18	LN:78077248\n"
        "@SQ	SN:chr19	LN:59128983\n"
        "@SQ	SN:chr20	LN:63025520\n"
        "@SQ	SN:chr21	LN:48129895\n"
        "@SQ	SN:chr22	LN:51304566\n"
        "@SQ	SN:chrX	LN:155270560\n"
        "@SQ	SN:chrY	LN:59373566\n"
        "@PG	ID:bwa	PN:bwa	VN:0.7.7-r441	CL:bwa mem -t 2 ref.fa left.fq.gz right.fq.gz\n"
        "@PG	ID:samtools	PN:samtools	PP:bwa	VN:1.15.1	CL:samtools view -s 0.01 -b HG002.mate_pair.sorted.bam\n"
        "@PG	ID:samtools.1	PN:samtools	PP:samtools	VN:1.15.1	CL:samtools view -h sampled.bam\n"
        "D00360:64:HBAP3ADXX:2:1115:16968:20035	invalidPosition	chrM	1	6	61H40M	=	8953	9053	GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCT	IGGDGGEHEA=BDFFFFCECEECE>CCCBCCCCC=80>B<	NM:i:0	MD:Z:40	AS:i:40	XS:i:0	SA:Z:chrM,16511,+,61M40S,6,1;\n"
    };

    SECTION("Creating Test file") {
        auto ofs = std::ofstream{tmp / "file.sam", std::ios::binary};
        ofs << test_data;
    }

    SECTION("Read from std::filesystem::path") {
        auto reader = ivio::sam::reader{{tmp / "file.sam"}};
        REQUIRE_THROWS([&]() {
            auto vec = std::vector(begin(reader), end(reader));
        }());
    }

    SECTION("cleanup - deleting temp folder") {
        std::filesystem::remove_all(tmp);
    }
}
