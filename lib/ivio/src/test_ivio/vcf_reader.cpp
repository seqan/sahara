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

TEST_CASE("reading vcf files", "[vcf][reader]") {
    auto tmp = std::filesystem::temp_directory_path() / "ivio_test";
    std::filesystem::create_directory(tmp);

    auto expected = std::vector<ivio::vcf::record> {
        ivio::vcf::record{.chrom = "20", .pos = 14370,   .id = "rs6054257", .ref = "G",   .alts = "A",      .qual = 29.f, .filters = "PASS", .infos = "NS=3;DP=14;AF=0.5;DB;H2",           .formats = "GT:GQ:DP:HQ", .samples = "0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,."},
        ivio::vcf::record{.chrom = "20", .pos = 17330,   .id = ".",         .ref = "T",   .alts = "A",      .qual = 3.f,  .filters = "q10",  .infos = "NS=3;DP=11;AF=0.017",               .formats = "GT:GQ:DP:HQ", .samples = "0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3"},
        ivio::vcf::record{.chrom = "20", .pos = 1110696, .id = "rs6040355", .ref = "A",   .alts = "G,T",    .qual = 67.f, .filters = "PASS", .infos = "NS=2;DP=10;AF=0.333,0.667;AA=T;DB", .formats = "GT:GQ:DP:HQ", .samples = "1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4"},
        ivio::vcf::record{.chrom = "20", .pos = 1230237, .id = ".",         .ref = "T",   .alts = ".",      .qual = 47.f, .filters = "PASS", .infos = "NS=3;DP=13;AA=T",                   .formats = "GT:GQ:DP:HQ", .samples = "0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2"},
        ivio::vcf::record{.chrom = "20", .pos = 1234567, .id = "microsat1", .ref = "GTC", .alts = "G,GTCT", .qual = 50.f, .filters = "PASS", .infos = "NS=3;DP=9;AA=G",                    .formats = "GT:GQ:DP",    .samples = "0/1:35:4	0/2:17:2	1/1:40:3"},
    };

    auto expected_header = ivio::vcf::header {
        .table = {
            {R"(fileformat)", R"(VCFv4.3)"},
            {R"(fileDate)", R"(20090805)"},
            {R"(source)", R"(myImputationProgramV3.1)"},
            {R"(reference)", R"(file:///seq/references/1000GenomesPilot-NCBI36.fasta)"},
            {R"(contig)", R"(<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>)"},
            {R"(phasing)", R"(partial)"},
            {R"(INFO)", R"(<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">)"},
            {R"(INFO)", R"(<ID=DP,Number=1,Type=Integer,Description="Total Depth">)"},
            {R"(INFO)", R"(<ID=AF,Number=A,Type=Float,Description="Allele Frequency">)"},
            {R"(INFO)", R"(<ID=AA,Number=1,Type=String,Description="Ancestral Allele">)"},
            {R"(INFO)", R"(<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">)"},
            {R"(INFO)", R"(<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">)"},
            {R"(FILTER)", R"(<ID=q10,Description="Quality below 10">)"},
            {R"(FILTER)", R"(<ID=s50,Description="Less than 50% of samples have data">)"},
            {R"(FORMAT)", R"(<ID=GT,Number=1,Type=String,Description="Genotype">)"},
            {R"(FORMAT)", R"(<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">)"},
            {R"(FORMAT)", R"(<ID=DP,Number=1,Type=Integer,Description="Read Depth">)"},
            {R"(FORMAT)", R"(<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">)"},
        },
        .genotypes = {"NA00001", "NA00002", "NA00003"},
    };

    auto test_data = std::string {
        R"(##fileformat=VCFv4.3)""\n"
        R"(##fileDate=20090805)""\n"
        R"(##source=myImputationProgramV3.1)""\n"
        R"(##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta)""\n"
        R"(##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>)""\n"
        R"(##phasing=partial)""\n"
        R"(##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">)""\n"
        R"(##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">)""\n"
        R"(##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">)""\n"
        R"(##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">)""\n"
        R"(##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">)""\n"
        R"(##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">)""\n"
        R"(##FILTER=<ID=q10,Description="Quality below 10">)""\n"
        R"(##FILTER=<ID=s50,Description="Less than 50% of samples have data">)""\n"
        R"(##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">)""\n"
        R"(##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">)""\n"
        R"(##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">)""\n"
        R"(##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">)""\n"
        R"(#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003)""\n"
        R"(20	14370	rs6054257	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2	GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.)""\n"
        R"(20	17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3)""\n"
        R"(20	1110696	rs6040355	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4)""\n"
        R"(20	1230237	.	T	.	47	PASS	NS=3;DP=13;AA=T	GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2)""\n"
        R"(20	1234567	microsat1	GTC	G,GTCT	50	PASS	NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3)""\n"
    };

    SECTION("Creating Test file") {
        auto ofs = std::ofstream{tmp / "file.vcf", std::ios::binary};
        ofs << test_data;
    }

    SECTION("Read from std::filesystem::path") {
        auto reader = ivio::vcf::reader{{tmp / "file.vcf"}};
        CHECK(reader.header().table == expected_header.table);
        CHECK(reader.header().genotypes == expected_header.genotypes);
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(vec == expected);
    }

    SECTION("Read from std::ifstream") {
        auto ifs = std::ifstream{tmp / "file.vcf"};
        auto reader = ivio::vcf::reader{{ifs}};
        CHECK(reader.header().table == expected_header.table);
        CHECK(reader.header().genotypes == expected_header.genotypes);
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(vec == expected);
    }

    SECTION("Read from std::stringstream") {
        auto ss = std::stringstream{test_data};
        auto reader = ivio::vcf::reader{{ss}};
        CHECK(reader.header().table == expected_header.table);
        CHECK(reader.header().genotypes == expected_header.genotypes);
        auto vec = std::vector(begin(reader), end(reader));
        static_assert(std::same_as<decltype(vec), decltype(expected)>, "vec and expected should have the exact same type");
        CHECK(vec == expected);
    }

    SECTION("cleanup - deleting temp folder") {
        std::filesystem::remove_all(tmp);
    }
}
