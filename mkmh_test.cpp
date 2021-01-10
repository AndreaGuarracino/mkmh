#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "mkmh.hpp"
#include "rkmh.hpp"
#include <fstream>
#include <string>
#include <algorithm>


using namespace std;
using namespace mkmh;

TEST_CASE("Reverse complement function works", "[reverse_complement]") {
    string t = "ACTGGCC";
    string rev = reverse_complement(t);

    SECTION("reverse_complement works on C++ strings.") {
        REQUIRE(rev == "GGCCAGT");
    }


    char k[6] = "AGGTC";
    char *ret = new char[6];
    char *retret = new char[6];
    reverse_complement(k, ret, 5);

    SECTION("reverse_complement returns expected string") {
        REQUIRE(strcmp(ret, "GACCT") == 0);
        REQUIRE(strlen(ret) == 5);
    }


    reverse_complement(ret, retret, 5);
    SECTION("reverse_complement does not affect its sequence pointer") {
        REQUIRE(ret == ret);
    }

    SECTION("reverse_complement, when applied twice, returns the original input string") {
        REQUIRE(*retret == *k);
    }


}

TEST_CASE("Canonical function works", "[canonical]") {
    string t = "ACTGGCNNNN";
    SECTION("canonical(string) catches Ns in a DNA string") {
        REQUIRE(canonical(t) == false);
    }

    char k[27] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    SECTION("canonical(char* , len) catches non-DNA letters") {
        REQUIRE(canonical(k, 26) == false);
    }

    char o[8] = "ACCCCTG";
    SECTION("canonical(char*, len) doesn't flag valid uppercase DNA characters") {
        REQUIRE(canonical(o, 7) == true);
    }


    char low[8] = "acccctg";
    SECTION("canonical(char*, len) doesn't flag valid lowercase DNA characters") {
        REQUIRE(canonical(low, 7) == true);
    }
}

TEST_CASE("Upper works for strings and chars") {
    char c[10] = "actgtgccc";
    char noncon[4] = "aBd";
    string d = "ABCDEFG";


}

TEST_CASE("v_set removes duplicates and returns a vector", "[v_set]") {

}

TEST_CASE("kmerize works as expected for strings", "[kmerize(string, ...)]") {

}

TEST_CASE("kmerize functions for char* work as expected", "[kmerize(char*, ...)]") {

}

TEST_CASE("minimizers behave as expected", "[minimizers]") {

}

TEST_CASE("Calc_hashes functions produce the right hashes", "[calc_hashes]") {

    char o[8] = "ACCCCTG";
    char t[8] = "ACCCCTG";

    SECTION("Hashes from calc_hashes for char* are consistent with those from calc_hash") {
        int numhashes = 7 - 4;
        hash_t *h = new hash_t[numhashes];
        calc_hashes((const char *) o, 4, h, numhashes);

        for (int i = 0; i < numhashes; i++) {
            cerr << string(o + i, o + i + 4) << " : " << *(h + i) << " : " << calc_hash(o + i, 4) << endl;
        }

        bool trip = false;
        for (int i = 0; i < 7 - 4; i++) {
            trip = *(h + i) != calc_hash(o + i, 4);
        }
        REQUIRE(trip == false);
    }

    SECTION("calc_hashes functions all return the same hashes") {
        vector <hash_t> x = calc_hashes((const char *) o, 7, 4);

        int num = 7 - 4;
        vector <hash_t> y(num);
        hash_t *h = y.data();
        calc_hashes((const char *) o, 4, h, num);

        string zstr(o);
        vector <hash_t> z = calc_hashes(zstr, 4);

        vector<int> kmers;
        kmers.push_back(4);
        vector <hash_t> multis = calc_hashes(zstr, kmers);

        vector <hash_t> matched = calc_hashes((const char *) t, 7, 4);

        REQUIRE(std::mismatch(x.begin(), x.end(), y.begin(), y.end()).first == x.end());
        REQUIRE(std::mismatch(x.begin(), x.end(), z.begin(), z.end()).first == x.end());
        REQUIRE(std::mismatch(y.begin(), y.end(), z.begin(), z.end()).first == y.end());
        REQUIRE(std::mismatch(x.begin(), x.end(), multis.begin(), multis.end()).first == x.end());
        REQUIRE(std::mismatch(x.begin(), x.end(), matched.begin(), matched.end()).first == x.end());
    }

    SECTION("calc_hashes works with multiple kmer sizes") {
        vector<int> kmers;
        kmers.push_back(4);
    }

    SECTION("calc_hashes returns identical hashes for forward and reverse compliment") {
        char s[8] = "AAAAAAA";
        char t[8] = "TTTTTTT";
        vector <hash_t> s_hashes = calc_hashes(s, 7, 4);
        vector <hash_t> t_hashes = calc_hashes(t, 7, 4);
        REQUIRE(std::mismatch(s_hashes.begin(), s_hashes.end(), t_hashes.begin(), t_hashes.end()).first ==
                s_hashes.end());
    }

}

TEST_CASE("Sort and minhashes return the expect", "sort / minhash") {
    string x = "ACTGGCTTGCC";
    string y = "GGCAAGCCAGT";


}

TEST_CASE("Calc_hash family of functions work correctly", "[calc_hash()]") {
    string x = "ACTGGCTTGCC";
    string y = "GGCAAGCCAGT";

    SECTION("Hashes of forward and reverse-complement sequences are equal") {
        hash_t c_x = calc_hash(x);
        hash_t c_y = calc_hash(y);
        REQUIRE(c_x == c_y);

        REQUIRE(calc_hash("AAAAAA") == calc_hash("TTTTTT"));
    }

    SECTION("Hashes of calc_hash and calc_hashes are equivalent") {
        vector <hash_t> x_hashes = calc_hashes(x, 10);
        vector <hash_t> comp_hashes;
        for (int i = 0; i < x.length() - 10; i++) {
            comp_hashes.push_back(calc_hash(x.substr(i, i + 10)));
        }

        REQUIRE(std::mismatch(x_hashes.begin(), x_hashes.end(), comp_hashes.begin(), comp_hashes.end()).first ==
                x_hashes.end());
    }

    SECTION("Non-canonical bases cause a sequence to hash to zero") {
        string z = "ACGTNTTA";
        REQUIRE(calc_hash(z) == 0);
    }
}

TEST_CASE("hash_intersection family of functions work correctly", "[hash_intersection]") {

    hash_t *x = new hash_t[4];
    hash_t *y = new hash_t[6];
    int num;

    x[0] = 0;
    x[1] = 2;
    x[2] = 20938475420;
    x[3] = 987728;

    y[0] = 0;
    y[1] = 1;
    y[2] = 0;
    y[3] = 20938475420;
    y[4] = 10;
    y[5] = 987728;

    SECTION("fastest hash-intersection works") {
        hash_intersection_size(x, 4, y, 6, num);
        REQUIRE(num == 2);
    }

}

TEST_CASE("kmer_to_integer", "[kmer_to_integer]") {
    char s[8] = "ATAGAAA";
    char s_p[8] = "ATAGAAA";
    char non_s[13] = "ATAGAATTTTAA";
    char fail[12] = "ATAGANNNNAA";

    hash_t *z = new hash_t[1];
    hash_t *z_p = new hash_t[1];
    hash_t x;

    bool r = kmer_to_integer(s, 7, z[0]);
    kmer_to_integer(s_p, 7, z_p[0]);
    bool shouldfail = kmer_to_integer(fail, 11, x);
    hash_t nonz;

    kmer_to_integer(non_s, 12, nonz);
    REQUIRE(shouldfail == false);
    REQUIRE(x == 0);
    //cout << (hash_t) z << endl;
    REQUIRE(r == true);
    REQUIRE(z[0] == z_p[0]);
    REQUIRE(nonz != z[0]);

    delete[] z;
    delete[] z_p;
}


TEST_CASE("dirty_seq_to_hahes", "[dirty_seq_to_hahes]") {
    std::string dirty_seq = "ATTCCAXXXTTCCAATCGATTCCATTCCATTCAATTCCATTXXXCCACTCGGGTTGATTCCATTCCATTCCATTCCATTTAATTCCATTCCACTCTGGTTGATTCCATTCCATTCCATTCCA";

    std::vector <hash_t> hashes = rkmh::hash_sequence(dirty_seq.c_str(), dirty_seq.length(), 17);

    std::vector<mkmh::hash_t> hashes_to_check { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 149891068869527339, 168382979207254341, 265041207590036277, 396750732942071689, 396750732942071689, 462790877938481053, 747973565502263643, 845564525554761269, 874009950677834167, 947680952221525423, 1021111316835440685, 1021111316835440685, 1060677449931109975, 1118486967359146196, 1189674342521077837, 1472582835660144245, 1525396889846938985, 1528643647800622190, 1528643647800622190, 1569623208124008282, 1761348516187543570, 1761348516187543570, 1761348516187543570, 2377464206934081074, 2515733927996799239, 2764378512359959862, 2807469306412738643, 2825702366505580794, 3025572705888784997, 3249998615155623264, 3406714372864218566, 3479694609920169651, 3481951538478207575, 3739868762558139823, 3761251575281293195, 3946398098715183825, 4329170880601917952, 4456068341283124241, 4474381755821901961, 4904900660178402163, 5182246861356579551, 5182246861356579551, 5207041342950107931, 5407824226294554335, 5730288989241871124, 5730288989241871124, 5783832666340056974, 5786794638415951585, 6715540742195326870, 6888513602129946462, 6888513602129946462, 8020997838156559485, 8141620841386204796, 8381821265184171007, 8638352122713257417, 8811304748114462666, 8811304748114462666, 9005894505542228239, 9336796222033361739, 9590797608944518166, 9684358218647197431, 9821931412729383740, 10079800293474304882, 10282242885806663493, 10734597772072675339, 10812385764316197214, 11087805979949599047, 12052992867156394409, 12802981501696578844, 12931068024011070043, 12947430219925530077, 13962010917989428056, 14413852054642462893, 14784589235147396869, 15519400483688987475, 15519400483688987475, 15519400483688987475 };

    SECTION("The graph is as expected when sorted") {
        REQUIRE(hashes.size() == hashes_to_check.size());

        for (uint64_t i = 0; i < hashes.size(); ++i) {
            REQUIRE(hashes[i] == hashes_to_check[i]);
        }
    }
}
