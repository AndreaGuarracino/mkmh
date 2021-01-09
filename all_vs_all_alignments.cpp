#include <string>
#include <iostream>

#include "rkmh.hpp"

using namespace std;
using namespace mkmh;

int main(int argc, char **argv) {
    std::vector<std::string* > seqs;

    {
        std::ifstream in("data/sequences.fa");
        std::string line;
        do {
            std::getline(in, line);
            std::getline(in, line);
            seqs.push_back(new std::string(line));
        } while (in.good());
        in.close();
    }

    std::vector <std::vector<mkmh::hash_t>> hashes;
    hashes.resize(seqs.size());
    std::vector<int> hash_lengths;
    hash_lengths.resize(seqs.size());

    int kmer = 17;
    rkmh::hash_sequences(seqs, hashes, hash_lengths, kmer);

    std::cerr << "sequences: " << seqs.size() << std::endl;

    for (uint64_t i = 0; i < seqs.size(); ++i) {
        for (uint64_t j = i + 1; j < seqs.size(); ++j) {
            rkmh::compare(hashes[i], hashes[j], kmer);
        }
    }

    return 0;
}
