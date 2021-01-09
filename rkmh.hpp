#ifndef MKMH_RKMH_HPP
#define MKMH_RKMH_HPP

#include <vector>
#include <string>
#include "mkmh.hpp"
#include <cmath>
#include <algorithm>

namespace rkmh {
    using namespace std;
    using namespace mkmh;

    inline void hash_sequences(std::vector<std::string *> &seqs,
                               std::vector<std::vector<mkmh::hash_t>> &hashes,
                               std::vector<int> &hash_lengths,
                               int kmer) {
        for (int i = 0; i < seqs.size(); ++i) {
            if (seqs[i] != nullptr) {
                hashes[i] = calc_hashes(seqs[i]->c_str(), seqs[i]->length(), kmer);
                std::sort(hashes[i].begin(), hashes[i].end());
                hash_lengths[i] = hashes[i].size();
            }
        }
    }

    inline std::vector<hash_t> hash_sequence(const char* seq, const int len, const int k) {
        std::vector<hash_t> hashes = calc_hashes(seq, len, k);
        std::sort(hashes.begin(), hashes.end());
        return hashes;
    }

    inline std::vector<hash_t> hash_sequence(const char* seq, const int len, const int k, const uint64_t sketch_size) {
        std::vector<hash_t> hashes = hash_sequence(seq, len, k);
        if (hashes.size() > sketch_size) {
            hashes.erase(hashes.begin() + sketch_size, hashes.end());
        }
        return hashes;
    }

    inline double compare(std::vector<mkmh::hash_t> alpha, std::vector<mkmh::hash_t> beta, int kmer_size) {
        int i = 0;
        int j = 0;

        uint64_t common = 0;
        uint64_t denom;

        while (i < alpha.size() && alpha[i] == 0) {
            i++;
        }
        while (j < beta.size() && beta[j] == 0) {
            j++;
        }
        denom = i + j;

        //todo early stopping
        while (i < alpha.size() && j < beta.size()) {
            if (alpha[i] == beta[j]) {
                i++;
                j++;
                common++;
            } else if (alpha[i] > beta[j]) {
                j++;
            } else {
                i++;
            }

            denom++;
        }

        // complete the union operation
        denom += alpha.size() - i;
        denom += beta.size() - j;

        //std::cerr << "common " << common << std::endl;
        //std::cerr << "denom " << denom << std::endl;

        double distance;

        //todo put a flag for denom: take the smallest between alpha.size, beta.size
        double jaccard = double(common) / denom;

        if (common == denom) // avoid -0
        {
            distance = 0;
        } else if (common == 0) // avoid inf
        {
            distance = 1.;
        } else {
            //distance = log(double(common + 1) / (denom + 1)) / log(1. / (denom + 1));
            distance = -log(2 * jaccard / (1. + jaccard)) / kmer_size;

            if (distance > 1) {
                distance = 1;
            }
        }

        return distance;
    }
}

#endif //MKMH_RKMH_HPP
