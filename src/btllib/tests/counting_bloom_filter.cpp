#include "../include/btllib/counting_bloom_filter.hpp"

#include "helpers.hpp"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <string>

int
main()
{
  std::cerr << "Testing CountingBloomFilter" << std::endl;
  btllib::CountingBloomFilter8 cbf(1024 * 1024, 3);

  cbf.insert({ 1, 10, 100 });
  cbf.insert({ 1, 10, 100 });
  cbf.insert({ 100, 200, 300 });

  assert(cbf.contains({ 1, 10, 100 }) == 2);
  assert(cbf.contains({ 100, 200, 300 }) == 1);
  assert(cbf.contains({ 1, 20, 100 }) == 0);

  auto filename = get_random_name(64);
  cbf.write(filename);

  btllib::CountingBloomFilter8 cbf2(filename);

  assert(cbf2.contains({ 1, 10, 100 }) == 2);
  assert(cbf2.contains({ 100, 200, 300 }) == 1);
  assert(cbf2.contains({ 1, 20, 100 }) == 0);

  std::remove(filename.c_str());

  std::string seq = "CACTATCGACGATCATTCGAGCATCAGCGACTG";
  std::string seq2 = "GTAGTACGATCAGCGACTATCGAGCTACGAGCA";
  assert(seq.size() == seq2.size());

  btllib::KmerCountingBloomFilter8 kmer_counting_bf(seq.size() / 2,
                                                    1024 * 1024);
  kmer_counting_bf.insert(seq);
  assert(kmer_counting_bf.contains(seq) == (seq.size() - seq.size() / 2 + 1));
  assert(kmer_counting_bf.contains(seq2) <= 1);

  return 0;
}