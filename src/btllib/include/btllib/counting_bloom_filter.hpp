#ifndef BTLLIB_BLOOM_FILTER_HPP
#define BTLLIB_BLOOM_FILTER_HPP

#include "bloom_filter.hpp"
#include "nthash.hpp"
#include "status.hpp"

#include "vendor/cpptoml.hpp"

#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace btllib {

static const char* const COUNTING_BLOOM_FILTER_MAGIC_HEADER =
  "BTLCountingBloomFilter_v2";

template<typename T>
class CountingBloomFilter
{

public:
  CountingBloomFilter() {}
  CountingBloomFilter(size_t bytes, unsigned hash_num);
  CountingBloomFilter(const std::string& path);

  void insert(const std::vector<uint64_t>& hashes);
  void insert(const uint64_t* hashes);

  T contains(const std::vector<uint64_t>& hashes) const;
  T contains(const uint64_t* hashes) const;

  size_t get_bytes() const { return bytes; }
  uint64_t get_pop_cnt() const;
  unsigned get_hash_num() const { return hash_num; }
  double get_fpr() const;

  /**
   * Write bloom filter data to a file
   * @param path output filepath
   */
  void write(const std::string& path);

private:
  std::vector<unsigned char> bytearray;
  size_t bytes = 0;
  size_t counters = 0;
  unsigned hash_num = 0;
};

template<typename T>
class KmerCountingBloomFilter : public CountingBloomFilter<T>
{

public:
  KmerCountingBloomFilter(unsigned k, size_t bytes, unsigned hash_num = 4);

  void insert(const std::string& seq);
  void insert(const char* seq, size_t seq_len);

  uint64_t contains(const std::string& seq) const;
  uint64_t contains(const char* seq, size_t seq_len) const;

private:
  unsigned k;
};

using CountingBloomFilter8 = CountingBloomFilter<uint8_t>;
using CountingBloomFilter16 = CountingBloomFilter<uint16_t>;
using CountingBloomFilter32 = CountingBloomFilter<uint32_t>;

using KmerCountingBloomFilter8 = KmerCountingBloomFilter<uint8_t>;
using KmerCountingBloomFilter16 = KmerCountingBloomFilter<uint16_t>;
using KmerCountingBloomFilter32 = KmerCountingBloomFilter<uint32_t>;

template<typename T>
inline CountingBloomFilter<T>::CountingBloomFilter(size_t bytes,
                                                   unsigned hash_num)
  : bytes(std::ceil(bytes / sizeof(uint64_t)) * sizeof(uint64_t))
  , counters(bytes / sizeof(T))
  , hash_num(hash_num)
{
  bytearray.resize(bytes);
}

template<typename T>
inline CountingBloomFilter<T>::CountingBloomFilter(const std::string& path)
{
  std::ifstream file(path);

  std::string magic_with_brackets =
    std::string("[") + COUNTING_BLOOM_FILTER_MAGIC_HEADER + "]";

  std::string line;
  std::getline(file, line);
  if (line != magic_with_brackets) {
    log_error(
      std::string("Magic string does not match (likely version mismatch)\n") +
      "Your magic string:\t" + line + "\n" + "BloomFilter magic string:\t" +
      magic_with_brackets);
    std::exit(EXIT_FAILURE);
  }

  /* Read bloom filter line by line until it sees "[HeaderEnd]"
  which is used to mark the end of the header section and
  assigns the header to a char array*/
  std::string toml_buffer(line + '\n');
  bool header_end_found = false;
  while (bool(std::getline(file, line))) {
    toml_buffer.append(line + '\n');
    if (line == "[HeaderEnd]") {
      header_end_found = true;
      break;
    }
  }
  if (!header_end_found) {
    log_error("Pre-built bloom filter does not have the correct header end.");
    std::exit(EXIT_FAILURE);
  }

  // Send the char array to a stringstream for the cpptoml parser to parse
  std::istringstream toml_stream(toml_buffer);
  cpptoml::parser toml_parser(toml_stream);
  auto header_config = toml_parser.parse();

  // Obtain header values from toml parser and assign them to class members
  auto table = header_config->get_table(COUNTING_BLOOM_FILTER_MAGIC_HEADER);
  bytes = *table->get_as<size_t>("bytes");
  hash_num = *table->get_as<unsigned>("hash_num");
  counters = bytes / sizeof(T);
  check_error(sizeof(T) * CHAR_BIT != *table->get_as<size_t>("counter_bits"),
              "CountingBloomFilter" + std::to_string(sizeof(T) * CHAR_BIT) +
                " tried to load a file of CountingBloomFilter" +
                std::to_string(*table->get_as<size_t>("counter_bits")));

  bytearray.resize(bytes);
  file.read((char*)bytearray.data(), bytes);
}

template<typename T>
inline void
CountingBloomFilter<T>::insert(const std::vector<uint64_t>& hashes)
{
  insert(hashes.data());
}

template<typename T>
inline void
CountingBloomFilter<T>::insert(const uint64_t* hashes)
{
  // Update flag to track if increment is done on at least one counter
  bool update_done = false;
  T new_val;
  T min_val = contains(hashes);
  while (!update_done) {
    // Simple check to deal with overflow
    new_val = min_val + 1;
    if (min_val > new_val) {
      return;
    }
    for (size_t i = 0; i < hash_num; ++i) {
      if (__sync_bool_compare_and_swap( // NOLINT
            &(((T*)(bytearray.data()))[hashes[i] % counters]),
            min_val,
            new_val)) { // NOLINT
        update_done = true;
      }
    }
    // Recalculate minval because if increment fails, it needs a new minval to
    // use and if it doesnt hava a new one, the while loop runs forever.
    if (!update_done) {
      min_val = contains(hashes);
    }
  }
}

template<typename T>
inline T
CountingBloomFilter<T>::contains(const std::vector<uint64_t>& hashes) const
{
  return contains(hashes.data());
}

template<typename T>
inline T
CountingBloomFilter<T>::contains(const uint64_t* hashes) const
{
  T min = ((T*)(bytearray.data()))[hashes[0] % counters];
  for (size_t i = 1; i < hash_num; ++i) {
    size_t idx = hashes[i] % counters;
    if (((T*)(bytearray.data()))[idx] < min) {
      min = ((T*)(bytearray.data()))[idx];
    }
  }
  return min;
}

template<typename T>
inline uint64_t
CountingBloomFilter<T>::get_pop_cnt() const
{
  uint64_t pop_cnt = 0;
#pragma omp parallel for reduction(+ : pop_cnt)
  for (size_t i = 0; i < counters; ++i) {
    if (((T*)(bytearray.data()))[i]) {
      ++pop_cnt;
    }
  }
  return pop_cnt;
}

template<typename T>
inline double
CountingBloomFilter<T>::get_fpr() const
{
  return std::pow(double(get_pop_cnt()) / double(bytes), double(hash_num));
}

template<typename T>
inline void
CountingBloomFilter<T>::write(const std::string& path)
{
  std::ofstream file(path.c_str(), std::ios::out | std::ios::binary);

  /* Initialize cpptoml root table
    Note: Tables and fields are unordered
    Ordering of table is maintained by directing the table
    to the output stream immediately after completion  */
  auto root = cpptoml::make_table();

  /* Initialize bloom filter section and insert fields
      and output to ostream */
  auto header = cpptoml::make_table();
  header->insert("bytes", bytes);
  header->insert("hash_num", hash_num);
  header->insert("counter_bits", size_t(sizeof(T) * CHAR_BIT));
  root->insert(COUNTING_BLOOM_FILTER_MAGIC_HEADER, header);
  file << *root << "[HeaderEnd]\n";

  file.write((char*)bytearray.data(), bytes);
}

template<typename T>
inline KmerCountingBloomFilter<T>::KmerCountingBloomFilter(unsigned k,
                                                           size_t bytes,
                                                           unsigned hash_num)
  : CountingBloomFilter<T>(bytes, hash_num)
  , k(k)
{}

template<typename T>
inline void
KmerCountingBloomFilter<T>::insert(const std::string& seq)
{
  insert(seq.c_str(), seq.size());
}

template<typename T>
inline void
KmerCountingBloomFilter<T>::insert(const char* seq, size_t seq_len)
{
  NtHash nthash(seq, seq_len, k, CountingBloomFilter<T>::get_hash_num());
  while (nthash.roll()) {
    CountingBloomFilter<T>::insert(nthash.hashes());
  }
}

template<typename T>
inline uint64_t
KmerCountingBloomFilter<T>::contains(const std::string& seq) const
{
  return contains(seq.c_str(), seq.size());
}

template<typename T>
inline uint64_t
KmerCountingBloomFilter<T>::contains(const char* seq, size_t seq_len) const
{
  uint64_t count = 0;
  NtHash nthash(seq, seq_len, k, CountingBloomFilter<T>::get_hash_num());
  while (nthash.roll()) {
    count += CountingBloomFilter<T>::contains(nthash.hashes());
  }
  return count;
}

} // namespace btllib

#endif
