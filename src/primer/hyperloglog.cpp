//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// hyperloglog.cpp
//
// Identification: src/primer/hyperloglog.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/hyperloglog.h"

namespace bustub {

/** @brief Parameterized constructor. */
template <typename KeyType>
HyperLogLog<KeyType>::HyperLogLog(int16_t n_bits) : cardinality_(0) {
  if (n_bits <= 0) {
    registers_.resize(1);
  } else {
    registers_.resize(1 << n_bits);
  }
}

/**
 * @brief Function that computes binary.
 *
 * @param[in] hash
 * @returns binary of a given hash
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeBinary(const hash_t &hash) const -> std::bitset<BITSET_CAPACITY> {
  return std::bitset<BITSET_CAPACITY>(hash);
}

/**
 * @brief Function that computes leading zeros.
 *
 * @param[in] bset - binary values of a given bitset
 * @returns leading zeros of given binary set
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::PositionOfLeftmostOne(const std::bitset<BITSET_CAPACITY> &bset) const -> uint64_t {
  for (uint64_t i = BITSET_CAPACITY; i > 0; i--) {
    if (bset.test(i - 1)) {
      return i;
    }
  }
  return 0;
}

/**
 * @brief Adds a value into the HyperLogLog.
 *
 * @param[in] val - value that's added into hyperloglog
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::AddElem(KeyType val) -> void {
  hash_t hash = CalculateHash(val);
  std::bitset<BITSET_CAPACITY> binary = ComputeBinary(hash);

  // Get the number of buckets (registers)
  size_t num_buckets = registers_.size();

  // Calculate the number of bits needed for bucket index
  size_t bits_for_index = 0;
  size_t temp = num_buckets;
  while (temp > 1) {
    temp >>= 1;
    bits_for_index++;
  }

  // Extract bucket index from the MOST significant bits
  size_t bucket_index = 0;
  for (size_t i = 0; i < bits_for_index; i++) {
    if (binary.test(BITSET_CAPACITY - 1 - i)) {
      bucket_index |= (1ULL << (bits_for_index - 1 - i));
    }
  }

  // Create a bitset for the remaining bits (after removing the index bits from MSB)
  std::bitset<BITSET_CAPACITY> remaining_bits;
  for (size_t i = bits_for_index; i < BITSET_CAPACITY; i++) {
    remaining_bits.set(BITSET_CAPACITY - 1 - i, binary.test(BITSET_CAPACITY - 1 - i));
  }

  // Get position of leftmost 1 in remaining bits
  uint64_t position = PositionOfLeftmostOne(remaining_bits);
  // Calculate rho: if position is 0, means all zeros, rho = BITSET_CAPACITY - bits_for_index + 1
  // Otherwise, rho = (BITSET_CAPACITY - bits_for_index) - position + 1
  uint64_t rho = (BITSET_CAPACITY - bits_for_index) - position + 1;

  // Update the register with the maximum value
  std::bitset<BITSET_CAPACITY> new_value(rho);
  if (new_value.to_ullong() > registers_[bucket_index].to_ullong()) {
    registers_[bucket_index] = new_value;
  }
}

/**
 * @brief Function that computes cardinality.
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeCardinality() -> void {
  size_t num_buckets = registers_.size();
  double sum = 0.0;

  for (size_t i = 0; i < num_buckets; i++) {
    uint64_t register_value = registers_[i].to_ullong();
    if (register_value >= 64) {
      // Avoid undefined behavior for shift >= 64
      sum += 0.0;
    } else {
      sum += 1.0 / (1ULL << register_value);
    }
  }

  double raw_estimate = CONSTANT * num_buckets * num_buckets / sum;
  cardinality_ = static_cast<size_t>(raw_estimate);
}

template class HyperLogLog<int64_t>;
template class HyperLogLog<std::string>;

}  // namespace bustub
