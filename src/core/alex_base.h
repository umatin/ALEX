// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

/* This file contains the classes for linear models and model builders, helpers
 * for the bitmap,
 * cost model weights, statistic accumulators for collecting cost model
 * statistics,
 * and other miscellaneous functions
 */

#pragma once

#include "log_regression.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include <bitset>
#include <cassert>
#ifdef _WIN32
#include <intrin.h>
#include <limits.h>
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

#ifdef _MSC_VER
#define forceinline __forceinline
#elif defined(__GNUC__)
#define forceinline inline __attribute__((__always_inline__))
#elif defined(__CLANG__)
#if __has_attribute(__always_inline__)
#define forceinline inline __attribute__((__always_inline__))
#else
#define forceinline inline
#endif
#else
#define forceinline inline
#endif

namespace alex {

/*** Linear model and model builder ***/

// Forward declaration
template <class T>
class LinearModelBuilder;

// Linear regression model
template <class T>
class LinearModel {
 public:
  double a_ = 0;  // slope
  double b_ = 0;  // intercept

  LinearModel() = default;
  LinearModel(double a, double b) : a_(a), b_(b) {}
  explicit LinearModel(const LinearModel& other) : a_(other.a_), b_(other.b_) {}

  void expand(double expansion_factor) {
    a_ *= expansion_factor;
    b_ *= expansion_factor;
  }

  inline int predict(T key) const {
    return static_cast<int>(a_ * static_cast<double>(key) + b_);
  }

  inline double predict_double(T key) const {
    return a_ * static_cast<double>(key) + b_;
  }
};

template <class T>
class LinearModelBuilder {
 public:
  LinearModel<T>* model_;
  
  std::vector<double> x_data;
  std::vector<double> y_data;

  explicit LinearModelBuilder<T>(LinearModel<T>* model) : model_(model) {}

  inline void add(T x, int y) {
    x_data.push_back(static_cast<double>(x));
    y_data.push_back(static_cast<double>(y));
    min_y = std::min(min_y, y);
    max_y = std::max(max_y, y);
  }

  void build() {
    if (y_data.size() == 0) {
      model_->a_ = 0;
      model_->b_ = 0;
      return;
    }
    if (y_data.size() == 1) {
      model_->a_ = 0;
      model_->b_ = static_cast<double>(y_data[0]);
      return;
    }
    int fr,to;
    create_regression_tournament_selection<FastDiscreteLogNorm,true,true>(x_data, y_data, fr, to, std::min<int>(20, (int)log2(x_data.size())), min_y, max_y);
    fit_line(x_data, y_data, fr, to, model_->a_, model_->b_);
    //std::cout << model_->a_ << model_->b_ << std::endl;
  }

 private:
  int min_y = INT32_MAX;
  int max_y = 0;
  T x_min_ = std::numeric_limits<T>::max();
  T x_max_ = std::numeric_limits<T>::lowest();
  double y_min_ = std::numeric_limits<double>::max();
  double y_max_ = std::numeric_limits<double>::lowest();
};

/*** Comparison ***/

struct AlexCompare {
  template <class T1, class T2>
  bool operator()(const T1& x, const T2& y) const {
    static_assert(
        std::is_arithmetic<T1>::value && std::is_arithmetic<T2>::value,
        "Comparison types must be numeric.");
    return x < y;
  }
};

/*** Helper methods for bitmap ***/

// Extract the rightmost 1 in the binary representation.
// e.g. extract_rightmost_one(010100100) = 000000100
inline uint64_t extract_rightmost_one(uint64_t value) {
  return value & -static_cast<int64_t>(value);
}

// Remove the rightmost 1 in the binary representation.
// e.g. remove_rightmost_one(010100100) = 010100000
inline uint64_t remove_rightmost_one(uint64_t value) {
  return value & (value - 1);
}

// Count the number of 1s in the binary representation.
// e.g. count_ones(010100100) = 3
inline int count_ones(uint64_t value) {
  return static_cast<int>(_mm_popcnt_u64(value));
}

// Get the offset of a bit in a bitmap.
// word_id is the word id of the bit in a bitmap
// bit is the word that contains the bit
inline int get_offset(int word_id, uint64_t bit) {
  return (word_id << 6) + count_ones(bit - 1);
}

/*** Cost model weights ***/

// Intra-node cost weights
double kExpSearchIterationsWeight = 20;
double kShiftsWeight = 0.5;

// TraverseToLeaf cost weights
double kNodeLookupsWeight = 20;
double kModelSizeWeight = 5e-7;

/*** Stat Accumulators ***/

struct DataNodeStats {
  double num_search_iterations = 0;
  double num_shifts = 0;
};

// Used when stats are computed using a sample
struct SampleDataNodeStats {
  double log2_sample_size = 0;
  double num_search_iterations = 0;
  double log2_num_shifts = 0;
};

// Accumulates stats that are used in the cost model, based on the actual vs
// predicted position of a key
class StatAccumulator {
 public:
  virtual ~StatAccumulator() = default;
  virtual void accumulate(int actual_position, int predicted_position) = 0;
  virtual double get_stat() = 0;
  virtual void reset() = 0;
};

// Mean log error represents the expected number of exponential search
// iterations when doing a lookup
class ExpectedSearchIterationsAccumulator : public StatAccumulator {
 public:
  void accumulate(int actual_position, int predicted_position) override {
    cumulative_log_error_ +=
        std::log2(std::abs(predicted_position - actual_position) + 1);
    count_++;
  }

  double get_stat() override {
    if (count_ == 0) return 0;
    return cumulative_log_error_ / count_;
  }

  void reset() override {
    cumulative_log_error_ = 0;
    count_ = 0;
  }

 public:
  double cumulative_log_error_ = 0;
  int count_ = 0;
};

// Mean shifts represents the expected number of shifts when doing an insert
class ExpectedShiftsAccumulator : public StatAccumulator {
 public:
  explicit ExpectedShiftsAccumulator(int data_capacity)
      : data_capacity_(data_capacity) {}

  // A dense region of n keys will contribute a total number of expected shifts
  // of approximately
  // ((n-1)/2)((n-1)/2 + 1) = n^2/4 - 1/4
  // This is exact for odd n and off by 0.25 for even n.
  // Therefore, we track n^2/4.
  void accumulate(int actual_position, int) override {
    if (actual_position > last_position_ + 1) {
      long long dense_region_length = last_position_ - dense_region_start_idx_ + 1;
      num_expected_shifts_ += (dense_region_length * dense_region_length) / 4;
      dense_region_start_idx_ = actual_position;
    }
    last_position_ = actual_position;
    count_++;
  }

  double get_stat() override {
    if (count_ == 0) return 0;
    // first need to accumulate statistics for current packed region
    long long dense_region_length = last_position_ - dense_region_start_idx_ + 1;
    long long cur_num_expected_shifts =
        num_expected_shifts_ + (dense_region_length * dense_region_length) / 4;
    return cur_num_expected_shifts / static_cast<double>(count_);
  }

  void reset() override {
    last_position_ = -1;
    dense_region_start_idx_ = 0;
    num_expected_shifts_ = 0;
    count_ = 0;
  }

 public:
  int last_position_ = -1;
  int dense_region_start_idx_ = 0;
  long long num_expected_shifts_ = 0;
  int count_ = 0;
  int data_capacity_ = -1;  // capacity of node
};

// Combines ExpectedSearchIterationsAccumulator and ExpectedShiftsAccumulator
class ExpectedIterationsAndShiftsAccumulator : public StatAccumulator {
 public:
  ExpectedIterationsAndShiftsAccumulator() = default;
  explicit ExpectedIterationsAndShiftsAccumulator(int data_capacity)
      : data_capacity_(data_capacity) {}

  void accumulate(int actual_position, int predicted_position) override {
    cumulative_log_error_ +=
        std::log2(std::abs(predicted_position - actual_position) + 1);

    if (actual_position > last_position_ + 1) {
      long long dense_region_length = last_position_ - dense_region_start_idx_ + 1;
      num_expected_shifts_ += (dense_region_length * dense_region_length) / 4;
      dense_region_start_idx_ = actual_position;
    }
    last_position_ = actual_position;

    count_++;
  }

  double get_stat() override {
    assert(false);  // this should not be used
    return 0;
  }

  double get_expected_num_search_iterations() {
    if (count_ == 0) return 0;
    return cumulative_log_error_ / count_;
  }

  double get_expected_num_shifts() {
    if (count_ == 0) return 0;
    long long dense_region_length = last_position_ - dense_region_start_idx_ + 1;
    long long cur_num_expected_shifts =
        num_expected_shifts_ + (dense_region_length * dense_region_length) / 4;
    return cur_num_expected_shifts / static_cast<double>(count_);
  }

  void reset() override {
    cumulative_log_error_ = 0;
    last_position_ = -1;
    dense_region_start_idx_ = 0;
    num_expected_shifts_ = 0;
    count_ = 0;
  }

 public:
  double cumulative_log_error_ = 0;
  int last_position_ = -1;
  int dense_region_start_idx_ = 0;
  long long num_expected_shifts_ = 0;
  int count_ = 0;
  int data_capacity_ = -1;  // capacity of node
};

/*** Miscellaneous helpers ***/

// https://stackoverflow.com/questions/364985/algorithm-for-finding-the-smallest-power-of-two-thats-greater-or-equal-to-a-giv
inline int pow_2_round_up(int x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return x + 1;
}

// https://stackoverflow.com/questions/994593/how-to-do-an-integer-log2-in-c
inline int log_2_round_down(int x) {
  int res = 0;
  while (x >>= 1) ++res;
  return res;
}

// https://stackoverflow.com/questions/1666093/cpuid-implementations-in-c
class CPUID {
  uint32_t regs[4];

 public:
  explicit CPUID(unsigned i, unsigned j) {
#ifdef _WIN32
    __cpuidex((int*)regs, (int)i, (int)j);
#else
    asm volatile("cpuid"
                 : "=a"(regs[0]), "=b"(regs[1]), "=c"(regs[2]), "=d"(regs[3])
                 : "a"(i), "c"(j));
#endif
  }

  const uint32_t& EAX() const { return regs[0]; }
  const uint32_t& EBX() const { return regs[1]; }
  const uint32_t& ECX() const { return regs[2]; }
  const uint32_t& EDX() const { return regs[3]; }
};

// https://en.wikipedia.org/wiki/CPUID#EAX=7,_ECX=0:_Extended_Features
bool cpu_supports_bmi() {
  return static_cast<bool>(CPUID(7, 0).EBX() & (1 << 3));
}
}