// SPDX-FileCopyrightText: 2022 Moritz Bültmann <moritz.bueltmann@gmx.de>
// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef SRC_STL_ALGORITHMS_HPP_
#define SRC_STL_ALGORITHMS_HPP_
/** \file stl_algorithms.hpp
 *  \brief This file contains algorithms for C++ STL objects.
 */
// Includes
#include <algorithm>
#include <utility>
#include <vector>
namespace stl_algorithm {
// _____________________________________________________________________________
/** \brief This function inserts a whole list of elements to a std::vector
 *
 *  This algorithm goes through the vector once (linear complexity).
 *  The sorting algorithm is N*log(N), but is berformed on the much smaller
 *  insertion index array.
 */
template <typename T>
void insert(
    std::vector<std::pair<size_t, T>>& indices, std::vector<T>& v) {
  // Resize the array to new size
  size_t old_size{v.size()};
  v.resize(v.size() + indices.size());
  std::sort(indices.begin(), indices.end());
  // Next spot to fill; starts at the kast position of the resized vector
  long int next{static_cast<int>(v.size()) - 1};
  long int section_begin;
  long int section_end;
  for (auto it = indices.rbegin(); it != indices.rend(); ++it) {  // reverse
    section_end = it->first;
    if (it != indices.rbegin()) { section_begin = (it - 1)->first - 1; }
    else { section_begin = old_size - 1; }
    for (long int i = section_begin; i >= section_end; --i) {
      v[next] = v[i];
      --next;
    }
    v[next] = it->second;
    --next;
  }
}
// _____________________________________________________________________________
/** \brief This function removes a whole list of elements from a std::vector
 *
 *  This algorithm goes through the vector once (linear complexity).
 *  The sorting algorithm is N*log(N), but is berformed on the much smaller
 *  deletion index array.
 */
template <typename T>
void erase(std::vector<size_t>& indices, std::vector<T>& v) {
  // Next spot to fill; starts at first element to be erased
  size_t next{indices[0]};
  size_t section_begin;
  size_t section_end;
  std::sort(indices.begin(), indices.end());
  for(auto it = indices.begin(); it != indices.end(); ++it) {
    // This loop iterates over the data between two erased elements.
    // The edge cases are covered by the if-condition.
    section_begin = *it + 1;
    if (it != indices.end() - 1) { section_end = *(it + 1); }
    else { section_end = v.size(); }
    for (size_t i = section_begin; i < section_end; ++i) {
      v[next] = v[i];
      ++next;
    }
  }
  v.resize(v.size() - indices.size());
}
// _____________________________________________________________________________
}
#endif  // SRC_STL_ALGORITHMS_HPP_
