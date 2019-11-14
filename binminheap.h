// Copyright 2019 Matt Gara
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/*
 * Binary Minimum Heap
 */

#ifndef BINMINHEAP_H_
#define BINMINHEAP_H_

#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <boost/type_traits.hpp>
#include <boost/static_assert.hpp>

namespace binminheap {

/** A wrapper class to wrap a map-type interface (ex. std::map, or
boost::unordered_map) into a compatible heapscratch dataype.

    XXX: This converter when used with a map_type in the BinaryMinHeap will not
    delete elements from it. In fact, this is a limitation of the implementation
    of the BinaryMinHeap. This will not necessarily mean that efficiency will
    be impacted, in general, it should be as efficient or even more so than
    deletion, however, memory consumption is not optimal.

*/
template <typename map_type> class HeapscratchConverter {
private:
  map_type &map;
  size_t infinity;

public:
  /**
     Sets the underlying map type.

     @param map The underlying map type used internally.
   */
  HeapscratchConverter(map_type &map)
      : map(map), infinity(std::numeric_limits<std::size_t>::max()) {
    BOOST_STATIC_ASSERT(
        (boost::is_same<typename map_type::mapped_type, size_t>::value));
  }

  /**
     Wraps around the underlying [] const implementation of the map. Returns
     infinity if not found which is the necessary behaviour for correct
     performance when used with BinMinHeap.

     @param map The key to choose from the map.

     @return A constant reference to the value in the map.

   */
  template <typename K> const size_t &operator[](const K &key) const {
    typename map_type::const_iterator fit = map.find(key);
    if (fit == map.end()) {
      return infinity;
    }
    return fit->second;
  }

  /**
     Exposes the underlying [] map operator.

     @param map The key to choose from the map.

     @return A reference to the value in the map.

   */
  template <typename K> size_t &operator[](const K &key) { return map[key]; }
};

/** A binary minimum heap templated implementation.

  Improves upon std::priority_queue by supplying updating of priorities in
  heap dynamically. For efficiency, this implementation requires an external
  buffer (refered to as a heapscratch) to maintain a mapping of keys to position
  in heap. If an upperbound on the keys is known, an array (std::vector) is
  suggested as testing shows it to be the most efficient. If no upperbound is
  known a-priori an unordered map (hash table) is suggested, mapping key type to
  std::size_t.

    XXX: Note that due to interface differences of the access operator [] in a
    map like std:map or boost::unordered_map and an array-like data structure
    std::vector, some ground-work must be done to consolidate the use of the
    map structure for a heapscratch. See the binary heap test for an example of
    usage of both. Note that std::numeric_limits<size_t>::max() is a special
    value (infinity) indicating the lack of presence of this element from the
    heap.

    XXX: Note that the key type must be convertable to the std::size_t datatype
    in order for the default heapscratch type (std::vector) to make sense. This
    can be seemlessly accomplished by implementing the cast operator to size_t
    for the key data type supplied.


  */
template <typename P, typename K,
          typename heapscratch_type = std::vector<size_t> >
class BinaryMinHeap {
private:
  std::vector<std::pair<P, K> > data;
  heapscratch_type &key2idx; /* Assume that K is convertable to size_t */

  inline size_t leftChildIdx(size_t k) const {
    // Compute (2 * k + 1) faster
    return (k << 1) + 1;
  }

  inline size_t rightChildIdx(size_t k) const {
    // Compute (2 * k + 2) faster
    return (k + 1) << 1;
  }

  inline size_t parentIdx(size_t k) const {
    // Compute ((k-1) / 2) faster
    return (k - 1) >> 1;
  }

  inline bool lessThan(size_t idx1, size_t idx2) const {
    register size_t dataSize = data.size();
    if (idx1 >= dataSize && idx2 >= dataSize) {
      return false;
    }
    if (idx1 < dataSize && idx2 >= dataSize) {
      return true;
    }
    if (idx1 >= dataSize && idx2 < dataSize) {
      return false;
    }
    return data[idx1].first < data[idx2].first;
  }

  void swapInHeap(size_t idx1, size_t idx2) {
    const K &s1 = data[idx1].second;
    const K &s2 = data[idx2].second;
    assert((size_t)s1 != std::numeric_limits<size_t>::max() &&
           (size_t)s2 != std::numeric_limits<size_t>::max() &&
           data[idx1].second == s1 && data[idx2].second == s2);
    key2idx[s1] = idx2;
    key2idx[s2] = idx1;
    std::swap(data[idx1], data[idx2]);
  }

  void bubbleUp(size_t idx) {
    if (!idx) {
      return;
    }
    size_t parent = parentIdx(idx);
    if (lessThan(idx, parent)) {
      swapInHeap(idx, parent);
      bubbleUp(parent);
    }
  }

  void bubbleDown(size_t idx) {

    if (idx >= data.size()) {
      return;
    }

    size_t left = leftChildIdx(idx);
    size_t right = rightChildIdx(idx);
    size_t minIdx = lessThan(left, right) ? left : right;

    if (lessThan(minIdx, idx)) {
      swapInHeap(idx, minIdx);
      bubbleDown(minIdx);
    }
  }

public:
  /**
     Does nothing to populate heap. Simply sets the underlying heapscratch data
     structure.

     @param heapscratch The heapscratch to use internally.
   */
  BinaryMinHeap(heapscratch_type &heapscratch, const size_t sizehint = 0)
      : key2idx(heapscratch) {
    if (sizehint) {
      data.reserve(sizehint);
    }
  }

  /**
     Construct a heap from a list described by two iterators. Keys in range
     are assumed to be unique, if not behaviour is undefined.

     This construction allows for filling a heap in O(n) time
     as is described in The Algorithm Design Manual (Steven Skiena).

     @param heapscratch The heapscratch to use internally.

     @param begin The starting iterator of the list to populate with.

     @param end The ending iterator of the list to populate with.
  */
  template <typename Iterator>
  BinaryMinHeap(heapscratch_type &heapscratch, Iterator begin, Iterator end)
      : key2idx(heapscratch) {
#ifndef NDEBUG
    std::vector<K> uniqkeys;
    for (Iterator it = begin; it != end; ++it) {
      uniqkeys.push_back(it->second);
    }
    std::sort(uniqkeys.begin(), uniqkeys.end());
    typename std::vector<K>::const_iterator ait =
        adjacent_find(uniqkeys.begin(), uniqkeys.end());
    assert(ait == uniqkeys.end() &&
           "Can not have duplicate keys in range when populating heap.");
#endif
    data.resize(std::distance(begin, end));
    std::copy(begin, end, data.begin());
    size_t k = 0;
    for (Iterator it = begin; it != end; ++it, ++k) {
      key2idx[it->second] = k;
    }
    for (int k = data.size() - 1; k >= 0; k--) {
      bubbleDown(k);
    }
  }

  void pare() {
    std::vector<std::pair<P, K> > tmp;
    data.swap(tmp);
  }

  void clear() { data.clear(); }

  /**
     Get pointer to underlying array implementation.

     This method is only provided as a convenience function to iterate over
     the underlying data in the heap efficiently. It is strongly
     suggested that the underlying data is *not* modified unless
     it is well understood what is to be done will not break the heap
     condition.

     @return A pointer to the start of the internal data array.
  */
  const std::pair<P, K> *values() const { return &data[0]; }

  /**
      Get the size of the heap.

      @return The size of the heap.
  */
  size_t size() { return data.size(); }

  /**
      Find key in the heap.

      Finds the requested key in the heap if it exists. If so,
      the (current) priority of the key in the heap is also
      returned.

      @param key The key to find in the heap.

      @param priority The priority of the found key, only if key is found,
     otherwise this is untouched.

      @return True if key is found, false otherwise.

   */
  bool findKey(const K &key, P &priority) const {
    size_t idx = key2idx[key];
    if (idx == std::numeric_limits<size_t>::max()) {
      return false;
    }
    assert(data[idx].second == key);
    priority = data[idx].first;
    return true;
  }

  /**
      Push priority, key pair onto heap. Keys pushed must be unique, this is not
     checked for
      if debugging is disabled.

      @param keyval The pair of priority and key values to be pushed onto heap.

   */
  void push(const std::pair<P, K> &keyval) {
#ifndef NDEBUG
    const heapscratch_type &hs = key2idx;
    size_t idx = hs[keyval.second];
    assert(idx == std::numeric_limits<size_t>::max() &&
           "Non unique key attempted to be inserted.");
#endif
    data.push_back(keyval);
    size_t idxInserted = data.size() - 1;
    key2idx[keyval.second] = idxInserted;
    bubbleUp(idxInserted);
  }

  /**
      Get top priority, key pair in heap.

      @return The top (minimum) priority, key pair in heap.

   */
  const std::pair<P, K> &top() const { return data[0]; }

  /**
      Delete the top priority, key pair in heap.
   */
  void pop() {
    assert(data.size() > 0 && "Can not pop an empty heap.");
    key2idx[data[0].second] = std::numeric_limits<size_t>::max();
    data[0] = data[data.size() - 1];
    if (data.size() > 1) {
      key2idx[data[data.size() - 1].second] = 0;
    }
    data.erase(data.end() - 1);
    bubbleDown(0);
  }

  /**
      Update the priority of a given key in the heap.

      Updating keys through this interface guarantees that the heap
      condition will not be violated.

      @param key The key to update.

      @param newpriority The new priority to use when updating.
   */
  void updateKey(const K &key, const P &newpriority) {
    size_t idx = key2idx[key];
    if (idx == std::numeric_limits<size_t>::max()) {
      return;
    }
    assert(data[idx].second == key);
    P oldpriority = data[idx].first;
    data[idx].first = newpriority;
    if (oldpriority < newpriority) {
      bubbleDown(idx);
    } else if (newpriority < oldpriority) {
      bubbleUp(idx);
    } else {
      // oldpriority == newpriority NOTHING TO DO
    }
  }

  template <typename Pp, typename Kk, typename Ss>
  friend std::ostream &operator<<(std::ostream &os,
                                  const BinaryMinHeap<Pp, Kk, Ss> &mh);
};

template <typename Pp, typename Kk, typename Ss>
std::ostream &operator<<(std::ostream &os,
                         const BinaryMinHeap<Pp, Kk, Ss> &mh) {

  boost::unordered_map<std::pair<size_t, size_t>, size_t> rc2idx;
  boost::unordered_map<size_t, std::pair<size_t, size_t> > idx2rc;
  size_t numlevels = std::log(mh.data.size()) / std::log(2) + 1;
  for (int level = int(numlevels) - 1; level >= 0; level--) {
    size_t start = (1 << level) - 1;
    size_t num2print = 1 << level;

    if (level == (int)numlevels - 1) {

      for (size_t k = 0; k < num2print; k++) {
        size_t idx = start + k;
        std::pair<size_t, size_t> rc = std::make_pair(2 * k, 0);
        idx2rc[idx] = rc;
        rc2idx[rc] = idx;
      }

    } else {

      for (size_t k = 0; k < num2print; k++) {
        size_t idx = start + k;
        if (idx >= mh.data.size()) {
          break;
        }
        size_t lcidx = mh.leftChildIdx(idx);
        size_t rcidx = mh.rightChildIdx(idx);
        std::pair<size_t, size_t> lftrc = idx2rc[lcidx];
        std::pair<size_t, size_t> rgtrc = idx2rc[rcidx];
        std::pair<size_t, size_t> rc =
            std::make_pair((lftrc.first + rgtrc.first) / 2, lftrc.second + 1);
        idx2rc[idx] = rc;
        rc2idx[rc] = idx;
      }
    }
  }

  boost::unordered_map<std::pair<size_t, size_t>, size_t>::iterator it;

  size_t maxCol, maxRow;
  maxCol = 0;
  maxRow = 0;
  for (it = rc2idx.begin(); it != rc2idx.end(); it++) {
    size_t r = it->first.first;
    size_t c = it->first.second;
    if (r > maxRow) {
      maxRow = r;
    }
    if (c > maxCol) {
      maxCol = c;
    }
  }

  const size_t padsize = 3;
  for (size_t r = 0; r <= maxRow; r++) {
    for (size_t c = 0; c <= maxCol; c++) {
      boost::unordered_map<std::pair<size_t, size_t>, size_t>::iterator fit;
      fit = rc2idx.find(std::make_pair(r, c));
      bool found = fit != rc2idx.end();
      size_t idx = 0;
      if (found) {
        idx = fit->second;
      }
      bool withinRange = found && idx < mh.data.size();
      if (!withinRange) {
        os << std::setfill(' ') << " " << std::setw(padsize) << ' ' << "  "
           << std::setw(padsize) << ' ' << "  ";
      } else {
        os << std::setfill(' ') << "(" << std::setw(padsize)
           << mh.data[idx].first << ", " << std::setw(padsize)
           << mh.data[idx].second << ") ";
      }
    }
    os << std::endl;
  }
  return os;
}
}
#endif
