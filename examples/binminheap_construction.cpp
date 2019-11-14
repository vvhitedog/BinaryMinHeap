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
//

// A simple example showing three different ways to initialize a heap, using
// different underlying datastructures
#include <binminheap.h>
#include <map>
#include <boost/unordered_map.hpp>

using namespace binminheap;
using namespace std;

int main(int argc, char *argv[]) {

  // Ex 1
  {
    // Since we use a std::vector for the heapscratch,
    // we must make sure that the highest key can be
    // mapped in the vector, hence the need for 4 elements
    // instead of 3.
    cout << "example #1 using heapscratch type std::vector: " << endl;
    std::vector<size_t> heapscratch(4, std::numeric_limits<size_t>::max());
    BinaryMinHeap<double, int> heap(heapscratch);
    heap.push(make_pair(.1, 2));
    heap.push(make_pair(.1, 3));
    heap.push(make_pair(.1, 1));
    cout << "heap: " << endl << heap << endl;
    heap.pop();
    cout << "heap: " << endl << heap << endl;
  }

  // Ex 2
  {
    cout << "example #2 using heapscratch type boost::unordered_map: " << endl;
    typedef boost::unordered_map<int, size_t> heapscratch_map;
    typedef HeapscratchConverter<heapscratch_map> heapscratch_type;
    heapscratch_map heapscratchMap;
    heapscratch_type heapscratch(heapscratchMap);
    BinaryMinHeap<double, int, heapscratch_type> heap(heapscratch);
    heap.push(make_pair(.1, 2));
    heap.push(make_pair(.1, 3));
    heap.push(make_pair(.1, 1));
    cout << "heap: " << endl << heap << endl;
    heap.pop();
    cout << "heap: " << endl << heap << endl;
  }

  // Ex 3
  {
    cout << "example #3 using heapscratch type boost::unordered_map with O(n) "
            "constructor: " << endl;
    typedef boost::unordered_map<int, size_t> heapscratch_map;
    typedef HeapscratchConverter<heapscratch_map> heapscratch_type;
    heapscratch_map heapscratchMap;
    heapscratch_type heapscratch(heapscratchMap);
    map<double, int> prioritykeys;
    for (size_t i = 0; i < 10; i++) {
      prioritykeys[2 * i] = i;
    }
    BinaryMinHeap<double, int, heapscratch_type> heap(
        heapscratch, prioritykeys.begin(), prioritykeys.end());
    cout << "heap: " << endl << heap << endl;
    heap.pop();
    cout << "heap: " << endl << heap << endl;
  }

  return 0;
}
