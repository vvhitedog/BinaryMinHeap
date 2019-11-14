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

// A simple example showing how to update keys in a BinaryMinHeap
#include <binminheap.h>
#include <map>

using namespace binminheap;
using namespace std;

int main(int argc, char *argv[]) {

  cout << "example updating keys in a BinaryMinHeap: " << endl;
  std::vector<size_t> heapscratch(4, std::numeric_limits<size_t>::max());
  BinaryMinHeap<double, int> heap(heapscratch);
  heap.push(make_pair(.2, 2));
  heap.push(make_pair(.1, 3));
  heap.push(make_pair(.3, 1));
  cout << "heap: " << endl << heap << endl;
  // Now we update the key 3 to be the largest (priority = 1.)
  heap.updateKey(3,1.);
  cout << "heap: " << endl << heap << endl;
}
