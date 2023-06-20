#ifndef KURE_KMERINDEX_H
#define KURE_KMERINDEX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <stdint.h>
#include <iostream>
#include <numeric>
#include <limits>

#include "common.h"
#include "Kmer.hpp"
#include "CompactedDBG.hpp"


struct KmerIndex {
  KmerIndex(const ProgramOptions& opt) : k(opt.k), {
  }

  ~KmerIndex() {}

  void BuildDistinguishingGraph(const ProgramOptions& opt, std::ofstream& out);
  int k; // k-mer size used
  int num_trans; // number of targets
};

#endif // KURE_KMERINDEX_H
