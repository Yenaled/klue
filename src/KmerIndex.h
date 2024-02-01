#ifndef KLUE_KMERINDEX_H
#define KLUE_KMERINDEX_H

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
  KmerIndex(const ProgramOptions& opt) : k(opt.k), num_trans(0) {
  }

  ~KmerIndex() {}

  void BuildDistinguishingGraph(const ProgramOptions& opt, const std::vector<std::string>& transfasta, bool reconstruct=false, const std::string& tmp_dir="");
  void BuildReconstructionGraph(const ProgramOptions& opt, const std::string& tmp_dir="");
  int k; // k-mer size used
  int num_trans; // number of targets
};

#endif // KLUE_KMERINDEX_H
