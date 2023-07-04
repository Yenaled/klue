#ifndef KURE_COMMON_H
#define KURE_COMMON_H

#define KURE_VERSION "0.29.0"

// NOTE: MAKE SURE THIS FILE GETS INCLUDED FIRST IN ALL OTHER FILES AND BEFORE ANY EXTERNAL LIBRARIES

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include "kseq.h"

//#if defined(_MSVC_LANG)
//#define KALLISTO_CPP_VERSION _MSVC_LANG
//#else
//#define KURE_CPP_VERSION __cplusplus
//#endif
//#if KURE_CPP_VERSION < 201703L
//#include "robin_hood.h"
//#define u_map_ robin_hood::unordered_flat_map
//#define u_set_ robin_hood::unordered_set
//#else
//#include "unordered_dense.h"
//#define u_map_ ankerl::unordered_dense::map
//#define u_set_ ankerl::unordered_dense::set
//#endif

#define u_map_ std::unordered_map
#define u_set_ std::unordered_set


#ifdef _WIN64
typedef unsigned int uint;
#endif

struct ProgramOptions {
  bool verbose;
  int threads;
  std::string distinguish_output_fasta;
  int k;
  int g;
  bool stream_out;
  bool distinguish_all_but_one_color;
  bool distinguish_union;
  int distinguish_range_begin;
  int distinguish_range_end;
  std::vector<int> kmer_multiplicity;
  std::vector<int> inner;
  std::vector<std::string> transfasta;

ProgramOptions() :
  verbose(false),
  threads(1),
  k(31),
  g(0),
  stream_out(false),
  distinguish_all_but_one_color(false),
  distinguish_union(false),
  distinguish_range_begin(0),
  distinguish_range_end(0)
  {}
};

std::string pretty_num(size_t num);
std::string pretty_num(int64_t num);
std::string pretty_num(int num);


#ifdef KURE_USE_ZLIB_NG
#include "zlib-ng/zlib.h"
#else
#include <zlib.h>
#endif

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif


#endif // KURE_COMMON_H
