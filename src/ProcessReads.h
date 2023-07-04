#ifndef KURE_PROCESSREADS_H
#define KURE_PROCESSREADS_H

#include "common.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

class MasterProcessor;

int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt);

class SequenceReader {
public:
  
  SequenceReader(const ProgramOptions& opt) :
  readbatch_id(-1) {};
  SequenceReader() : state(false), readbatch_id(-1) {};
  virtual ~SequenceReader() {}
  
  virtual bool empty() = 0;
  virtual void reset();
  virtual void reserveNfiles(int n) = 0;
  virtual bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                              std::vector<std::pair<const char*, int>>& names,
                              std::vector<std::pair<const char*, int>>& quals,
                              std::vector<uint32_t>& flags,
                              int &readbatch_id,
                              bool full=false,
                              bool comments=false) = 0;
  
  
public:
  bool state; // is the file open
  int readbatch_id = -1;
};

class FastqSequenceReader : public SequenceReader {
public:
  
  FastqSequenceReader(const ProgramOptions& opt) : SequenceReader(opt),
  current_file(0), files(opt.transfasta) {
    SequenceReader::state = false;
    nfiles = opt.transfasta.size();
    reserveNfiles(nfiles);
  }
  FastqSequenceReader() : SequenceReader(), 
  current_file(0) {};
  FastqSequenceReader(FastqSequenceReader &&o);
  ~FastqSequenceReader();
  
  bool empty();  
  void reset();
  void reserveNfiles(int n);
  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      int &readbatch_id,
                      bool full=false,
                      bool comments=false);
  
public:
  int nfiles = 1;
  uint32_t numreads = 0;
  std::vector<gzFile> fp;
  std::vector<int> l;
  std::vector<int> nl;
  std::vector<std::string> files;
  int current_file;
  std::vector<kseq_t*> seq;
  int interleave_nfiles;
};

class MasterProcessor {
public:
  MasterProcessor (const ProgramOptions& opt)
    : opt(opt), numreads(0), bufsize(1ULL<<23), curr_readbatch_id(0) { 
    
    SR = new FastqSequenceReader(opt);
    verbose = opt.verbose;
    nfiles = opt.transfasta.size();
    auto f = opt.distinguish_output_fasta;
    out = fopen(f.c_str(), "wb");
  }
  
  ~MasterProcessor() {
    fclose(out);
    delete SR;
  }
  
  std::mutex reader_lock;
  std::vector<std::mutex> parallel_reader_locks;
  bool parallel_read;
  std::mutex writer_lock;
  std::condition_variable cv;
  
  FILE* out;
  bool verbose;
  
  SequenceReader *SR;
  std::vector<FastqSequenceReader> FSRs;
  
  const ProgramOptions& opt;
  int64_t numreads;
  size_t bufsize;
  int nfiles;
  int curr_readbatch_id;
  
  void processReads();
  void update(int n,
              std::vector<std::pair<const char*, int>>& seqs,
              std::vector<std::pair<const char*, int>>& names,
              std::vector<std::pair<const char*, int>>& quals,
              std::vector<uint32_t>& flags,
              int readbatch_id);  
};

class ReadProcessor {
public:
  ReadProcessor(const ProgramOptions& opt, MasterProcessor& mp);
  ReadProcessor(ReadProcessor && o);
  ~ReadProcessor();
  char *buffer;
  
  size_t bufsize;
  MasterProcessor& mp;
  int64_t numreads;
  
  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;
  
  bool full;
  bool comments;
  
  /*std::vector<std::vector<int>> newIDs;
   std::vector<std::vector<int>> IDs;*/
  
  void operator()();
  void processBuffer();
  void clear();
};

std::string pretty_num(size_t num);

#endif // KURE_PROCESSREADS_H
