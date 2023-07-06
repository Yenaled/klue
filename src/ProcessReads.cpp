#include <fstream>
#include <limits>
#include <iomanip>
#include "common.h"
#include "ProcessReads.h"
#include "kseq.h"
#include <unordered_set>
#include "StringSearch.h"

std::string pretty_num(size_t num) {
  auto s = std::to_string(num);
  auto ret = std::string("");
  
  if (s.size() <= 3) {
    return s;
  }
  
  int remainder = s.size() % 3;
  if (remainder == 0) {
    remainder = 3;
  }
  
  size_t start_pos = 0;
  while (start_pos + remainder < s.size() - 1) {
    ret += s.substr(start_pos, remainder) + ",";
    start_pos += remainder;
    remainder = 3;
  }
  
  ret += s.substr(start_pos, 3);
  
  return ret;
}


//methods


int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt) {
  std::ios_base::sync_with_stdio(false);
  size_t numreads = 0;
  
  if (MP.verbose) {
    std::cerr << "* processing the contigs ..."; std::cerr.flush();
  }
  MP.processContigs();
  numreads = MP.numreads;
  if (MP.verbose) {
    std::cerr << std::endl << "done " << std::endl;
  }
  if (MP.verbose) {
    std::cerr << "* processed " << pretty_num(numreads) << " contigs";
    std::cerr << std::endl;
  }
  
  if (MP.verbose) {
    std::cerr << "* processing the contigs ..."; std::cerr.flush();
  }
  MP.processReads();
  numreads = MP.numreads;
  if (MP.verbose) {
    std::cerr << std::endl << "done " << std::endl;
  }
  if (MP.verbose) {
    std::cerr << "* processed " << pretty_num(numreads) << " sequences";
    std::cerr << std::endl;
  }

  return numreads;
}


/** -- read processors -- **/

void MasterProcessor::processContigs() {

  // start worker threads
  
  std::vector<std::thread> workers;
  
  for (int i = 0; i < opt.threads; i++) {
    workers.emplace_back(std::thread(ReadProcessor(opt,*this)));
  }
    
  // let the workers do their thing
  for (int i = 0; i < opt.threads; i++) {
    workers[i].join(); //wait for them to finish
  }
  delete inSR;
  inSR = nullptr;
}


void MasterProcessor::processReads() {
  
  readSeqs = true;
  
  // start worker threads
  
  std::vector<std::thread> workers;
  
  for (int i = 0; i < opt.threads; i++) {
    workers.emplace_back(std::thread(ReadProcessor(opt,*this)));
  }
  
  // let the workers do their thing
  for (int i = 0; i < opt.threads; i++) {
    workers[i].join(); //wait for them to finish
  }
  delete SR;
  SR = nullptr;
}

void MasterProcessor::update(int n, 
                             std::vector<std::pair<const char*, int>>& seqs,
                             std::vector<std::pair<const char*, int>>& names,
                             std::vector<std::pair<const char*, int>>& quals,
                             std::vector<uint32_t>& flags,
                             int readbatch_id) {
  // acquire the writer lock
  std::unique_lock<std::mutex> lock(this->writer_lock);
  

  while (readbatch_id != curr_readbatch_id) {
    cv.wait(lock, [this, readbatch_id]{ return readbatch_id == curr_readbatch_id; });
  }
  writeOutput(seqs, names, quals, flags);
  
  numreads += n;
  curr_readbatch_id++;
  lock.unlock(); // releases the lock
  cv.notify_all(); // Alert all other threads to check their readbatch_id's!
}

void MasterProcessor::writeOutput(std::vector<std::pair<const char*, int>>& seqs,
                                  std::vector<std::pair<const char*, int>>& names,
                                  std::vector<std::pair<const char*, int>>& quals,
                                  std::vector<uint32_t>& flags) {
  // Write out fastq
  int incf, jmax;
  incf = nfiles-1;
  jmax = nfiles;
  
  std::vector<const char*> s(jmax, nullptr);
  std::vector<const char*> n(jmax, nullptr);
  std::vector<const char*> nl(jmax, nullptr);
  std::vector<const char*> q(jmax, nullptr);
  std::vector<int> l(jmax,0);
  char start_char = '>';

  size_t readnum = 0;
  for (int i = 0; i + incf < seqs.size(); i++, readnum++) {
    for (int j = 0; j < jmax; j++) {
      // TODO:
    }
    i += incf;
  }
}

ReadProcessor::ReadProcessor(const ProgramOptions& opt, MasterProcessor& mp) : 
  mp(mp), numreads(0) {
  // initialize buffer
  bufsize = mp.bufsize;
  buffer = new char[bufsize];
  seqs.reserve(bufsize/50);
  clear();
}

ReadProcessor::ReadProcessor(ReadProcessor && o) :
  bufsize(o.bufsize),
  mp(o.mp),
  numreads(o.numreads),
  seqs(std::move(o.seqs)),
  names(std::move(o.names)),
  quals(std::move(o.quals)),
  flags(std::move(o.flags)),
  full(o.full),
  comments(o.comments) {
  buffer = o.buffer;
  o.buffer = nullptr;
  o.bufsize = 0;
}

ReadProcessor::~ReadProcessor() {
  if (buffer != nullptr) {
    delete[] buffer;
    buffer = nullptr;
  }
}

void ReadProcessor::operator()() {
  while (true) {
    int readbatch_id;
    std::vector<std::string> umis;
    FastqSequenceReader *SR;
    bool full = false;
    bool comments = false;
    if (mp.readSeqs) {
      SR = mp.SR;
    } else {
      full = true; // We need to get the name
      SR = mp.inSR;
    }
    // grab the reader lock
    {
      std::lock_guard<std::mutex> lock(mp.reader_lock);
      if (SR->empty()) {
        // nothing to do
        return;
      } else {
        // get new sequences
        SR->fetchSequences(buffer, bufsize, seqs, names, quals, flags, readbatch_id, full, comments);
      }
      // release the reader lock
    }
    
    // process our sequences
    if (mp.readSeqs) {
      processBuffer();
    } else {
      processBufferContigs();
    }
    
    // update the results, MP acquires the lock
    int nfiles = SR->nfiles;
    mp.update(seqs.size() / nfiles, seqs, names, quals, flags, readbatch_id);
    clear();
  }
}

void ReadProcessor::processBuffer() {
  // actually process the sequence
  
  int incf, jmax, nfiles;
  nfiles = mp.nfiles;
  incf = nfiles-1;
  jmax = nfiles;
  
  std::vector<const char*> s(jmax, nullptr);
  std::vector<int> l(jmax,0);

  for (int i = 0; i + incf < seqs.size(); i++) {
    for (int j = 0; j < jmax; j++) {
      s[j] = seqs[i+j].first;
      l[j] = seqs[i+j].second;
    }
    i += incf;
    numreads++;
    if (numreads > 0 && numreads % 1000000 == 0 && mp.verbose) {
      numreads = 0; // reset counter
      std::cerr << '\r' << (mp.numreads/1000000) << "M reads processed";
      std::cerr << "         ";
      std::cerr.flush();
    }
  }
}

void ReadProcessor::processBufferContigs() {
  // actually process the sequence
  
  int incf, jmax, nfiles;
  nfiles = 1; //nfiles = mp.nfiles;
  incf = nfiles-1;
  jmax = nfiles;

  for (int i = 0; i + incf < seqs.size(); i++) {
    for (int j = 0; j < jmax; j++) {
      // Debug:
      std::cout << std::string(names[i+j].first, names[i+j].second) << " " << std::string(s[j], l[j]) << std::endl;
    }
    i += incf;
    numreads++;
    if (numreads > 0 && numreads % 1000000 == 0 && mp.verbose) {
      numreads = 0; // reset counter
      std::cerr << '\r' << (mp.numreads/1000000) << "M contigs processed";
      std::cerr << "         ";
      std::cerr.flush();
    }
  }
  
  exit(0); // TODO: remove
}

void ReadProcessor::clear() {
  memset(buffer,0,bufsize);
}

/** -- sequence readers -- **/

void SequenceReader::reset() {
  state = false;
  readbatch_id = -1;
}

FastqSequenceReader::~FastqSequenceReader() {
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
  }
  
  for (auto &s : seq) {
    kseq_destroy(s);
  }
}


bool FastqSequenceReader::empty() {
  return (!state && current_file >= files.size());
}

void FastqSequenceReader::reset() {
  SequenceReader::reset();
  
  for (auto &f : fp) {
    if (f) {
      gzclose(f);
    }
    f = nullptr;
  }
  
  for (auto &ll : l) {
    ll = 0;
  }
  for (auto &nll : nl) {
    nll = 0;
  }
  
  current_file = 0;
  for (auto &s : seq) {
    kseq_destroy(s);
    s = nullptr;
  }
}

void FastqSequenceReader::reserveNfiles(int n) {
  fp.resize(nfiles);
  l.resize(nfiles, 0);
  nl.resize(nfiles, 0);
  seq.resize(nfiles, nullptr);
}

// returns true if there is more left to read from the files
bool FastqSequenceReader::fetchSequences(char *buf, const int limit, std::vector<std::pair<const char *, int> > &seqs,
                                         std::vector<std::pair<const char *, int> > &names,
                                         std::vector<std::pair<const char *, int> > &quals,
                                         std::vector<uint32_t>& flags,
                                         int& read_id,
                                         bool full,
                                         bool comments) {
  
  std::string line;
  readbatch_id += 1; // increase the batch id
  read_id = readbatch_id; // copy now because we are inside a lock
  seqs.clear();
  if (full) {
    names.clear();
    quals.clear();
  }
  flags.clear();
  
  int bufpos = 0;
  int count = 0; // for interleaving
  int pad = nfiles;
  while (true) {
    if (!state) { // should we open a file
      if (current_file >= files.size()) {
        // nothing left
        return false;
      } else {
        // close the current files
        for (auto &f : fp) {
          if (f) {
            gzclose(f);
          }
        }
        
        // open the next one
        for (int i = 0; i < nfiles; i++) {
          fp[i] = files[0] == "-" && nfiles == 1 && files.size() == 1 ? gzdopen(fileno(stdin), "r") : gzopen(files[current_file+i].c_str(), "r");
          seq[i] = kseq_init(fp[i]);
          l[i] = kseq_read(seq[i]);
          
        }
        current_file+=nfiles;
        state = true; 
      }
    }
    // the file is open and we have read into seq1 and seq2
    bool all_l = true;
    int bufadd = nfiles;
    for (int i = 0; i < nfiles; i++) {
      all_l = all_l && l[i] >= 0;
      bufadd += l[i]; // includes seq
    }
    if (all_l) {      
      // fits into the buffer
      if (full) {
        for (int i = 0; i < nfiles; i++) {
          nl[i] = seq[i]->name.l + (comments ? seq[i]->comment.l+1 : 0);
          bufadd += nl[i]; // includes name
          //bufadd += l[i]; // include qual
        }
        bufadd += 2*pad;
      }
      
      if (bufpos+bufadd< limit) {
        if (interleave_nfiles != 0) { // Hack to allow interleaving
          // assert(nfiles == 1);
          if (bufpos+bufadd >= limit-262144 && count % interleave_nfiles == 0) {
            return true;
          }
          count++;
        }
        
        for (int i = 0; i < nfiles; i++) {
          char *pi = buf + bufpos;
          memcpy(pi, seq[i]->seq.s, l[i]+1);
          bufpos += l[i]+1;
          seqs.emplace_back(pi,l[i]);
          
          if (full && !comments) {
            //pi = buf + bufpos;
            //memcpy(pi, seq[i]->qual.s,l[i]+1);
            //bufpos += l[i]+1;
            //quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            memcpy(pi, seq[i]->name.s, nl[i]+1);
            bufpos += nl[i]+1;
            names.emplace_back(pi, nl[i]);
          } else if (full && comments) {
            //pi = buf + bufpos;
            //memcpy(pi, seq[i]->qual.s,l[i]+1);
            //bufpos += l[i]+1;
            //quals.emplace_back(pi,l[i]);
            pi = buf + bufpos;
            memcpy(pi, seq[i]->name.s, (nl[i]-(seq[i]->comment.l+1)));
            names.emplace_back(pi, nl[i]);
            bufpos += (nl[i]-(seq[i]->comment.l+1));
            pi = buf + bufpos;
            const char* blank_space = " ";
            memcpy(pi, blank_space, 1);
            bufpos += 1;
            pi = buf + bufpos;
            memcpy(pi, seq[i]->comment.s, seq[i]->comment.l+1);
            bufpos += seq[i]->comment.l+1;
          }
        }
        
        numreads++;
        flags.push_back((current_file-nfiles) / nfiles); // flags.push_back(numreads-1);
      } else {
        if (interleave_nfiles != 0) {
          std::cerr << "Error: There was an error processing interleaved FASTQ input. Exiting..." << std::endl;
          exit(1);
        }
        return true; // read it next time
      }
      
      // read for the next one
      for (int i = 0; i < nfiles; i++) {
        l[i] = kseq_read(seq[i]);
      }        
    } else {
      state = false; // haven't opened file yet
    }
  }
}

// move constructor

FastqSequenceReader::FastqSequenceReader(FastqSequenceReader&& o) :
  nfiles(o.nfiles),
  numreads(o.numreads),
  fp(std::move(o.fp)),
  l(std::move(o.l)),
  nl(std::move(o.nl)),
  files(std::move(o.files)),
  current_file(o.current_file),
  seq(std::move(o.seq)) {
  
  o.fp.resize(nfiles);
  o.l.resize(nfiles, 0);
  o.nl.resize(nfiles, 0);
  o.seq.resize(nfiles, nullptr);
  o.state = false;
}

