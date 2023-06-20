#include <algorithm>
#include <random>
#include <sstream>
#include <ctype.h>
#include <unordered_set>
#include <functional>
#include "common.h"
#include "KmerIndex.h"
#include <iostream>
#include <unordered_map>
#include <string>
#include "ColoredCDBG.hpp"


// other helper functions
// pre: u is sorted
bool isUnique(const std::vector<int>& u) {
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] == u[j]) {
      return false;
    }
  }
  return true;
}

std::vector<int> unique(const std::vector<int>& u) {
  std::vector<int> v;
  v.reserve(u.size());
  v.push_back(u[0]);
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] != u[j]) {
      v.push_back(u[j]);
    }
  }
  return v;
}

const char Dna(int i) {
  static const char *dna = "ACGT";
  return dna[i & 0x03];
}

int hamming(const char *a, const char *b) {
  int h = 0;
  while (*a != 0 && *b != 0) {
    if (*a != *b) {
      h++;
    }
    a++;
    b++;
  }
  return h;
}

std::string generate_tmp_file(std::string seed) {
  std::string base = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  std::string tmp_file = ".kure.";
  srand((unsigned int)std::hash<std::string>{}(seed));
  int pos;
  while(tmp_file.length() < 32) {
    pos = ((rand() % (base.size() - 1)));
    tmp_file += base.substr(pos, 1);
  }
  return tmp_file;
}

void KmerIndex::BuildDistinguishingGraph(const ProgramOptions& opt, std::ofstream& out) {
  k = opt.k;
  std::cerr << "[build] k-mer length: " << k << std::endl;
  size_t ncolors = 0;
  std::string out_file = opt.distinguish_output_fasta;
  std::vector<std::string> tmp_files;
  for (auto& fasta : opt.transfasta) {
    std::cerr << "[build] loading fasta file " << fasta
              << std::endl;
    tmp_files.push_back(generate_tmp_file(opt.distinguish_output_fasta + fasta));
  }
  
  std::vector<std::ofstream*> ofs; // Store pointers to circumvent certain compiler bugs where ofstream is non-movable
  for (auto tmp_file : tmp_files) ofs.push_back(new std::ofstream(tmp_file));
  num_trans = 0;
  
  // read fasta file using kseq
  gzFile fp = 0;
  kseq_t *seq;
  int l = 0;
  std::mt19937 gen(42);
  int countNonNucl = 0;
  int countTrim = 0;
  int countUNuc = 0;
  
  int i = 0;
  for (auto& fasta : opt.transfasta) {
    fp = opt.transfasta.size() == 1 && opt.transfasta[0] == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fasta.c_str(), "r");
    seq = kseq_init(fp);
    while (true) {
      l = kseq_read(seq);
      if (l <= 0) {
        break;
      }
      int trimNonNuclStart = 0;
      int trimNonNuclEnd = 0;
      int runningValidNuclLength = 0;
      bool finishTrimStart = false;
      std::string str = seq->seq.s;
      auto n = str.size();
      for (auto i = 0; i < n; i++) {
        char c = str[i];
        c = ::toupper(c);
        if (c=='U') {
          str[i] = 'T';
          countUNuc++;
        }
        if (c !='A' && c != 'C' && c != 'G' && c != 'T' && c != 'U') {
          countNonNucl++;
          if (runningValidNuclLength >= k && finishTrimStart) trimNonNuclEnd = i-1; // We trim the sequence end from the last valid k-mer onward
          runningValidNuclLength = 0;
        } else { // Valid nucleotide
          runningValidNuclLength++;
          if (!finishTrimStart) {
            if (runningValidNuclLength >= k) {
              finishTrimStart = true;
              trimNonNuclStart = i; // We trim the sequence beginning until we encounter k valid nucleotides (first valid k-mer)
            }
          }
          if (runningValidNuclLength >= k && finishTrimStart) trimNonNuclEnd = i; // We trim the sequence end from the last valid k-mer onward
        }
      }
      std::transform(str.begin(), str.end(),str.begin(), ::toupper);
      str = (trimNonNuclEnd == 0 ? str.substr(trimNonNuclStart) : str.substr(trimNonNuclStart,trimNonNuclEnd+1-trimNonNuclStart));
      countTrim += n - str.length();
      if (str.length() >= k) {
        *(ofs[i]) << ">" << num_trans++ << "\n" << str << std::endl;
      }
      //target_lens_.push_back(seq->seq.l);
      //std::string name(seq->name.s);
      //size_t p = name.find(' ');
      //if (p != std::string::npos) {
      //  name = name.substr(0,p);
      //}
    }
    
    gzclose(fp);
    fp=0;
    i++;
  }
  
  for (auto& of : ofs) (*of).close(); // Close files now that we've outputted everything
  for (auto& of : ofs) delete of; // Free pointer memory
  
  if (countNonNucl > 0) {
    std::cerr << "[build] warning: counted " << countNonNucl << " non-ACGUT characters in the input sequence" << std::endl;
  }
  
  if (countTrim > 0) {
    std::cerr << "[build] warning: trimmed " << countTrim << " characters from ends of input sequences" << std::endl;
  }
  
  if (countUNuc > 0) {
    std::cerr << "[build] warning: replaced " << countUNuc << " U characters with Ts" << std::endl;
  }
  
  CCDBG_Build_opt c_opt;
  c_opt.k = k;
  c_opt.nb_threads = opt.threads;
  c_opt.build = true;
  c_opt.clipTips = false;
  c_opt.deleteIsolated = false;
  c_opt.verbose = true;
  c_opt.filename_ref_in = tmp_files;
  
  if (opt.g > 0) { // If minimizer length supplied, override the default
    c_opt.g = opt.g;
  } else { // Define minimizer length defaults
    int g = k-8;
    if (k <= 13) {
      g = k-2;
    } else if (k <= 17) {
      g = k-4;
    } else if (k <= 19) {
      g = k-6;
    }
    c_opt.g = g;
  }
  
  std::cerr << "[build] Building colored graph" << std::endl;
  ColoredCDBG<void> ccdbg = ColoredCDBG<void>(k, c_opt.g);
  ccdbg.buildGraph(c_opt);
  ccdbg.buildColors(c_opt);
  std::cerr << "[build] Extracting k-mers from graph" << std::endl;
  std::ofstream of(out_file); // Write color contigs into another file
  size_t max_threads_read = opt.threads;
  std::vector<std::vector<std::pair<const UnitigColors*, const UnitigMap<DataAccessor<void>, DataStorage<void>, false> > > > unitigs_v(max_threads_read);
  size_t n = 0;
  const size_t thresh_size = 50000; // Max number of unitigs across all threads
  std::mutex mutex_unitigs; // Lock for multithreading writing output FASTA file
  std::vector<std::thread> workers; // Worker threads
  uint32_t rb = std::max(opt.distinguish_range_begin,0); // range begin filter
  uint32_t re = opt.distinguish_range_end == 0 ? rb : std::max(opt.distinguish_range_end,0); // range end filter
  if (rb == 0 && re == 0) re = std::numeric_limits<uint32_t>::max();
  for (const auto& unitig : ccdbg) {
    const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
    const UnitigMap<DataAccessor<void>, DataStorage<void>, false> unitig_ = unitig;
    unitigs_v[n % unitigs_v.size()].push_back(std::make_pair(uc, unitig_));
    n++;
    if (unitigs_v[unitigs_v.size()-1].size() >= thresh_size || n >= ccdbg.size()) {
      for (size_t u_i = 0; u_i < unitigs_v.size(); u_i++) {
        workers.emplace_back(
          [&, u_i] {
            std::ostringstream oss;
            for (auto unitig_x : unitigs_v[u_i]) {
              auto uc = unitig_x.first;
              auto& unitig = (unitig_x.second);
              UnitigColors::const_iterator it_uc = uc->begin(unitig);
              UnitigColors::const_iterator it_uc_end = uc->end();
              std::map<int, std::set<int>> k_map;
              for (; it_uc != it_uc_end; ++it_uc) {
                int color = it_uc.getColorID();
                k_map[color].insert(it_uc.getKmerPosition());
                // DEBUG:
                // std::cout << color << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).rep().toString() << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).toString() << " " << it_uc.getKmerPosition() << " " << unitig.strand << std::endl;
              }
              std::set<int> positions_to_remove;
              if (!opt.distinguish_all_but_one_color && !opt.distinguish_union) {
                int i_ = 0;
                for (const auto& k_elem : k_map) {
                  int j_ = 0;
                  for (const auto& k_elem2 : k_map) {
                    if (j_ > i_ && k_elem.first != k_elem2.first) {
                      std::set<int> intersect;
                      std::set<int> set_result;
                      std::set_intersection(k_elem.second.begin(), k_elem.second.end(), k_elem2.second.begin(), k_elem2.second.end(), std::inserter(intersect, intersect.begin())); //if (k_elem2.second.count(k_elem1.second)) // check if set intersection with k_elem2
                      std::set_union(positions_to_remove.begin(), positions_to_remove.end(), intersect.begin(), intersect.end(), std::inserter(set_result,set_result.begin()));
                      positions_to_remove = std::move(set_result);
                    }
                    j_++;
                  }
                  i_++;
                }
              } else if (!opt.distinguish_union) {
                int i_ = 0;
                if (k_map.size() == tmp_files.size()) {
                  for (const auto& k_elem : k_map) {
                    i_++;
                    if (positions_to_remove.size() == 0) {
                      positions_to_remove = k_elem.second;
                    } else {
                      std::set<int> set_result;
                      std::set_intersection(positions_to_remove.begin(), positions_to_remove.end(), k_elem.second.begin(), k_elem.second.end(), std::inserter(set_result,set_result.begin()));
                      positions_to_remove = std::move(set_result);
                    }
                  }
                }
              }
              for (const auto& k_elem : k_map) {
                int curr_pos = -1;
                std::string colored_contig = "";
                auto color = k_elem.first;
                //std::string contig_metadata = " :" + unitig.dist + "," + unitig.len + "," + unitig.size + "," + unitig.strand;
                for (const auto &pos : k_elem.second) {
                  if (!positions_to_remove.count(pos)) {
                    std::string km = unitig.getUnitigKmer(pos).toString();
                    if (curr_pos == -1) { // How to correspond color?
                      colored_contig = km;
                    } else if (pos == curr_pos+1) {
                      colored_contig += km[km.length()-1];
                    } else {
                      if (colored_contig.length() >= rb && colored_contig.length() <= re) oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n";
                      colored_contig = km;
                    }
                    curr_pos = pos;
                  }
                }
                if (colored_contig != "") {
                  if (colored_contig.length() >= rb && colored_contig.length() <= re) oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n";
                }
              }
            }
            std::unique_lock<std::mutex> lock(mutex_unitigs);
            of << oss.str();
          }
        );
      }
      for (auto& t : workers) t.join();
      workers.clear();
      for (size_t u_i = 0; u_i < unitigs_v.size(); u_i++) {
        unitigs_v[u_i].clear();
      }
    }
  }
  ccdbg.clear(); // Free memory associated with the colored compact dBG
  ncolors = tmp_files.size(); // Record the number of "colors"
  for (auto tmp_file : tmp_files) std::remove(tmp_file.c_str()); // Remove temp files needed to make colored graph
  of.close();
  tmp_files.clear();
  ofs.clear();
}

