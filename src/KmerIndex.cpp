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

/*
#include <set>
#include <map>
#include <vector>
#include <thread>
#include <mutex>
*/


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
  std::string tmp_file = ".klue.";
  srand((unsigned int)std::hash<std::string>{}(seed));
  int pos;
  while(tmp_file.length() < 32) {
    pos = ((rand() % (base.size() - 1)));
    tmp_file += base.substr(pos, 1);
  }
  return tmp_file;
}

void KmerIndex::BuildReconstructionGraph(const ProgramOptions& opt) {
  // Read in FASTAs and output each color into a new separate temp file
  std::vector<std::string> transfasta = opt.transfasta;
  std::vector<std::string> tmp_files;
  std::vector<std::ofstream*> ofs; // Store pointers to circumvent certain compiler bugs where ofstream is non-movable
  gzFile fp = 0;
  kseq_t *seq;
  int l = 0;
  num_trans = 0;
  size_t range_discard = 0;
  uint32_t rb = std::max(opt.distinguish_range_begin,0); // range begin filter
  uint32_t re = opt.distinguish_range_end == 0 ? rb : std::max(opt.distinguish_range_end,0); // range end filter
  if (rb == 0 && re == 0) re = std::numeric_limits<uint32_t>::max();
  for (auto& fasta : transfasta) {
    std::cerr << "[build] Preparing FASTA file: " << fasta << std::endl;
    fp = transfasta.size() == 1 && transfasta[0] == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fasta.c_str(), "r");
    seq = kseq_init(fp);
    while (true) {
      l = kseq_read(seq);
      if (l <= 0) {
        break;
      }
      std::string name = seq->name.s;
      std::string str = seq->seq.s;
      int color;
      try {
        color = std::stoi(name);
      } catch (std::exception const & e) {
        std::cerr << "Error: Non-numerical name found in " << fasta << std::endl;
        exit(1);
      }
      if (color < 0 || color > 4096) {
        std::cerr << "Error: Invalid number name, " << std::to_string(color) << ", found in " << fasta << std::endl;
        exit(1);
      }
      if (color >= tmp_files.size()) {
        tmp_files.resize(color+1);
        ofs.resize(color+1);
        for (int i = 0; i < tmp_files.size(); i++) {
          if (tmp_files[i].empty()) {
            tmp_files[i] = generate_tmp_file(opt.distinguish_output_fasta + fasta + std::to_string(i));
            ofs[i] = new std::ofstream(tmp_files[i]); // Store pointers to circumvent certain compiler bugs where ofstream is non-movable
          }
        }
      }
      if (str.length() < k) {
        continue;
      }
      if (str.length() >= rb && str.length() <= re) {
        *(ofs[color]) << ">" << std::to_string(color) << "\n" << str << "\n";
        num_trans++;
      } else {
        range_discard++;
      }
    }
    gzclose(fp);
    fp=0;
  }
  for (auto& of : ofs) (*of).close(); // Close files now that we've outputted everything
  for (auto& of : ofs) delete of; // Free pointer memory
  if (range_discard > 0) {
    std::cerr << "[build] Number of input sequences filtered out due to length: " << range_discard << std::endl;
  }
  BuildDistinguishingGraph(opt, tmp_files, true); // TODO: modify to handle temporary files
}

struct SharedPositionInfo {
    std::set<int> positions;
};

std::mutex mtx;
bool test_distinguish_combinations = false; /* DEBUG: change to ProgramOptions */

void findSharedPositions(const std::vector<int>& colors,
    const std::map<int, std::set<int>>& k_map,
    std::map<std::vector<int>, SharedPositionInfo>& shared_positions,
    std::map<std::vector<int>, std::set<int>>& positions_to_remove_map,
    std::map<std::vector<int>, std::set<int>>& aggregated_positions,
    std::mutex& positions_mutex,
    std::mutex& aggregated_mutex) {

    std::set<int> common_positions = k_map.at(colors[0]);

    for (size_t i = 1; i < colors.size(); ++i) {
        std::set<int> temp;
        std::set_intersection(common_positions.begin(), common_positions.end(),
            k_map.at(colors[i]).begin(), k_map.at(colors[i]).end(),
            std::inserter(temp, temp.begin()));
        common_positions.swap(temp);
    }

    if (!common_positions.empty()) {
        std::lock_guard<std::mutex> lock(mtx);

        for (auto& entry : shared_positions) {
            const std::vector<int>& subset_colors = entry.first;
            if (subset_colors.size() <= colors.size() && std::includes(colors.begin(), colors.end(),
                subset_colors.begin(), subset_colors.end())) {
                SharedPositionInfo& info = entry.second;
                std::set<int>& positions = info.positions;
                std::set<int> excluded;
                for (int pos : common_positions) {
                    if (positions.erase(pos)) {
                        excluded.insert(pos);
                    }
                }
                if (!excluded.empty()) {
                    std::lock_guard<std::mutex> positions_lock(positions_mutex);
                    positions_to_remove_map[subset_colors].insert(excluded.begin(), excluded.end());
                }
            }
        }

        SharedPositionInfo& info = shared_positions[colors];
        std::set<int>& positions = info.positions;
        std::set<int> unique_positions;
        std::set_difference(common_positions.begin(), common_positions.end(),
            positions.begin(), positions.end(),
            std::inserter(unique_positions, unique_positions.begin()));
        positions.insert(unique_positions.begin(), unique_positions.end());

        if (!aggregated_positions.count(colors)) {
            std::lock_guard<std::mutex> aggregated_lock(aggregated_mutex);
            aggregated_positions[colors] = common_positions;
        }
    }
}

void generateCombinations(const std::vector<int>& elements,
    int startIdx,
    int remaining,
    std::vector<int>& currentCombination,
    std::vector<std::vector<int>>& allCombinations) {

    if (remaining == 0) {
        allCombinations.emplace_back(currentCombination);
        return;
    }

    for (size_t i = startIdx; i <= elements.size() - remaining; ++i) {
        currentCombination.push_back(elements[i]);
        generateCombinations(elements, i + 1, remaining - 1, currentCombination, allCombinations);
        currentCombination.pop_back();
    }
}

// convert a vector of integers to a string
std::string to_string(const std::vector<int>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i < vec.size() - 1) {
            oss << "_";
        }
    }
    return oss.str();
}

void KmerIndex::BuildDistinguishingGraph(const ProgramOptions& opt, const std::vector<std::string>& transfasta, bool reconstruct) {
  k = opt.k;
  std::cerr << "[build] k-mer length: " << k << std::endl;
  size_t ncolors = 0;
  std::string out_file = opt.distinguish_output_fasta;
  std::vector<std::string> tmp_files;
  if (reconstruct) {
    tmp_files = transfasta;
  } else {
    for (auto& fasta : transfasta) {
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
    for (auto& fasta : transfasta) {
      if (opt.kmer_multiplicity[i] != 1) {
        i++;
        continue; // Don't do any processing or anything, use input file as-is
      }
      fp = transfasta.size() == 1 && transfasta[0] == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fasta.c_str(), "r");
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
                trimNonNuclStart = i+1-k; // We trim the sequence beginning until we encounter k valid nucleotides (first valid k-mer)
              }
            }
            if (runningValidNuclLength >= k && finishTrimStart) trimNonNuclEnd = i; // We trim the sequence end from the last valid k-mer onward
          }
        }
        std::transform(str.begin(), str.end(),str.begin(), ::toupper);
        str = (trimNonNuclEnd == 0 ? str.substr(trimNonNuclStart) : str.substr(trimNonNuclStart,trimNonNuclEnd+1-trimNonNuclStart));
        countTrim += n - str.length();
        if (str.length() >= k) {
          *(ofs[i]) << ">" << num_trans++ << "\n" << str << "\n";
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
    ofs.clear();
    
    if (countNonNucl > 0) {
      std::cerr << "[build] warning: counted " << countNonNucl << " non-ACGUT characters in the input sequence" << std::endl;
    }
    
    if (countTrim > 0) {
      std::cerr << "[build] warning: trimmed " << countTrim << " characters from ends of input sequences" << std::endl;
    }
    
    if (countUNuc > 0) {
      std::cerr << "[build] warning: replaced " << countUNuc << " U characters with Ts" << std::endl;
    }
  }
  
  CCDBG_Build_opt c_opt;
  c_opt.k = k;
  c_opt.nb_threads = opt.threads;
  c_opt.build = true;
  c_opt.clipTips = false;
  c_opt.deleteIsolated = false;
  c_opt.verbose = opt.verbose;
  for (int i = 0; i < tmp_files.size(); i++) {
    if (reconstruct || opt.kmer_multiplicity[i] == 1) {
      c_opt.filename_ref_in.push_back(tmp_files[i]);
    } else {
      c_opt.filename_seq_in.push_back(transfasta[i]);
    }
  }
  
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
  auto color_names = ccdbg.getColorNames();
  std::vector<int> color_map;
  color_map.resize(tmp_files.size());
  if (!reconstruct) {
    for (int i = 0; i < color_names.size(); i++) {
      auto color_name = color_names[i];
      for (int j = 0; j < tmp_files.size(); j++) {
        std::string fname = tmp_files[j];
        std::string real_fname = transfasta[j];
        if (opt.kmer_multiplicity[j] != 1) {
          fname = real_fname;
        }
        if (color_name == fname) {
          std::cerr << "        " << real_fname << ": " << j << " (multiplicity: " << opt.kmer_multiplicity[j] << ")" << std::endl;
          color_map[i] = j;
        }
      }
    }
  } else {
    for (int i = 0; i < color_names.size(); i++) {
      auto color_name = color_names[i];
      for (int j = 0; j < tmp_files.size(); j++) {
        if (color_name == tmp_files[j]) {
          color_map[i] = j;
        }
      }
    }
  }

  std::cerr << "[build] Extracting k-mers from graph" << std::endl;
  std::streambuf* buf = nullptr;
  std::ofstream of;
  if (!opt.stream_out) {
    of.open(out_file); // Write color contigs into another file
    buf = of.rdbuf();
  } else {
    buf = std::cout.rdbuf();
  }
  std::ostream o(buf);
  size_t max_threads_read = opt.threads;
  std::vector<std::vector<std::pair<const UnitigColors*, const UnitigMap<DataAccessor<void>, DataStorage<void>, false> > > > unitigs_v(max_threads_read);
  size_t n = 0;
  const size_t thresh_size = 50000; // Max number of unitigs across all threads
  std::mutex mutex_unitigs; // Lock for multithreading writing output FASTA file
  std::vector<std::thread> workers; // Worker threads
  uint32_t rb = std::max(opt.distinguish_range_begin,0); // range begin filter
  uint32_t re = opt.distinguish_range_end == 0 ? rb : std::max(opt.distinguish_range_end,0); // range end filter
  if (rb == 0 && re == 0) re = std::numeric_limits<uint32_t>::max();
  int range_discard = 0;
  int num_written = 0;
  // TODO: Reconstruct below
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
            int _num_written = 0;
            int _range_discard = 0;
            for (auto unitig_x : unitigs_v[u_i]) {
              auto uc = unitig_x.first;
              auto& unitig = (unitig_x.second);
              UnitigColors::const_iterator it_uc = uc->begin(unitig);
              UnitigColors::const_iterator it_uc_end = uc->end();
              std::map<int, std::set<int>> k_map;
              for (; it_uc != it_uc_end; ++it_uc) {
                int color = color_map[it_uc.getColorID()];
                k_map[color].insert(it_uc.getKmerPosition());
                // DEBUG:
                // std::cout << color << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).rep().toString() << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).toString() << " " << it_uc.getKmerPosition() << " " << unitig.strand << std::endl;
              }
              std::set<int> positions_to_remove;
              std::map<std::vector<int>, std::set<int>> positions_to_remove_map;
              std::map<std::vector<int>, std::set<int>> aggregated_positions;
              std::map<std::vector<int>, std::set<int>> new_k_map;
              if (!opt.distinguish_all_but_one_color && !opt.distinguish_union && !test_distinguish_combinations) {
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
              else if (test_distinguish_combinations) {
                  int i_ = 0;
                  std::vector<std::vector<int>> allCombinations;
                  std::vector<int> elements;
                  for (const auto& k_elem : k_map) {
                      elements.push_back(k_elem.first);
                  }
                  for (size_t i = 1; i <= elements.size(); ++i) {
                      std::vector<int> currentCombination;
                      generateCombinations(elements, 0, static_cast<int>(i), currentCombination, allCombinations);
                  }
                  std::map<std::vector<int>, SharedPositionInfo> shared_positions;
                  std::vector<std::thread> threads;
                  threads.reserve(allCombinations.size());
                  std::mutex positions_mutex;
                  std::mutex aggregated_mutex;
                  for (const auto& colors : allCombinations) {
                      threads.emplace_back(findSharedPositions, std::cref(colors), std::cref(k_map),
                          std::ref(shared_positions), std::ref(positions_to_remove_map),
                          std::ref(aggregated_positions), std::ref(positions_mutex), std::ref(aggregated_mutex));
                  }
                  for (auto& thread : threads) {
                      thread.join();
                  }
              }
              // k_map modification
              for (const auto& entry : k_map) {
                  std::vector<int> key_vector = { entry.first };
                  new_k_map[key_vector] = entry.second;
              }
              k_map.clear();
              k_map.~map();
              if (test_distinguish_combinations) { new_k_map = aggregated_positions; }
              for (const auto& k_elem : new_k_map) {
                int curr_pos = -1;
                std::string colored_contig = "";
                auto color = k_elem.first;
                if (test_distinguish_combinations) { positions_to_remove = positions_to_remove_map[combination]; }
                //std::string contig_metadata = " :" + unitig.dist + "," + unitig.len + "," + unitig.size + "," + unitig.strand;
                for (const auto &pos : k_elem.second) {
                  if (!positions_to_remove.count(pos)) {
                    std::string km = unitig.getUnitigKmer(pos).toString();
                    if (curr_pos == -1) { // How to correspond color?
                      colored_contig = km;
                    } else if (pos == curr_pos+1) {
                      colored_contig += km[km.length()-1];
                    } else {
                      if (colored_contig.length() >= rb && colored_contig.length() <= re) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                      else _range_discard++;
                      colored_contig = km;
                    }
                    curr_pos = pos;
                  }
                }
                if (colored_contig != "") {
                  if (colored_contig.length() >= rb && colored_contig.length() <= re) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                  else _range_discard++;
                }
              }
            }
            std::unique_lock<std::mutex> lock(mutex_unitigs);
            o << oss.str();
            num_written += _num_written;
            range_discard += _range_discard;
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
  o.flush();
  ccdbg.clear(); // Free memory associated with the colored compact dBG
  ncolors = tmp_files.size(); // Record the number of "colors"
  for (auto tmp_file : tmp_files) std::remove(tmp_file.c_str()); // Remove temp files needed to make colored graph
  if (!opt.stream_out) {
    of.close();
  }
  tmp_files.clear();
  if (opt.verbose) {
    if (range_discard > 0) std::cerr << "[build] Number of output sequences filtered out due to length: " << range_discard << std::endl;
    std::cerr << "[build] Number of output sequences written: " << num_written << std::endl;
  }
}

