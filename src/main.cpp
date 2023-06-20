#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>
#include <time.h>
#include <algorithm>
#include <limits>

#include <cstdio>

#include "common.h"
//#include "ProcessReads.h"
#include "KmerIndex.h"
#include <CompactedDBG.hpp>

//#define ERROR_STR "\033[1mError:\033[0m"
#define ERROR_STR "Error:"

using namespace std;

int my_mkdir(const char *path, mode_t mode) {
  #ifdef _WIN64
  return mkdir(path);
  #else
  return mkdir(path,mode);
  #endif
}

bool checkFileExists(std::string fn) {
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0;
}


void ParseOptionsDistinguish(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  const char *opt_string = "o:k:m:t:z:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"all", no_argument, &distinguish_all_flag, 1},
    {"all-but-one", no_argument, &distinguish_all_but_one_flag, 1},
    // short args
    {"output", required_argument, 0, 'o'},
    {"kmer-size", required_argument, 0, 'k'},
    {"min-size", required_argument, 0, 'm'},
    {"threads", required_argument, 0, 't'},
    {"distinguish-range", required_argument, 0, 'z'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 'o': {
      opt.distinguish_output_fasta = optarg;
      break;
    }
    case 'k': {
      stringstream(optarg) >> opt.k;
      break;
    }
    case 'm': {
      stringstream(optarg) >> opt.g;
      break;
    }
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'z': {
      std::string range_input_str;
      stringstream(optarg) >> range_input_str;
      stringstream ss(range_input_str);
      std::string range_val;
      int i = 0;
      char delimeter = '-';
      if (range_input_str.find(',') < range_input_str.length()) {
        delimeter = ','; // If string contains commas, use commas as delimeter
      }
      while (std::getline(ss, range_val, delimeter)) { 
        if (i == 0) opt.distinguish_range_begin = std::atoi(range_val.c_str());
        if (i == 1) opt.distinguish_range_end = std::atoi(range_val.c_str());
        i++;
      }
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
  
  if (distinguish_all_flag) {
    opt.distinguish_union = true;
  }
  if (distinguish_all_but_one_flag) {
    opt.distinguish_all_but_one_color = true;
  }

  for (int i = optind; i < argc; i++) {
    opt.transfasta.push_back(argv[i]);
  }
}

bool CheckOptionsDistinguish(ProgramOptions& opt) {

  bool ret = true;
  
  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
      opt.threads = n;
    }
  }

  if (opt.k <= 1 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error: invalid k-mer length " << opt.k << ", minimum is 3 and maximum is " << (MAX_KMER_SIZE -1) << endl;
    ret = false;
  }

  if (opt.k % 2 == 0) {
    cerr << "Error: k needs to be an odd number" << endl;
    ret = false;
  }

  if (opt.transfasta.empty()) {
    cerr << "Error: no FASTA files specified" << endl;
    ret = false;
  } else {

    for (auto& fasta : opt.transfasta) {
      // we want to generate the index, check k, index and transfasta
      struct stat stFileInfo;
      auto intStat = stat(fasta.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: FASTA file not found " << fasta << endl;
        ret = false;
      }
    }
  }

  if (opt.distinguish_output_fasta.empty()) {
    cerr << "Error: need to specify output FASTA file name" << endl;
    ret = false;
  }
  
  if (opt.distinguish_union && distinguish_all_but_one_color) {
    cerr << "Error: Cannot use multiple distinguish options" << endl;
    ret = false;
  }
  
  if (opt.transfasta.size() == 1) {
    opt.distinguishUseInput = true; // Use pre-generated FASTA input
  }

  if (opt.g != 0) {
    if (opt.g <= 2 || opt.g > opt.k - 2) {
      cerr << "Error: invalid minimizer size " << opt.g << ", minimum is 3 and maximum is k - 2" << endl;
      ret = false;
    }
  }

  return ret;
}

void PrintCite() {
  cout << "When using this program in your research, please cite" << endl << endl
       << "  [todo]" << endl
       << endl;
}

void PrintVersion() {
  cout << "kure, version " << KURE_VERSION << endl;
}

void usage() {
  cout << "kure " << KURE_VERSION << endl << endl
       << "Usage: kure <CMD> [arguments] .." << endl << endl
       << "Where <CMD> can be one of:" << endl << endl
       << "    distinguish   Extracts distinguishing contigs from FASTA/FASTQ files "<< endl
       << "    version       Prints version information" << endl
       << "    cite          Prints citation information" << endl << endl
       << "Running kure <CMD> without arguments prints usage information for <CMD>"<< endl << endl;
}

void usageDistinguish() {
  cout << "kure " << KURE_VERSION << endl
       << "Extracts distinguishing contigs from FASTA/FASTQ files" << endl << endl
       << "Usage: kure distinguish [arguments] FASTA-files" << endl << endl
       << "Required argument:" << endl
       << "-o, --output=STRING        Filename for the output FASTA" << endl << endl
       << "Optional argument:" << endl
       << "-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: " << (MAX_KMER_SIZE-1) << ")" << endl
       << "-t, --threads=INT           Number of threads to use (default: 1)" << endl
       << "-m, --min-size=INT          Length of minimizers (default: automatically chosen)" << endl
       << endl;

}

std::string argv_to_string(int argc, char *argv[]) {
  std::string res;
  for (int i = 0; i < argc; ++i) {
    res += argv[i];
    if (i + 1 < argc) {
      res += " ";
    }
  }

  return res;
}

std::string get_local_time() {
  time_t rawtime;
  struct tm * timeinfo;

  time( &rawtime );
  timeinfo = localtime( &rawtime );
  std::string ret(asctime(timeinfo));

  // chomp off the newline
  return ret.substr(0, ret.size() - 1);
}

int main(int argc, char *argv[]) {
  std::cout.sync_with_stdio(false);
  setvbuf(stdout, NULL, _IOFBF, 1048576);

  if (argc < 2) {
    usage();
    exit(1);
  } else {
    auto start_time(get_local_time());
    ProgramOptions opt;
    string cmd(argv[1]);
    if (cmd == "version") {
      PrintVersion();
    } else if (cmd == "cite") {
      PrintCite();
    } else if (cmd == "distinguish") {
      cerr << endl;
      if (argc==2) {
        usageIndex();
        return 0;
      }
      ParseOptionsDistinguish(argc-1,argv+1,opt);
      if (!CheckOptionsDistinguish(opt)) {
        usageIndex();
        exit(1);
      } else {
        Kmer::set_k(opt.k);
        KmerIndex index(opt);
        std::ofstream out;
        out.open(opt.distinguish_output_fasta, std::ios::out | std::ios::binary);
        index.BuildDistinguishingGraph(opt, out);
      }
      cerr << endl;
    } else {
      cerr << "Error: invalid command " << cmd << endl;
      usage();
      exit(1);
    }

  }

  fflush(stdout);

  return 0;
}
