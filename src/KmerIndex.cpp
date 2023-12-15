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
#include <ExpressionParser.h>
#include <iomanip>
#include <queue>
#include <algorithm>

// other helper functions
// pre: u is sorted
bool isUnique(const std::vector<int>& u) {
    for (int j = 1; j < u.size(); j++) {
        if (u[j - 1] == u[j]) {
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
        if (u[j - 1] != u[j]) {
            v.push_back(u[j]);
        }
    }
    return v;
}

const char Dna(int i) {
    static const char* dna = "ACGT";
    return dna[i & 0x03];
}

int hamming(const char* a, const char* b) {
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
    while (tmp_file.length() < 32) {
        pos = ((rand() % (base.size() - 1)));
        tmp_file += base.substr(pos, 1);
    }
    return tmp_file;
}

// Perform DFS traversal starting from unitig head until no neighbors are found or cycle completed
/*
std::string forwardDFS(ColoredCDBG<void>& ccdbg, UnitigMap<DataAccessor<void>, DataStorage<void>>& um, int k, int depth, int MAX_DEPTH, std::unordered_set<std::string>& visited) {
    std::string assembledSequence = "";

    std::string currentTail = um.getUnitigTail().toString();

    if (depth > MAX_DEPTH || um.isEmpty || visited.find(currentTail) != visited.end()) { return ""; }

    visited.insert(currentTail);

    std::vector<int> colors;
    UnitigColors::const_iterator it_uc_ = um.getData()->getUnitigColors(um)->begin(um);
    for (; it_uc_ != um.getData()->getUnitigColors(um)->end(); ++it_uc_) { colors.push_back(it_uc_.getColorID()); }

    for (const auto& c : colors) {
        for (char nextChar : "ACGT") {
            std::string next_kmer = currentTail.substr(1, k - 1) + nextChar;
            Kmer neighbor(next_kmer.c_str());
            UnitigMap<DataAccessor<void>, DataStorage<void>> neighbor_um = ccdbg.find(neighbor, false);

            if (!neighbor_um.isEmpty && std::find(visited.begin(), visited.end(), neighbor_um.getUnitigKmer(neighbor_um.dist).toString()) == visited.end()) {
                std::vector<int> neighbor_colors;
                UnitigColors::const_iterator neighbor_it_uc_ = neighbor_um.getData()->getUnitigColors(neighbor_um)->begin(neighbor_um);
                for (; neighbor_it_uc_ != neighbor_um.getData()->getUnitigColors(neighbor_um)->end(); ++neighbor_it_uc_) { neighbor_colors.push_back(neighbor_it_uc_.getColorID()); }
                if (std::find(neighbor_colors.begin(), neighbor_colors.end(), c) != neighbor_colors.end() &&
                    currentTail.substr(1, k - 1) == neighbor_um.getUnitigHead().toString().substr(0, k - 1)) {
                    visited.insert(neighbor_um.getUnitigTail().toString()); // Keep track of visited kmer strings
                    assembledSequence += forwardDFS(ccdbg, neighbor_um, k, depth + 1, MAX_DEPTH, visited);
                }
            }
        }
    }

    return assembledSequence;
}
*/

// Perform DFS traversal starting from unitig head until no neighbors are found or cycle completed
std::string forwardDFS(ColoredCDBG<void>& ccdbg, UnitigMap<DataAccessor<void>, DataStorage<void>>& um, int k, std::unordered_set<std::string>& visited) {
    std::string assembledSequence = "";
    std::string currentTail = um.getUnitigTail().toString();
    visited.insert(currentTail);

    std::vector<int> colors;
    UnitigColors::const_iterator it_uc_ = um.getData()->getUnitigColors(um)->begin(um);
    for (; it_uc_ != um.getData()->getUnitigColors(um)->end(); ++it_uc_) { colors.push_back(it_uc_.getColorID()); }

    for (const auto& c : colors) {
        for (char nextChar : "ACGT") {
            std::string next_kmer = currentTail.substr(1, k - 1) + nextChar;
            Kmer neighbor(next_kmer.c_str());
            UnitigMap<DataAccessor<void>, DataStorage<void>> neighbor_um = ccdbg.find(neighbor, false);

            if (!neighbor_um.isEmpty && std::find(visited.begin(), visited.end(), neighbor_um.getUnitigKmer(neighbor_um.dist).toString()) == visited.end()) {
                std::vector<int> neighbor_colors;
                UnitigColors::const_iterator neighbor_it_uc_ = neighbor_um.getData()->getUnitigColors(neighbor_um)->begin(neighbor_um);
                for (; neighbor_it_uc_ != neighbor_um.getData()->getUnitigColors(neighbor_um)->end(); ++neighbor_it_uc_) { neighbor_colors.push_back(neighbor_it_uc_.getColorID()); }
                if (std::find(neighbor_colors.begin(), neighbor_colors.end(), c) != neighbor_colors.end() &&
                    currentTail.substr(1, k - 1) == neighbor_um.getUnitigHead().toString().substr(0, k - 1)) {
                    visited.insert(neighbor_um.getUnitigHead().toString()); // Keep track of visited kmer strings or tail?
                    assembledSequence += forwardDFS(ccdbg, neighbor_um, k, visited);
                }
            }
        }
    }
    return assembledSequence;
}

// Updated forward DFS (void)
/*
void extendUnitig(const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current, const int& current_color, std::string current_str, const std::unordered_set<int>& superset_colors, const int& k, std::unordered_set<std::string>& visited) {
    if (visited.find(current.getUnitigHead().rep().toString()) != visited.end()) { return; }
    visited.insert(current.getUnitigHead().rep().toString());

    bool hasValidSuccessor = false; 
    for (const auto& next : current.getSuccessors()) {
        UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
        UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
        std::unordered_set<int> colors;
        for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
        for (const auto& cc : colors) {
            if (std::find(superset_colors.begin(), superset_colors.end(), cc) != superset_colors.end() && cc == current_color && next.strand == current.strand) {
                hasValidSuccessor = true;
                std::string output;
                if (next.getUnitigKmer(next.dist).rep() == next.getUnitigKmer(next.dist).twin()) { output = next.getUnitigKmer(next.dist).rep().toString(); }
                else { output = next.getUnitigKmer(next.dist).twin().toString(); }
                if (next.getUnitigHead().rep().toString().substr(0, k - next.len) != current.getUnitigTail().toString().substr(current.getUnitigTail().toString().size() - k)) { output = next.getUnitigKmer(next.dist).toString(); }
                std::string extended_str = current_str + output.substr(output.length() - next.len);
                extendUnitig(next, current_color, extended_str, superset_colors, k, visited);
            }
        }
    }
    if (!hasValidSuccessor) { std::cout << "extended unitig: " << current_str << "\n"; } // at leaf node
}
*/

// Updated forward DFS (string)
// Perform DFS traversal starting from unitig tail
std::string extendUnitig(const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current, const int& current_color, std::string current_str, const std::unordered_set<int>& superset_colors, const int& k, std::unordered_set<std::string>& visited) {
    if (visited.find(current.getUnitigHead().rep().toString()) != visited.end()) { return "";  }
    visited.insert(current.getUnitigHead().rep().toString());

    std::string result = "";
    bool hasValidSuccessor = false;
    for (const auto& next : current.getSuccessors()) {
        UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
        UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
        std::unordered_set<int> colors;
        for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
        for (const auto& cc : colors) {
            if (std::find(superset_colors.begin(), superset_colors.end(), cc) != superset_colors.end() && cc == current_color && next.strand == current.strand) {
                hasValidSuccessor = true;
                std::string output;
                if (next.getUnitigKmer(next.dist).rep() == next.getUnitigKmer(next.dist).twin()) { output = next.getUnitigKmer(next.dist).rep().toString(); }
                else { output = next.getUnitigKmer(next.dist).twin().toString(); }
                if (next.getUnitigHead().rep().toString().substr(0, k - next.len) != current.getUnitigTail().toString().substr(current.getUnitigTail().toString().size() - k)) { output = next.getUnitigKmer(next.dist).toString(); }
                std::string extended_str = current_str + output.substr(output.length() - next.len);
                result += extendUnitig(next, current_color, extended_str, superset_colors, k, visited);
            }
        }
    }
    if (!hasValidSuccessor) { result += current_str + " "; } // at leaf node
    return result;
}

// Begin set operations
// Split user-inputted set operation commands by " "
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Recursively compute set operations (union and intersection are nodes) that correspond to specific sets (leaves)
std::set<int> computeSetOperation(const Node* root, const std::map<int, std::set<int>>& k_map) {
    if (!root) { return {}; }
    if (root->value != 'U' && root->value != 'I' && root->value != '\\') {
        auto it = k_map.find(root->value - 'A'); // 'A' maps to 0, 'B' to 1, etc.
        return (it != k_map.end()) ? it->second : std::set<int>{};
    }
    auto leftSet = computeSetOperation(root->left, k_map);
    auto rightSet = computeSetOperation(root->right, k_map);
    std::set<int> result;
    if (root->value == 'U') {
        std::set_union(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == 'I') {
        std::set_intersection(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == '\\') {
        std::set_difference(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    return result;
}
// End set operations

std::string stringOpt(size_t n, std::string option) {
    std::string result;
    if (option == "default") {
        for (size_t i = 0; i < n; ++i) {
            result += 'A' + i;
            result += "\\(";
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    result += 'A' + j;
                    if (j != n - 1 && !(i == n - 1 && j == n - 2)) { result += 'U'; }
                }
            }
            result += ')';
            if (i != n - 1) { result += ' '; }
        }
    }
    else if (option == "all-but-one") {
        for (size_t i = 0; i < n; ++i) {
            result += 'A' + i;
            result += "\\(";
            for (size_t j = 0; j < n; ++j) {
                result += 'A' + j;
                if (j != n - 1) { result += 'I'; }
            }
            result += ')';
            if (i != n - 1) { result += ' '; }
        }
    }

    return result;
}

void KmerIndex::BuildReconstructionGraph(const ProgramOptions& opt) {
    // Read in FASTAs and output each color into a new separate temp file
    std::vector<std::string> transfasta = opt.transfasta;
    std::vector<std::string> tmp_files;
    std::vector<std::ofstream*> ofs; // Store pointers to circumvent certain compiler bugs where ofstream is non-movable
    gzFile fp = 0;
    kseq_t* seq;
    int l = 0;
    num_trans = 0;
    size_t range_discard = 0;
    uint32_t rb = std::max(opt.distinguish_range_begin, 0); // range begin filter
    uint32_t re = opt.distinguish_range_end == 0 ? rb : std::max(opt.distinguish_range_end, 0); // range end filter
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
            }
            catch (std::exception const& e) {
                std::cerr << "Error: Non-numerical name found in " << fasta << std::endl;
                exit(1);
            }
            if (color < 0 || color > 4096) {
                std::cerr << "Error: Invalid number name, " << std::to_string(color) << ", found in " << fasta << std::endl;
                exit(1);
            }
            if (color >= tmp_files.size()) {
                tmp_files.resize(color + 1);
                ofs.resize(color + 1);
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
            }
            else {
                range_discard++;
            }
        }
        gzclose(fp);
        fp = 0;
    }
    for (auto& of : ofs) (*of).close(); // Close files now that we've outputted everything
    for (auto& of : ofs) delete of; // Free pointer memory
    if (range_discard > 0) {
        std::cerr << "[build] Number of input sequences filtered out due to length: " << range_discard << std::endl;
    }
    BuildDistinguishingGraph(opt, tmp_files, true); // TODO: modify to handle temporary files
}

void KmerIndex::BuildDistinguishingGraph(const ProgramOptions& opt, const std::vector<std::string>& transfasta, bool reconstruct) {
    k = opt.k;
    std::cerr << "[build] k-mer length: " << k << std::endl;
    size_t ncolors = 0;
    std::string out_file = opt.distinguish_output_fasta;
    std::vector<std::string> tmp_files;
    if (reconstruct) {
        tmp_files = transfasta;
    }
    else {
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
        kseq_t* seq;
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
                    if (c == 'U') {
                        str[i] = 'T';
                        countUNuc++;
                    }
                    if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'U') {
                        countNonNucl++;
                        if (runningValidNuclLength >= k && finishTrimStart) trimNonNuclEnd = i - 1; // We trim the sequence end from the last valid k-mer onward
                        runningValidNuclLength = 0;
                    }
                    else { // Valid nucleotide
                        runningValidNuclLength++;
                        if (!finishTrimStart) {
                            if (runningValidNuclLength >= k) {
                                finishTrimStart = true;
                                trimNonNuclStart = i + 1 - k; // We trim the sequence beginning until we encounter k valid nucleotides (first valid k-mer)
                            }
                        }
                        if (runningValidNuclLength >= k && finishTrimStart) trimNonNuclEnd = i; // We trim the sequence end from the last valid k-mer onward
                    }
                }
                std::transform(str.begin(), str.end(), str.begin(), ::toupper);
                str = (trimNonNuclEnd == 0 ? str.substr(trimNonNuclStart) : str.substr(trimNonNuclStart, trimNonNuclEnd + 1 - trimNonNuclStart));
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
            fp = 0;
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
        }
        else {
            c_opt.filename_seq_in.push_back(transfasta[i]);
        }
    }

    if (opt.g > 0) { // If minimizer length supplied, override the default
        c_opt.g = opt.g;
    }
    else { // Define minimizer length defaults
        int g = k - 8;
        if (k <= 13) {
            g = k - 2;
        }
        else if (k <= 17) {
            g = k - 4;
        }
        else if (k <= 19) {
            g = k - 6;
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
    }
    else {
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
    }
    else {
        buf = std::cout.rdbuf();
    }
    std::ostream o(buf);
    size_t max_threads_read = opt.threads;
    std::vector<std::vector<std::pair<const UnitigColors*, const UnitigMap<DataAccessor<void>, DataStorage<void>, false> > > > unitigs_v(max_threads_read);
    size_t n = 0;
    const size_t thresh_size = 50000; // Max number of unitigs across all threads
    std::mutex mutex_unitigs; // Lock for multithreading writing output FASTA file
    std::vector<std::thread> workers; // Worker threads
    uint32_t rb = std::max(opt.distinguish_range_begin, 0); // range begin filter
    uint32_t re = opt.distinguish_range_end == 0 ? rb : std::max(opt.distinguish_range_end, 0); // range end filter
    if (rb == 0 && re == 0) re = std::numeric_limits<uint32_t>::max();
    int range_discard = 0;
    int num_written = 0;
    std::string input_str = opt.input_set_operations;
    bool perform_set_operations = false;
    auto expressions = split(input_str, ' ');
    std::map<std::string, int> expr_to_int;
    if (!input_str.empty()) {
        perform_set_operations = true;
        ExpressionParser parser(input_str);  // create parser instance
        auto tokens = parser.tokenize(input_str);
        for (int i = 0; i < expressions.size(); ++i) { expr_to_int[expressions[i]] = i; }
        int maxWidth = 0;
        for (const auto& pair : expr_to_int) { maxWidth = std::max(maxWidth, static_cast<int>(pair.first.length())); }
        for (const auto& pair : expr_to_int) { std::cout << std::setw(maxWidth + 8) << pair.first << ": " << pair.second << "\n"; }
    }
    // TODO: Reconstruct below
    for (const auto& unitig : ccdbg) { // Iterate through all the unitigs in the de bruijn graph
        const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
        const UnitigMap<DataAccessor<void>, DataStorage<void>, false> unitig_ = unitig;
        unitigs_v[n % unitigs_v.size()].push_back(std::make_pair(uc, unitig_)); // unitigs_v = vector of vectors of unitigs (b/c each thread contains a vector of unitigs)
        n++;
        if (unitigs_v[unitigs_v.size() - 1].size() >= thresh_size || n >= ccdbg.size()) { // If we're ready to start processing the unitigs using a series of workers
            for (size_t u_i = 0; u_i < unitigs_v.size(); u_i++) { // u_i = specific batch of unitigs that a worker will act on
                workers.emplace_back(
                    [&, u_i] {
                        std::ostringstream oss;
                        int _num_written = 0;
                        int _range_discard = 0;
                        for (auto unitig_x : unitigs_v[u_i]) { // Go through unitigs in the batch labeled u_i
                            auto uc = unitig_x.first;
                            auto& unitig = (unitig_x.second);
                            UnitigColors::const_iterator it_uc = uc->begin(unitig);
                            UnitigColors::const_iterator it_uc_end = uc->end();
                            std::map<int, std::set<int>> k_map; // key = color; value = list of positions (i.e. k-mers) along the current unitig (note: a k-mer is a position along a unitig)
                            std::unordered_set<int> superset_colors;
                            for (; it_uc != it_uc_end; ++it_uc) {
                                superset_colors.insert(it_uc.getColorID());
                                int color = color_map[it_uc.getColorID()];
                                k_map[color].insert(it_uc.getKmerPosition());
                                // DEBUG:
                                // std::cout << color << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).rep().toString() << " " << unitig.getUnitigKmer(it_uc.getKmerPosition()).toString() << " " << it_uc.getKmerPosition() << " " << unitig.strand << std::endl;
                            }
                            // begin extend 
                            if (opt.extend) {
                                std::unordered_set<std::string> visited;
                                k_map.clear();
                                for (const auto& c : superset_colors) {
                                    // extendUnitig(unitig, c, unitig.getUnitigTail().toString(), superset_colors, k, visited); // void
                                    std::string result = extendUnitig(unitig, c, unitig.getUnitigTail().toString(), superset_colors, k, visited);
                                    if (!result.empty()) {
                                        for (const auto& s : split(result, ' ')) {
                                            // [DEBUG] print extended contigs (all)
                                            // std::cout << ">" << c << "\n" << s << "\n";
                                            for (auto iter = s.begin(); iter != s.end(); ++iter) {
                                                int pos = std::distance(s.begin(), iter);
                                                k_map[c].insert(pos);
                                            }
                                        }
                                    }
                                }
                                // distinguish (default : find k-mer positions exlusive to each color)
                                /*
                                std::string default_str = stringOpt(tmp_files.size(), "default"); // for n=8, generate "A\(AIBI...IH) B\(AIBI...IH) ... H\(AIBI...IH)"
                                auto all_expr = split(default_str, ' ');
                                for (const auto& expr : all_expr) {
                                    ExpressionParser expr_parser(expr);
                                    Node* root = expr_parser.parse();
                                    std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                    // write out what remains among the contigs
                                    for (const auto& k_elem : k_map) {
                                        int curr_pos = -1;
                                        auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                        std::string colored_contig = "";
                                        for (const auto& pos : k_elem.second) {
                                            if (set_operation_result.count(pos)) {
                                                std::string km = unitig.getUnitigKmer(pos).toString();
                                                if (curr_pos == -1) { colored_contig = km; }
                                                else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                else {
                                                    if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                    else _range_discard++;
                                                    colored_contig = km;
                                                }
                                                curr_pos = pos;
                                            }
                                        }
                                        if (colored_contig != "") {
                                            if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                            else _range_discard++;
                                        }
                                    }
                                }
                                */
                            }
                            // end extend
                            std::set<int> positions_to_remove; // Positions (i.e. k-mers) along the current unitig that will be cut out
                            std::map<std::vector<int>, int> result_map; // Key = colors; Value = Position (i.e. k-mer)
                            std::stringstream ss; // For --combinations outputting aggregated colors
                            int int_to_print = 1;
                            if (perform_set_operations) {
                                for (const auto& expr : expressions) {
                                    ExpressionParser expr_parser(expr);
                                    Node* root = expr_parser.parse();
                                    std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                    for (const auto& k_elem : k_map) {
                                        int curr_pos = -1;
                                        auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                        std::string colored_contig = "";
                                        for (const auto& pos : k_elem.second) {
                                            if (set_operation_result.count(pos)) {
                                                std::string km = unitig.getUnitigKmer(pos).toString();
                                                if (curr_pos == -1) { colored_contig = km; }
                                                else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                else {
                                                    if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                    else _range_discard++;
                                                    colored_contig = km;
                                                }
                                                curr_pos = pos;
                                            }
                                        }
                                        if (colored_contig != "") {
                                            if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                            else _range_discard++;
                                        }
                                    }
                                }
                                continue; // colored contigs already extracted, continue to next iteration
                            }
                            if (!opt.distinguish_union) { // If we don't specify --union (since if --union is specified, we don't actually need to do anything, remove any k-mers/positions, etc.)
                                if (!opt.distinguish_all_but_one_color && !opt.distinguish_combinations) { // Workflow: Find k-mers unique/exclusive to each color
                                    std::string default_str = stringOpt(tmp_files.size(), "default"); // for n=8, generate "A\(AIBI...IH) B\(AIBI...IH) ... H\(AIBI...IH)"
                                    auto all_expr = split(default_str, ' ');
                                    for (const auto& expr : all_expr) {
                                        ExpressionParser expr_parser(expr);
                                        Node* root = expr_parser.parse();
                                        std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                        // write out what remains among the contigs
                                        for (const auto& k_elem : k_map) {
                                            int curr_pos = -1;
                                            auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                            std::string colored_contig = "";
                                            for (const auto& pos : k_elem.second) {
                                                if (set_operation_result.count(pos)) {
                                                    std::string km = unitig.getUnitigKmer(pos).toString();
                                                    if (curr_pos == -1) { colored_contig = km; }
                                                    else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                    else {
                                                        if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                        else _range_discard++;
                                                        colored_contig = km;
                                                    }
                                                    curr_pos = pos;
                                                }
                                            }
                                            if (colored_contig != "") {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                else _range_discard++;
                                            }
                                        }
                                    }          
                                    continue; // colored contigs already extracted, continue to next iteration
                                }
                                else if (!opt.distinguish_union && !opt.distinguish_combinations) { // workflow: opt.distinguish_all_but_one_color (e.g. if we have 8 colors, for each color, output all k-mers except those that are 8-colored)
                                    // method 1
                                    std::string all_but_one_str = stringOpt(tmp_files.size(), "all-but-one"); // generate "A\(AIBI...IH) B\(AIBI...IH) ... H\(AIBI...IH)"
                                    auto all_expr = split(all_but_one_str, ' ');
                                    for (const auto& expr : all_expr) {
                                        ExpressionParser expr_parser(expr);
                                        Node* root = expr_parser.parse();
                                        std::set<int> set_operation_result = computeSetOperation(root, k_map); // set of positions to keep
                                        // write out what remains among the contigs
                                        for (const auto& k_elem : k_map) {
                                            int curr_pos = -1;
                                            auto color = expr[0] - 'A'; // A=0, B=1, C=2, ..., H=7
                                            std::string colored_contig = "";
                                            for (const auto& pos : k_elem.second) {
                                                if (set_operation_result.count(pos)) {
                                                    std::string km = unitig.getUnitigKmer(pos).toString();
                                                    if (curr_pos == -1) { colored_contig = km; }
                                                    else if (pos == curr_pos + 1) { colored_contig += km[km.length() - 1]; }
                                                    else {
                                                        if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                        else _range_discard++;
                                                        colored_contig = km;
                                                    }
                                                    curr_pos = pos;
                                                }
                                            }
                                            if (colored_contig != "") {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re && k_elem.first == color) { oss << ">" << std::to_string(color) << "\n" << colored_contig << "\n"; _num_written++; }
                                                else _range_discard++;
                                            }
                                        }
                                    }
                                    continue; // colored contigs already extracted, continue to next iteration
                                }
                                else if (!opt.distinguish_all_but_one_color) { // workflow: opt.distinguish_combinations (output every combination [i.e. each area in the Venn diagram], with colors separated by underscores)
                                    std::set<int> exclusive_positions;
                                    for (const auto& k_elem : k_map) { // Iterate over each color and get exclusive positions for that color
                                        bool is_exclusive = true;
                                        std::set<int> current_positions = k_elem.second;
                                        for (const auto& other_elem : k_map) { // Check if the positions of this color appear in other colors
                                            if (k_elem.first != other_elem.first) {
                                                std::set<int> intersection;
                                                std::set_intersection(current_positions.begin(), current_positions.end(),
                                                    other_elem.second.begin(), other_elem.second.end(),
                                                    std::inserter(intersection, intersection.begin()));
                                                if (!intersection.empty()) {
                                                    is_exclusive = false;
                                                    break;
                                                }
                                            }
                                        }
                                        if (is_exclusive) { exclusive_positions.insert(current_positions.begin(), current_positions.end()); }
                                    }
                                    positions_to_remove.insert(exclusive_positions.begin(), exclusive_positions.end()); // Add the exclusive positions to the positions_to_remove
                                    if (k_map.size() >= 2) {
                                        std::vector<int> consolidated_key;
                                        for (auto it = k_map.begin(); it != k_map.end(); ++it) { // Iterate over k_map and build the consolidated key
                                            consolidated_key.push_back(it->first);
                                        }
                                        int integer_value = *(k_map.begin()->second.begin());
                                        result_map[consolidated_key] = integer_value; // {0 1 2} : 0 (color : position)                                        
                                    }
                                    for (const auto& k_elem : result_map) {
                                        int curr_pos = -1;
                                        std::string colored_contig = "";
                                        auto color = k_elem.first;
                                        auto pos = k_elem.second;
                                        for (int i = 0; i < color.size(); ++i) {
                                            ss << color[i];
                                            if (i < color.size() - 1) {
                                                ss << "_";
                                            }
                                        }
                                    }
                                }
                            } // end if !opt.distinguish_union
                            
                            // Now, write out what remains among the contigs
                            for (const auto& k_elem : k_map) {
                                int curr_pos = -1;
                                std::string colored_contig = "";
                                auto color = k_elem.first;
                                std::string color_key = std::to_string(color);
                                //std::string contig_metadata = " :" + unitig.dist + "," + unitig.len + "," + unitig.size + "," + unitig.strand;
                                if (opt.distinguish_combinations) { color_key = ss.str(); int_to_print++; }
                                if (int_to_print == k_map.size()) {
                                    for (const auto& pos : k_elem.second) {
                                        if (!positions_to_remove.count(pos)) {
                                            std::string km = unitig.getUnitigKmer(pos).toString();
                                            if (curr_pos == -1) { // How to correspond color?
                                                colored_contig = km;
                                            }
                                            else if (pos == curr_pos + 1) {
                                                colored_contig += km[km.length() - 1];
                                            }
                                            else {
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re) { oss << ">" << color_key << "\n" << colored_contig << "\n"; _num_written++; }
                                                else _range_discard++;
                                                colored_contig = km;
                                            }
                                            curr_pos = pos;
                                        }
                                    }
                                    if (colored_contig != "") {
                                        if (colored_contig.length() >= rb && colored_contig.length() <= re) { oss << ">" << color_key << "\n" /*<< sequence_to_prepend*/ << colored_contig /*<< sequence_to_append*/ << "\n"; _num_written++; }
                                        else _range_discard++;
                                    }
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
