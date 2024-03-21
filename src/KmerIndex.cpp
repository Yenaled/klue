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
#include <cctype>

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

int my_mkdir_kmer_index(const char* path, mode_t mode) {
#ifdef _WIN64
    return mkdir(path);
#else
    return mkdir(path, mode);
#endif
}

std::string generate_tmp_file(std::string seed, std::string tmp_dir) {
    struct stat stFileInfo;
    auto intStat = stat(tmp_dir.c_str(), &stFileInfo);
    if (intStat == 0) {
        // file/dir exits
        if (!S_ISDIR(stFileInfo.st_mode)) {
            cerr << "Error: file " << tmp_dir << " exists and is not a directory" << endl;
            exit(1);
        }
    }
    else {
        // create directory
        if (my_mkdir_kmer_index(tmp_dir.c_str(), 0777) == -1) {
            cerr << "Error: could not create directory " << tmp_dir << endl;
            exit(1);
        }
    }
    std::string base = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::string tmp_file = "klue.";
    srand((unsigned int)std::hash<std::string>{}(seed));
    int pos;
    while (tmp_file.length() < 32) {
        pos = ((rand() % (base.size() - 1)));
        tmp_file += base.substr(pos, 1);
    }
    return tmp_dir + "/" + tmp_file;
}

// Perform DFS traversal starting from unitig tail
/*
std::string extendUnitig(const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current,
    std::unordered_set<std::string>& visited,
    const int& color,
    std::string current_str,
    const std::unordered_set<int>& superset_colors,
    const int& k) {

    std::string result;
    bool hasValidSuccessor = false;

    std::string colored_unitig = current.getUnitigHead().toString() + std::to_string(color);
    if (visited.find(colored_unitig) != visited.end() || current.isEmpty) { return result; }
    visited.insert(colored_unitig);

    auto successors = current.getSuccessors();
    if (successors.begin() == successors.end()) { return result; } // If there are no successors, mark as reached end

    for (const auto& next : successors) {
        if (visited.find(next.getUnitigHead().toString() + std::to_string(color)) == visited.end()) {
            UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
            UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
            std::unordered_set<int> colors;
            for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
            if (colors.find(color) != colors.end()) { // Check that next unitig color matches current traversal color
                hasValidSuccessor = true;
                std::string str;
                if (next.strand) {
                    for (int i = next.dist; i < next.len; ++i) {
                        str += next.getUnitigKmer(i).toString().substr(k - 1);
                    }
                }
                else {
                    for (int i = next.dist; i < next.len; ++i) {
                        str = next.getUnitigKmer(i).twin().toString().substr(k - 1) + str;
                    }
                }

                std::string extend_str = extendUnitig(next, visited, color, str, superset_colors, k);
                current_str += extend_str;
            }
        }
    }
    //if (hasValidSuccessor) { result += current_str; } // at leaf node
    result += current_str;
    return result;
}
*/

struct Extend {
    std::string result;
    int color;
};

Extend extendUnitig(const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current,
    std::unordered_set<std::string>& visited,
    const int& color,
    std::string current_str,
    const std::unordered_set<int>& superset_colors,
    const int& k) {

    Extend traversal_result;
    bool hasValidSuccessor = false;

    std::string colored_unitig = current.getUnitigHead().toString() + std::to_string(color);
    if (visited.find(colored_unitig) != visited.end() || current.isEmpty) {
        traversal_result.color = color;
        return traversal_result;
    }
    visited.insert(colored_unitig);

    auto successors = current.getSuccessors();
    if (successors.begin() == successors.end()) {
        traversal_result.color = color;
        return traversal_result;
    }

    for (const auto& next : successors) {
        if (visited.find(next.getUnitigHead().toString() + std::to_string(color)) == visited.end()) {
            UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
            UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
            std::unordered_set<int> colors;
            for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
            if (colors.find(color) != colors.end()) {
                hasValidSuccessor = true;
                std::string str;
                if (next.strand) {
                    for (int i = next.dist; i < next.len; ++i) {
                        str += next.getUnitigKmer(i).toString().substr(k - 1);
                    }
                }
                else {
                    for (int i = next.dist; i < next.len; ++i) {
                        str = next.getUnitigKmer(i).twin().toString().substr(k - 1) + str;
                    }
                }
                Extend extendResult = extendUnitig(next, visited, color, str, superset_colors, k);
                current_str += extendResult.result;
            }
        }
    }
    traversal_result.result = current_str;
    traversal_result.color = color;
    return traversal_result;
}

std::string generate_rev_comp(const std::string& sequence) {
    std::string rc;
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (*it) {
        case 'A': rc += 'T'; break;
        case 'T': rc += 'A'; break;
        case 'G': rc += 'C'; break;
        case 'C': rc += 'G'; break;
        }
    }
    return rc;
}

struct BubblePath {
    std::vector<std::string> variations;
    bool terminalNode = false;
    bool switchKmers = false;
};

BubblePath bubbleVariations(ColoredCDBG<void>& ccdbg,
    const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current,
    const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& terminal_kmer,
    const std::unordered_set<int>& superset_colors,
    int color,
    std::unordered_set<std::string>& visited,
    int k,
    std::string tag = "_0") {

    BubblePath pathResult;

    if (current == terminal_kmer) { // terminate traversal when we reach terminal node
        pathResult.terminalNode = true;
        return pathResult;
    }

    UnitigMap<DataAccessor<void>, DataStorage<void>, false> current_copy = current;
    UnitigMap<DataAccessor<void>, DataStorage<void>, false> terminal_kmer_copy = terminal_kmer;

    // If there are no successors, switch terminal_kmer and current to find variations
    if (current.getSuccessors().begin() == current.getSuccessors().end()) {
        current_copy = terminal_kmer;
        terminal_kmer_copy = current;
        pathResult.switchKmers = true;
    }

    // Traverse successors to find variations
    auto successors = current_copy.getSuccessors();
    int num_successors = 0;
    // int num_successors = std::distance(successors.begin(), successors.end());
    for (const auto& next : successors) {
        UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
        UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
        std::unordered_set<int> colors;
        for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
        if (colors.find(color) != colors.end()) { // check that next unitig color matches current traversal color
            num_successors++;
        }
    }
    int divergence_counter = 1;
    for (const auto& next : successors) {
        //if (visited.find(next.getUnitigHead().toString() + std::to_string(color)) == visited.end()) {
        UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
        UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
        std::unordered_set<int> colors;
        for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
        if (colors.find(color) != colors.end()) { // check that next unitig color matches current traversal color
            std::string str;
            if (next.strand) {
                // DEBUG
                //std::cout << "next: " << next.strand << " " << next.getUnitigHead().toString() << " -> " << next.getUnitigTail().toString() << "\n";
                for (int i = next.dist; i < next.len; ++i) {
                    str += next.getUnitigKmer(i).toString().substr(k - 1);
                }
            }
            else {
                // DEBUG
                //std::cout << "twin: " << !next.strand << " " << next.getUnitigTail().twin().toString() << " -> " << next.getUnitigHead().twin().toString() << "\n";
                for (int i = next.dist; i < next.len; ++i) {
                    str = next.getUnitigKmer(i).twin().toString().substr(k - 1) + str;
                }
            }
            // keep track of visited nodes (not needed for now)
            /*
            std::string colored_unitig = current_copy.getUnitigHead().toString() + std::to_string(color);
            if (visited.find(colored_unitig) != visited.end() || current.isEmpty) { return pathResult; }
            visited.insert(colored_unitig);
            */
            pathResult.variations.push_back(str);
            std::string temp_str;
            for (auto& p : pathResult.variations) {
                if (std::any_of(p.begin(), p.end(), [](unsigned char c) { return std::isdigit(c); })) { // if there are already strand tags, add to new string
                    temp_str = std::to_string(current_copy.strand) + p.substr(1, p.length() - 2) + std::to_string(next.strand);
                }
                else {
                    p = std::to_string(current_copy.strand) + p + std::to_string(next.strand);
                }
            }
            if (!temp_str.empty()) { pathResult.variations.push_back(temp_str); }
            temp_str.clear();
            // add tag to track diverging paths
            std::string divergence_tag = tag;
            if (num_successors > 1) {
                divergence_tag = tag.substr(0, tag.size() - 1) + std::to_string(divergence_counter++);
            }
            for (auto& p : pathResult.variations) {
                int tagged = p.find_last_of('_');
                if (tagged != std::string::npos) {
                    temp_str = p.substr(0, tagged) + divergence_tag;
                }
                else {
                    p += divergence_tag;
                }
            }
            BubblePath tempResult = bubbleVariations(ccdbg, next, terminal_kmer_copy, superset_colors, color, visited, k, divergence_tag);
            pathResult.variations.insert(pathResult.variations.end(), tempResult.variations.begin(), tempResult.variations.end());
            if (tempResult.terminalNode) {
                pathResult.terminalNode = true;
                break;
            }
        }
        //}
    }
    return pathResult;
}

struct Bubble {
    std::string bubble_left;
    std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::string>>> variations; // maybe try another data structure?
    std::string bubble_right;
};

Bubble exploreBubble(ColoredCDBG<void>& ccdbg,
    const UnitigMap<DataAccessor<void>, DataStorage<void>, false>& current,
    std::unordered_set<std::string>& visited,
    const std::unordered_set<int>& superset_colors,
    int color,
    int k) {

    auto successors = current.getSuccessors();
    auto predecessors = current.getPredecessors();
    if (successors.begin() == successors.end() && predecessors.begin() == predecessors.end()) { return {}; } // Terminal node

    std::string colored_unitig = current.getUnitigHead().toString() + std::to_string(color);
    if (visited.find(colored_unitig) != visited.end()) { return {}; } // Dont't revist nodes with the same color
    visited.insert(colored_unitig);

    Bubble bubblePath;

    // Explore predecessors to find bubble_left
    for (const auto& prev : predecessors) {
        UnitigColors::const_iterator it_prev = prev.getData()->getUnitigColors(prev)->begin(prev);
        UnitigColors::const_iterator it_prev_end = prev.getData()->getUnitigColors(prev)->end();
        std::unordered_set<int> colors_prev;
        for (; it_prev != it_prev_end; ++it_prev) { colors_prev.insert(it_prev.getColorID()); }
        if (colors_prev.find(color) != colors_prev.end() &&                                               // If color is found in set of valid predecessor colors
            (colors_prev.size() > 1 || prev.getPredecessors().begin() == prev.getPredecessors().end())) { // AND if there are multiple colors or no predecessors
            bubblePath.bubble_left = prev.getUnitigHead().toString(); // LATER adjust substr to get Head -> Tail
        }
    }

    // Explore successors to find bubble_right and variation
    for (const auto& next : successors) {
        UnitigColors::const_iterator it_next = next.getData()->getUnitigColors(next)->begin(next);
        UnitigColors::const_iterator it_next_end = next.getData()->getUnitigColors(next)->end();
        std::unordered_set<int> colors;
        for (; it_next != it_next_end; ++it_next) { colors.insert(it_next.getColorID()); }
        if (colors.find(color) != colors.end()) {
            if (colors.size() > 1 || next.getSuccessors().begin() == next.getSuccessors().end()) { // If multiple colors present or no successors
                bubblePath.bubble_right = next.getUnitigHead().toString(); // LATER: adjust substr to get Head -> Tail                
            }
            else {
                Bubble tempResult = exploreBubble(ccdbg, next, visited, superset_colors, color, k);
                if (!tempResult.bubble_right.empty() && bubblePath.bubble_right.empty()) { bubblePath.bubble_right = tempResult.bubble_right; }
                if (!tempResult.bubble_left.empty() && bubblePath.bubble_left.empty()) { bubblePath.bubble_left = tempResult.bubble_left; }
                std::unordered_map<std::string, std::string> bubble_map;
                // if bubble_left and bubble_right are not empty and not the same, they are a valid flanking pair
                if (!bubblePath.bubble_left.empty() && !bubblePath.bubble_right.empty() && bubblePath.bubble_left != bubblePath.bubble_right) {
                    Kmer bubble_left_kmer(bubblePath.bubble_left.c_str());
                    Kmer bubble_right_kmer(bubblePath.bubble_right.c_str());
                    std::unordered_set<std::string> variation_visited;
                    UnitigMap<DataAccessor<void>, DataStorage<void>, false> um_left = ccdbg.find(bubble_left_kmer, false);   // start node (left -> right)
                    UnitigMap<DataAccessor<void>, DataStorage<void>, false> um_right = ccdbg.find(bubble_right_kmer, false); // terminal node
                    // Found valid start and end unitigs, now find variation 
                    if (!um_left.isEmpty && !um_right.isEmpty) {
                        //std::cout << "found valid path at: " << bubblePath.bubble_left << " -> " << bubblePath.bubble_right << "\n";
                        for (const auto& c : superset_colors) {
                            BubblePath bubble_path = bubbleVariations(ccdbg, um_left, um_right, superset_colors, c, variation_visited, k);
                            if (bubble_path.switchKmers) {
                                auto temp = bubblePath.bubble_left;
                                bubblePath.bubble_left = bubblePath.bubble_right;
                                bubblePath.bubble_right = temp;

                                UnitigMap<DataAccessor<void>, DataStorage<void>, false> um_temp = um_right;
                                um_right = um_left;
                                um_left = um_temp;
                            }
                            int i = 0;
                            for (const auto& var : bubble_path.variations) {
                                int tag = var.find_last_of('_');
                                if (tag != std::string::npos) {
                                    bubblePath.variations[c][std::stoi(var.substr(tag + 1))][i] = var.substr(0, var.size() - 2); // remove "_X" tag
                                    i++;
                                }
                            }
                            bool concat = true; // flag to keep track of whether any concatenations were made in the last iteration
                            while (concat) {
                                // method 1 for concatenation (ok)
                                for (int i = 0; i < bubblePath.variations[c].size(); ++i) { // for each path
                                    concat = false; // reset flag at the beginning of each iteration
                                    for (int j = 0; j < bubblePath.variations[c][i].size(); ++j) { // for each variation within path
                                        if (bubblePath.variations[c][i][j].empty()) { continue; } // already concatenated
                                        for (int k = j + 1; k < bubblePath.variations[c][i].size(); ++k) { // for each subsequent variation
                                            if (bubblePath.variations[c][i][k].empty()) { continue; } // already concatenated
                                            if (bubblePath.variations[c][i][j].back() == bubblePath.variations[c][i][k].front()) {
                                                bubblePath.variations[c][i][j] += bubblePath.variations[c][i][k];
                                                //std::cout << "bubblePath.variations[" << c << "][" << i << "][" << j << "]: " << bubblePath.variations[c][i][j] << "\n";
                                                bubblePath.variations[c][i][k].clear();
                                                concat = true;
                                                break; // restart scanning from the beginning
                                            }
                                        }
                                    }
                                }
                                // method 2 for concatenation (incomplete)
                                /*
                                for (const auto& i_pair : bubblePath.variations[c]) {
                                    int i = i_pair.first;
                                    std::vector<int> sorted_keys;  // sort the keys to ensure a specific order (want to concatenate successive variations in a specific order)
                                    for (const auto& j_pair : i_pair.second) {
                                        sorted_keys.push_back(j_pair.first);
                                    }
                                    std::sort(sorted_keys.begin(), sorted_keys.end(), std::greater<int>());
                                    std::string* prev_ptr = nullptr;
                                    std::vector<int> keys_to_remove;
                                    for (int j : sorted_keys) {
                                        std::string& var_string = bubblePath.variations[c][i][j];
                                        std::cout << "j: " << j << " " << var_string << "\n";
                                        if (prev_ptr != nullptr && !prev_ptr->empty() && !var_string.empty() && prev_ptr->back() == var_string.front()) {
                                            *prev_ptr = var_string + *prev_ptr;
                                            std::cout << "prev_ptr: " << *prev_ptr << "\n";
                                            keys_to_remove.push_back(j);
                                            concat = true;
                                        }
                                        else {
                                            prev_ptr = &var_string; // update prev_ptr only if no concatenation
                                        }
                                    }
                                    for (int key : keys_to_remove) {
                                        bubblePath.variations[c][i].erase(key);
                                    }
                                }
                                */
                            }

                            // print variations using key-value iteration
                            /*
                            for (const auto& pair : bubblePath.variations[c]) { // for each {path : variation}
                                for (const auto& var : pair.second) { // for each variation
                                    if (!var.second.empty()) {
                                        //var.second.erase(remove_if(var.second.begin(), var.second.end(), ::isdigit), var.second.end()); // remove strand tags
                                        std::cout << "bubblePath.variations[" << c << "][" << i << "][" << j << "]: " << value << "\n";
                                    }
                                }
                            }
                            */
                            // print variations using range-based for loop
                            /*
                            for (int i = 0; i < bubblePath.variations[c].size(); ++i) { // for each path
                                for (int j = 0; j < bubblePath.variations[c][i].size(); ++j) { // for each variation within path
                                    std::cout << "bubblePath.variations[" << c << "][" << i << "][" << j << "]: " << bubblePath.variations[c][i][j] << "\n";
                                }
                            }
                            */

                            // prepend common prefix to all divergent paths
                            int cth = 0; int ith = 0; int jth = 0; // change names
                            bool prepend = false;
                            if (bubblePath.variations[c].find(0) != bubblePath.variations[c].end()) { // check if i=0 exists
                                // Iterate over each entry in i=0 to prepend it to successive i's
                                for (const auto& var : bubblePath.variations[c][0]) { // for each variation in path i=0
                                    std::string to_prepend = bubblePath.variations[c][0][var.first];
                                    for (auto& pair : bubblePath.variations[c]) { // for each {path : variation} in path > 1
                                        int i = pair.first;
                                        if (i == 0) { continue; }
                                        for (auto& j_pair : bubblePath.variations[c][i]) { // prepend to each variation in path i>0
                                            int j = j_pair.first;
                                            if (!to_prepend.empty() && !bubblePath.variations[c][i][j].empty()) {
                                                // METHOD 1 naive prepend
                                                /*
                                                bubblePath.variations[c][i][j] = to_prepend + bubblePath.variations[c][i][j];
                                                cth = c; ith = 0; jth = var.first;
                                                prepend = true;
                                                */
                                                // METHOD 2 check for strand match prepend
                                                // Take the first character (digit) of bubblePath.variations[c][i][j] and compare it to the last character (digit) of to_prepend;
                                                // If there is no match, remove the last three characters of to_pretend and comapare again.
                                                // Repeat until a match is found or to_prepend is empty
                                                while (!to_prepend.empty() && to_prepend.back() != bubblePath.variations[c][i][j].front()) {
                                                    // If no match, remove the last 3 characters or set to empty if not possible
                                                    // currently assuming that snp is cause for mismatch
                                                    if (to_prepend.size() > 3) {
                                                        to_prepend.erase(to_prepend.size() - 3);
													}
                                                    else {
                                                        to_prepend.clear();
														break; // Exit loop if to_prepend is or becomes empty
													}
                                                }
                                                if (!to_prepend.empty()) {
                                                    bubblePath.variations[c][i][j] = to_prepend + bubblePath.variations[c][i][j];
                                                    cth = c; ith = 0; jth = var.first;
                                                    prepend = true;
                                                }
                                            }
                                        }
                                    }
                                }                                
                            }
                            // Remove common prefix from divergent paths
                            if (prepend) { bubblePath.variations[cth][ith][jth].erase(); }
                            // Remove reverse complement duplicates
                            for (int i = 0; i < bubblePath.variations[c].size(); ++i) {
                                for (int j = 0; j < bubblePath.variations[c][i].size(); ++j) {
                                    bubblePath.variations[c][i][j].erase(remove_if(bubblePath.variations[c][i][j].begin(), bubblePath.variations[c][i][j].end(), ::isdigit), bubblePath.variations[c][i][j].end()); // remove strand tags
                                    // compare all but the last (k-1) characters and see if they are reverse complement
                                    for (int p = j + 1; p < bubblePath.variations[c][i].size(); ++p) {
                                        bubblePath.variations[c][i][p].erase(remove_if(bubblePath.variations[c][i][p].begin(), bubblePath.variations[c][i][p].end(), ::isdigit), bubblePath.variations[c][i][p].end()); // remove strand tags
                                        if (!bubblePath.variations[c][i][p].empty() && !bubblePath.variations[c][i][j].empty()) {
                                            if (generate_rev_comp(bubblePath.variations[c][i][p].substr(0, bubblePath.variations[c][i][p].length() - k)) == bubblePath.variations[c][i][j].substr(0, bubblePath.variations[c][i][j].length() - k)) {
                                                bubblePath.variations[c][i][j].erase(); // remove if duplicate (reverse complement)
                                            }
                                        }
                                    }
                                }
                            }
                            // Check that bubble_left and bubble_right are valid successor/predecessors of the variation
                            /*                            for (int i = 0; i < bubblePath.variations[c].size(); ++i) {
                                for (int j = 0; j < bubblePath.variations[c][i].size(); ++j) {
                                    bubblePath.variations[c][i][j] = bubblePath.variations[c][i][j].substr(0, bubblePath.variations[c][i][j].length() - 1);
                                    // currently finds SNPs or longer variations -- NO deletions
                                    // for deletions, do recursive search (try A,T,G,C appended to (k-1)mer and check if it exists in the graph)
                                    if (bubblePath.variations[c][i][j].length() > k) {
                                        
                                        // last (k-1) str.substr(str.length() - (k - 1));
                                        // first (k-1) str.substr(0, k-1);
                                        // case 1
                                        // bubble_left + variation + bubble_right (variation overlaps (k-1)mer with bubble_right), i.e. we want to prepend bubble_left
                                        if (bubblePath.variations[c][i][j].substr(bubblePath.variations[c][i][j].length() - (k - 1)) == bubblePath.bubble_right.substr(0, k - 1)) {
                                            Kmer km = Kmer(bubblePath.variations[c][i][j].substr(bubblePath.variations[c][i][j].length() - k).c_str()); // get last kmer
                                            UnitigMap<DataAccessor<void>, DataStorage<void>, false> um = ccdbg.find(km, false);
                                            if (!um.isEmpty) {
                                                for (const auto& next : um.getSuccessors()) {
                                                    std::cout << "bubblePath.variations[" << c << "][" << i << "][" << j << "]: " << bubblePath.variations[c][i][j] << "\n";
                                                    if(next.strand){ std::cout << "next: " << next.strand << " " << next.getUnitigHead().toString() << " -> " << next.getUnitigTail().toString() << "\n"; }
                                                    else { std::cout << "twin: " << !next.strand << " " << next.getUnitigTail().twin().toString() << " -> " << next.getUnitigHead().twin().toString() << "\n"; }
												}
                                            }
                                        }
                                        // case 2
                                        // bubble_left + variation + bubble_right (variation overlaps with bubble_left) i.e. we want to prepend bubble_right 
                                        if (bubblePath.variations[c][i][j].substr(0, k - 1) == bubblePath.bubble_left.substr(bubblePath.bubble_left.length() - (k - 1))){
                                            Kmer km = Kmer(bubblePath.variations[c][i][j].substr(0, k).c_str()); // get first kmer
                                            UnitigMap<DataAccessor<void>, DataStorage<void>, false> um = ccdbg.find(km, false);
                                            if (!um.isEmpty) {
                                                for (const auto& prev : um.getPredecessors()) {
                                                    std::cout << "bubblePath.variations[" << c << "][" << i << "][" << j << "]: " << bubblePath.variations[c][i][j] << "\n";
                                                    if (prev.strand) { std::cout << "prev: " << prev.strand << " " << prev.getUnitigHead().toString() << " -> " << prev.getUnitigTail().toString() << "\n"; }
                                                    else { std::cout << "twin: " << !prev.strand << " " << prev.getUnitigTail().twin().toString() << " -> " << prev.getUnitigHead().twin().toString() << "\n"; }
                                                }
                                            }
										}
                                    }
                                }
                            }
                            */

                            // Return result
                            for (int i = 0; i < bubblePath.variations[c].size(); ++i) {
                                for (int j = 0; j < bubblePath.variations[c][i].size(); ++j) {
                                    bubblePath.variations[c][i][j] = bubblePath.variations[c][i][j].substr(0, bubblePath.variations[c][i][j].length() - 1);
                                    // currently finds SNPs or longer variations -- NO deletions
                                    // for deletions, do recursive search (try A,T,G,C appended to (k-1)mer and check if it exists in the graph)
                                    if (bubblePath.variations[c][i][j].length() > k) { // should always be > k, but to prevent substr out of bounds in case -- NOT in the case of deletion
                                        bubblePath.variations[c][i][j].erase(remove_if(bubblePath.variations[c][i][j].begin(), bubblePath.variations[c][i][j].end(), ::isdigit), bubblePath.variations[c][i][j].end()); // remove strand tags
                                        //std::cout << "bubblePath.variations[" << c << "][" << i << "][" << j << "]: " << bubblePath.variations[c][i][j] << "\n";
                                        Kmer km = Kmer(bubblePath.variations[c][i][j].substr(0, k).c_str()); // get first kmer
                                        UnitigMap<DataAccessor<void>, DataStorage<void>, false> um = ccdbg.find(km, false);
                                        if (!um.isEmpty) {
                                            if (bubblePath.variations[c][i][j].substr(bubblePath.variations[c][i][j].length() - (k - 1)) == um_left.getUnitigTail().twin().toString().substr(1)) {
                                                // this means that um_left is successor, we want to prepend um_right to the variation
                                                if (um_right.strand) {
                                                    bubblePath.variations[c][i][j] = um_right.getUnitigTail().toString().substr(1) + bubblePath.variations[c][i][j];
                                                }
                                                else {
                                                    bubblePath.variations[c][i][j] = um_right.getUnitigHead().toString().substr(0, k - 1) + bubblePath.variations[c][i][j];
                                                }
                                            }
                                            else {
                                                if (um_left.strand) {
                                                    bubblePath.variations[c][i][j] = um_left.getUnitigTail().toString().substr(1) + bubblePath.variations[c][i][j];
                                                }
                                                else {
                                                    bubblePath.variations[c][i][j] = um_left.getUnitigHead().toString().substr(0, k - 1) + bubblePath.variations[c][i][j];
                                                }
                                            }
                                        }
                                        // fix
                                        else { bubblePath.variations[c][i][j].erase(); }
                                    }
                                    
                                }
                            }
                        }
                        return bubblePath;
                    }
                }
            }
        }
    }
    return bubblePath;
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

// switch to NAND boolean logic


/*

"not" : NOT(A) = A NAND A
"and" : A AND B = ( A NAND B ) NAND ( A NAND B )
"or"  : A OR B = ( A NAND A ) NAND ( B NAND B )
"nor" : A NOR B = [ ( A NAND A ) NAND ( B NAND B ) ] NAND [ ( A NAND A ) NAND ( B NAND B ) ]
"xor" : A XOR B = [ A NAND ( A NAND B ) ] NAND [ B NAND ( A NAND B ) ]
"xnor": A XNOR B = [ ( A NAND A ) NAND ( B NAND B ) ] NAND ( A NAND B )

*/
/*

// NAND operation
bool nand(bool a, bool b) {
    return !(a && b);
}

// Implement other operations using NAND
bool NOT(bool a) {
    return nand(a, a);
}

bool AND(bool a, bool b) {
    return NOT(nand(a, b));
}

bool OR(bool a, bool b) {
    return nand(NOT(a), NOT(b));
}

bool XOR(bool a, bool b) {
    return nand(nand(a, nand(a, b)), nand(b, nand(a, b)));
}

bool NOR(bool a, bool b) {
    return NOT(OR(a, b));
}

bool XNOR(bool a, bool b) {
    return NOT(XOR(a, b));
}
*/


// Recursively compute set operations (union and intersection are nodes) that correspond to specific sets (leaves)
std::set<int> computeSetOperation(const Node* root, const std::map<int, std::set<int>>& k_map) {
    if (!root) { return {}; }
    if (root->value != 'U' && root->value != 'I' && root->value != '\\' && root->value != 'N' && root->value != 'X') {
        auto it = k_map.find(root->value - 'A'); // 'A' maps to 0, 'B' to 1, etc.
        return (it != k_map.end()) ? it->second : std::set<int>{};
    }
    // Construct the universal set from all sets in k_map
    std::set<int> universalSet;
    for (const auto& pair : k_map) {
        universalSet.insert(pair.second.begin(), pair.second.end());
    }
    auto leftSet = computeSetOperation(root->left, k_map);
    auto rightSet = computeSetOperation(root->right, k_map);
    std::set<int> result;
    if (root->value == 'U') { // (Logical AND) union of sets A,B. User-inputted as "AUB"
        std::set_union(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == 'I') { // (Logical OR) intersection of sets A,B. User-inputted as "AIB"
        std::set_intersection(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == '\\') { // Difference of sets A,B. User-inputted as "A\B"
        std::set_difference(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(result, result.end()));
    }
    else if (root->value == 'N') { // Not both of sets A,B. User-inputted as "ANB"
        std::set<int> intersection;
        std::set_intersection(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(intersection, intersection.end()));
        std::set_difference(universalSet.begin(), universalSet.end(), intersection.begin(), intersection.end(), std::inserter(result, result.end()));
    }
    // DEBUG: add logical NOR operator -- the dual of NAND
    else if (root->value == 'X') { // Exclusive OR of sets A,B. User-inputted as "AXB"
        std::set<int> tempUnion, tempIntersection, xorResult;
        std::set_union(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(tempUnion, tempUnion.end()));
        std::set_intersection(leftSet.begin(), leftSet.end(), rightSet.begin(), rightSet.end(), std::inserter(tempIntersection, tempIntersection.end()));
        std::set_difference(tempUnion.begin(), tempUnion.end(), tempIntersection.begin(), tempIntersection.end(), std::inserter(xorResult, xorResult.end()));
        return xorResult;
    }
    return result;
}

std::string stringOpt(size_t n, std::string option) {
    std::string result;
    // Skip letters 'I', 'N', 'X', 'U' because they are used in set operations
    auto skipLetters = [&](char& ch) {
        if (ch >= 'I') { ch += 1; }
        if (ch >= 'N') { ch += 1; }
        if (ch >= 'X') { ch += 1; }
        if (ch >= 'U') { ch += 1; }
        };
    if (option == "default") {
        for (size_t i = 0; i < n; ++i) {
            char currentChar = 'A' + i;
            skipLetters(currentChar);
            result += currentChar;
            result += "\\(";
            for (size_t j = 0; j < n; ++j) {
                char innerChar = 'A' + j;
                skipLetters(innerChar);
                if (i != j) {
                    result += innerChar;
                    if (j < n - 1) {
                        if (!(i == n - 1 && j == n - 2)) {
                            result += 'U';
                        }
                    }
                }
            }
            result += ')';
            if (i != n - 1) { result += ' '; }
        }
    }
    else if (option == "all-but-one") {
        for (size_t i = 0; i < n; ++i) {
            char currentChar = 'A' + i;
            skipLetters(currentChar);
            result += currentChar;
            result += "\\(";
            for (size_t j = 0; j < n; ++j) {
                char innerChar = 'A' + j;
                skipLetters(innerChar);
                result += innerChar;
                if (j < n - 1) { result += 'I'; }
            }
            result += ')';
            if (i != n - 1) { result += ' '; }
        }
    }
    return result;
}
// End set operations

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
                        tmp_files[i] = generate_tmp_file(opt.distinguish_output_fasta + fasta + std::to_string(i), opt.tmp_dir);
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
            tmp_files.push_back(generate_tmp_file(opt.distinguish_output_fasta + fasta, opt.tmp_dir));
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
    // for bubble
    std::mutex mutex_bubbles;
    std::ofstream const_left_file("constant_left.fa", std::ofstream::out); // change to user-defined output file
    std::ofstream variation_file("variation.fa", std::ofstream::out);
    std::ofstream const_right_file("constant_right.fa", std::ofstream::out);
    // file validation
    if (!const_left_file.is_open() || !variation_file.is_open() || !const_right_file.is_open()) {
        std::cerr << "[WARNING]: Error opening output files." << std::endl;
        return; // exit
    }
    // continue with default processing
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
                        std::ostringstream const_left_stream, variation_stream, const_right_stream; // for bubble                        
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
                            std::string to_append = "";
                            // begin extend 
                            if (opt.extend) {
                                int color = *superset_colors.begin();
                                std::unordered_set<std::string> visited;
                                Extend traversal = extendUnitig(unitig, visited, color, "", superset_colors, k);
                                to_append = traversal.result;

                            }
                            //std::cout << "to_append: " << to_append << "\n";
                            // end extend
                            // begin bubble
                            if (opt.bubble) {
                                int color = *superset_colors.begin();
                                std::unordered_set<std::string> visited;
                                Bubble result = exploreBubble(ccdbg, unitig, visited, superset_colors, color, k);
                                for (const auto& c : superset_colors) {
                                    for (int i = 0; i < result.variations[c].size(); ++i) {
                                        for (int j = 0; j < result.variations[c][i].size(); ++j) {
                                            if (!result.bubble_left.empty() && !result.bubble_right.empty() && result.variations[c][i][j].length() > k - 1) {
                                                const_left_stream << ">" << c << "\n" << result.bubble_left << "\n";
                                                variation_stream << ">" << c << "\n" << result.variations[c][i][j] << "\n";
                                                const_right_stream << ">" << c << "\n" << result.bubble_right << "\n";
                                            }
                                        }
                                    }
                                }
                            }
                            // end bubble
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
                                        std::set<int> current_positions = k_elem.second;
                                        if (current_positions.size() == 1) { // Check if the k-mer is exclusive to one color
                                            // This k-mer is unique to one color, so process it accordingly
                                            exclusive_positions.insert(current_positions.begin(), current_positions.end());
                                            std::vector<int> single_color_key = { k_elem.first };
                                            int integer_value = *(current_positions.begin());
                                            result_map[single_color_key] = integer_value; // Add this single-color k-mer to the result map
                                        }
                                        else {
                                            bool is_exclusive = true;
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
                                    }
                                    positions_to_remove.insert(exclusive_positions.begin(), exclusive_positions.end()); // Add the exclusive positions to the positions_to_remove
                                    // updated color key method below
                                   /*
                                    if (k_map.size() >= 2) {
                                        std::vector<int> consolidated_key;
                                        for (auto it = k_map.begin(); it != k_map.end(); ++it) { // Iterate over k_map and build the consolidated key
                                            consolidated_key.push_back(it->first);
                                        }
                                        int integer_value = *(k_map.begin()->second.begin());
                                        result_map[consolidated_key] = integer_value; // {0 1 2} : 0 (color : position)
                                    }
                                    */
                                    std::vector<int> consolidated_key;
                                    for (auto it = k_map.begin(); it != k_map.end(); ++it) {
                                        consolidated_key.push_back(it->first);
                                    }
                                    int integer_value = *(k_map.begin()->second.begin());
                                    result_map[consolidated_key] = integer_value; // {0} or {0 1 2} : 0 (color : position)

                                    for (const auto& k_elem : result_map) {
                                        int curr_pos = -1;
                                        std::string colored_contig;
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
                                                if (colored_contig.length() >= rb && colored_contig.length() <= re) {
                                                    oss << ">" << color_key << "\n" << colored_contig + to_append << "\n"; _num_written++;
                                                }
                                                else _range_discard++;
                                                colored_contig = km;
                                            }
                                            curr_pos = pos;
                                        }
                                    }
                                    if (colored_contig != "") {
                                        if (colored_contig.length() >= rb && colored_contig.length() <= re) {
                                            // if opt.extend, then we need to append the traversal.result to the colored_contig
                                            oss << ">" << color_key << "\n" << colored_contig + to_append << "\n"; _num_written++;
                                        }
                                        else _range_discard++;
                                    }
                                }
                            }
                        }
                        {
                            std::unique_lock<std::mutex> lock(mutex_unitigs);
                            o << oss.str();
                        }
                        // write bubble results to respective files
                        {
                            std::unique_lock<std::mutex> lock(mutex_bubbles);
                            const_left_file << const_left_stream.str();
                            variation_file << variation_stream.str();
                            const_right_file << const_right_stream.str();
                        }
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
    const_left_file.flush(); variation_file.flush(); const_right_file.flush(); // for bubble
    // ADD close files

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