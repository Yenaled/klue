#ifndef KURE_STRINGSEARCH_H
#define KURE_STRINGSEARCH_H

#include "common.h"
#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <set>

struct ContigInfo {
  int16_t color; // Color ID of contig
  std::string s; // Holds more sequence info
  uint8_t rule; // How to process the sequence when found
  bool fwd; // Strandedness
  std::set<int16_t> colors_found; // How many colors (sequence files) it's found in
};

struct TrieNode {
  u_map_<char, TrieNode*> children;
  bool isEndOfWord;
  TrieNode* fail;
  std::vector<std::string> matchedWords;
  std::vector<int> positions;

  TrieNode() : isEndOfWord(false), fail(nullptr) {}

  ~TrieNode() {
    for (auto& kvp : children)
      delete kvp.second;
  }
};

class AhoCorasick {
private:
  TrieNode* root;
  int dictionary_index;
  bool initialized;

  TrieNode* createNewNode();
  void insert(const std::string& word, int index);
  void buildFailureLinks();
  void search(const char* corpus, size_t len, std::vector<ContigInfo*>& info_vec);

public:
  AhoCorasick();
  ~AhoCorasick();
  AhoCorasick( const AhoCorasick& ) = delete;
  AhoCorasick& operator=( const AhoCorasick& ) = delete;
  void searchInCorpus(const char* corpus, size_t len, std::vector<ContigInfo*>& info_vec);
  void add(const std::string& contig, uint16_t color);
  void init();
  u_map_<std::string, ContigInfo> infomap;
};

#endif // KURE_STRINGSEARCH_H

