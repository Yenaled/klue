#ifndef KURE_STRINGSEARCH_H
#define KURE_STRINGSEARCH_H

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>

struct ContigInfo {
  uint16_t color; // Color ID
  std::string s; // Holds more sequence info
  uint8_t rule; // How to process the sequence when found
  bool fwd; // Strandedness
  // ATAT = palindrome; could be ATATNCGCA or could be TGCGNATAT
};

struct TrieNode {
  std::unordered_map<char, TrieNode*> children;
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
  std::unordered_map<std::string, ContigInfo> infomap;
  int dictionary_index;
  bool initialized;

  TrieNode* createNewNode();
  void insert(const std::string& word, int index);
  void buildFailureLinks();
  void search(const char* corpus, size_t len);

public:
  AhoCorasick();
  ~AhoCorasick();
  void searchInCorpus(const char* corpus, size_t len);
  void add(const std::string& contig, uint16_t color);
  void init();
};

#endif // KURE_STRINGSEARCH_H

