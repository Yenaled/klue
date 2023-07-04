#ifndef KURE_STRINGSEARCH_H
#define KURE_STRINGSEARCH_H

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>

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

  TrieNode* createNewNode();
  void insert(const std::string& word, int index);
  void buildFailureLinks();
  void search(const char* corpus, size_t len);

public:
  AhoCorasick(const std::vector<std::string>& dictionary);
  ~AhoCorasick();
  void searchInCorpus(const char* corpus, size_t len);
};

#endif // KURE_STRINGSEARCH_H

