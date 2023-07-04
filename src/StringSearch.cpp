#include "StringSearch.h"

TrieNode* AhoCorasick::createNewNode() {
  TrieNode* node = new TrieNode;
  return node;
}

void AhoCorasick::insert(const std::string& word, int index) {
  TrieNode* curr = root;
  
  for (char c : word) {
    if (curr->children.find(c) == curr->children.end())
      curr->children[c] = createNewNode();
    
    curr = curr->children[c];
  }
  
  curr->isEndOfWord = true;
  curr->matchedWords.push_back(word);
  curr->positions.push_back(index);
}

void AhoCorasick::buildFailureLinks() {
  std::queue<TrieNode*> q;
  q.push(root);
  
  while (!q.empty()) {
    TrieNode* curr = q.front();
    q.pop();
    
    for (auto& kvp : curr->children) {
      char c = kvp.first;
      TrieNode* child = kvp.second;
      
      TrieNode* failure = curr->fail;
      
      while (failure && failure->children.find(c) == failure->children.end())
        failure = failure->fail;
      
      if (failure)
        child->fail = failure->children[c];
      else
        child->fail = root;
      
      q.push(child);
    }
  }
}

void AhoCorasick::search(const char* corpus, size_t len) {
  TrieNode* curr = root;
  
  for (int i = 0; i < len; ++i) {
    char c = corpus[i];
    
    while (curr && curr->children.find(c) == curr->children.end())
      curr = curr->fail;
    
    if (curr) {
      curr = curr->children[c];
      
      TrieNode* temp = curr;
      while (temp) {
        if (temp->isEndOfWord) {
          for (size_t j = 0; j < temp->matchedWords.size(); ++j) {
            const std::string& word = temp->matchedWords[j];
            int position = i - static_cast<int>(word.length()) + 1;
            
            // Debug:
            //std::cout << "Found occurrence of word: " << word << " at position: " << position << std::endl;
          }
        }
        temp = temp->fail;
      }
    } else {
      curr = root;
    }
  }
}

AhoCorasick::AhoCorasick(const std::vector<std::string>& dictionary) {
  root = createNewNode();
  
  for (int i = 0; i < dictionary.size(); ++i)
    insert(dictionary[i], i + 1);  // Adjust index by 1 to avoid multiplication by 0
  
  buildFailureLinks();
}

AhoCorasick::~AhoCorasick() {
  delete root;
}

void AhoCorasick::searchInCorpus(const char* corpus, size_t len) {
  search(corpus, len);
}

// Sample usage:

/*
 std::vector<std::string> dictionary = { "he", "she", "his", "hers","a","abc","ab" };
 std::string corpus = "ushershehisherabhehehesheisherhabcehe";
 AhoCorasick ac(dictionary);
 ac.searchInCorpus(corpus.c_str(), corpus.length());
 */
