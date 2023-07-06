#include "StringSearch.h"
#include <stdexcept>

std::string revcomp(const std::string s) {
  std::string r(s);
  std::transform(s.rbegin(), s.rend(), r.begin(), [](char c) {
    switch(c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'a': return 'T';
    case 'c': return 'G';
    case 'g': return 'C';
    case 't': return 'A';
    default: return 'N';
    }
    return 'N';
  });
  return r;
}

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

            auto it = infomap.find(word);
            if (it != infomap.end()) {
              auto info = it->second;
              if (strncmp(corpus+position, word.c_str(), word.length()) != 0) { // Sanity check
                throw std::runtime_error("String search corrupted; discrepancy with position");
              }
              // Debug:
              std::cout << "Word: " << word << ", Position: " << position << ", Color: " << info.color << ", String: " << info.s << ", Strand: " << std::to_string(info.fwd) << std::endl;
              if (info.fwd) {
                int middle_len = info.rule;
                int middle_pos = position+word.length();
                size_t end_pos = middle_pos+middle_len+info.s.length();
                if (end_pos <= len) {
                  if (middle_len == 0 || strncmp(corpus+(end_pos-info.s.length()), info.s.c_str(), info.s.length()) == 0) {
                    // Check if the middle is not zero
                    const char* c = corpus+middle_pos;
                    bool non_ATCG = false;
                    for (size_t i_c = 0; i_c < middle_len; c++, i_c++) {
                      // Debug:
                      // TODO: Check if a non-ATCG exists
                      if (*c != 'A' && *c != 'T' && *c != 'C' && *c != 'G' && *c != 'a' && *c != 't' && *c != 'c' && *c != 'g') {
                        non_ATCG = true;
                      }
                      std::cout << (*c);
                    }
                    // Debug:
                    if (!non_ATCG) {
                      std::cout << "\n^Success FWD match" << std::endl;
                    }
                    // TODO: Add to struct in hashmap
                  }
                }
              } else {
                // TODO: Handle the rev
              }
            } else {
              throw std::runtime_error("String search corrupted; discrepancy with hash map");
            }
          }
        }
        temp = temp->fail;
      }
    } else {
      curr = root;
    }
  }
}

AhoCorasick::AhoCorasick() {
  root = createNewNode();
  dictionary_index = 0;
  initialized = false;
}

void AhoCorasick::add(const std::string& contig, uint16_t color) {
  if (initialized) {
    throw std::runtime_error("String search already initialized; cannot add");
    return;
  }
  int flank;
  if (contig.size() % 2 == 0) {
    flank = contig.size() / 2;
  } else {
    flank = (contig.size() - 1) / 2;
  }
  std::string word = contig.substr(0, flank);
  ContigInfo info;
  info.color = color;
  info.s = contig.substr(contig.length()-flank);
  info.rule = (contig.size() % 2 != 0);
  info.fwd = true;
  infomap[word] = info;
  insert(word, dictionary_index + 1); // Adjust index by 1 to avoid multiplication by 0
  dictionary_index++;
  word = revcomp(word);
  info.s = revcomp(info.s);
  info.fwd = false;
  infomap[word] = info;
  insert(word, dictionary_index + 1); // Insert the reverse complement
  dictionary_index++;
}

void AhoCorasick::init() {
  if (initialized) {
    throw std::runtime_error("String search already initialized; cannot initialize again");
    return;
  }
  buildFailureLinks();
  initialized = true;
}

AhoCorasick::~AhoCorasick() {
  delete root;
}

void AhoCorasick::searchInCorpus(const char* corpus, size_t len) {
  if (!initialized) {
    throw std::runtime_error("String search not initialized; cannot do search");
    return;
  }
  search(corpus, len);
}
