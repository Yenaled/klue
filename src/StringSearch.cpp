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

void AhoCorasick::processWord(const char* corpus, size_t len, const std::string& word, int position, std::vector<ContigInfo*>& info_vec) {
  const std::string word_r = revcomp(word);
  bool lex = (word < word_r);
  auto it = infomap.find(lex ? word : word_r);
  if (it != infomap.end()) {
    auto& info = it->second;
    if (strncmp(corpus+position, word.c_str(), word.length()) != 0) { // Sanity check
      throw std::runtime_error("String search corrupted; discrepancy with position");
    }
    int flank;
    auto& contig = info.s;
    if (contig.size() % 2 == 0) {
      flank = contig.size() / 2;
    } else {
      flank = (contig.size() - 1) / 2;
    }
    auto info_s = contig.substr(contig.length()-flank);
    // Debug:
    std::cout << "Word: " << word << ", Position: " << position << ", Color: " << info.color << ", String: " << info_s << ", Strand: " << std::to_string(info.fwd) << std::endl;
    int middle_len = info.rule;
    bool is_palindrome = (word == word_r);
    bool found_fwd = (word == (lex ? word : word_r));
    bool success = false;
    if ((found_fwd && info.fwd) || (!found_fwd && !info.fwd) || is_palindrome) {
      int middle_pos = position+word.length();
      int end_pos = middle_pos+middle_len+info_s.length();
      if (end_pos <= len) {
        if (strncmp(corpus+(end_pos-info_s.length()), info_s.c_str(), info_s.length()) == 0) {
          const char* c = corpus+middle_pos;
          bool non_ATCG = false;
          for (size_t i_c = 0; i_c < middle_len; c++, i_c++) {
            // Debug:
            if (*c != 'A' && *c != 'T' && *c != 'C' && *c != 'G' && *c != 'a' && *c != 't' && *c != 'c' && *c != 'g') {
              non_ATCG = true;
            }
            std::cout << (*c); // Debug
          }
          // Debug:
          if (!non_ATCG) {
            success = true;
            std::cout << "\n^Success FWD match" << std::endl; // Debug
          }
        }
      }
    }
    // TODO: implement other features (other flanks); etc.
    if ((!((found_fwd && info.fwd) || (!found_fwd && !info.fwd)) || is_palindrome)) {
      int start_pos = position-middle_len-info_s.length();
      // Either we found it in the wrong way but it matches contig; e.g. we found TTTT, but AAAA is hashmap'd and is fwd which means we have to rev
      // Or we found it in the right way but it doesn't match contig; e.g. we found TTTT and TTTT is hashmap'd but it's not fwd
      // Either way: We need to reverse complement the other flank!
      if (start_pos >= 0) {
        std::string s_rev = revcomp(info_s);
        if (strncmp(corpus+start_pos, s_rev.c_str(), info_s.length()) == 0) {
          const char* c = corpus+start_pos+info_s.length();
          bool non_ATCG = false;
          for (size_t i_c = 0; i_c < middle_len; c++, i_c++) {
            if (*c != 'A' && *c != 'T' && *c != 'C' && *c != 'G' && *c != 'a' && *c != 't' && *c != 'c' && *c != 'g') {
              non_ATCG = true;
            }
            std::cout << (*c); // Debug
          }
          if (!non_ATCG) {
            success = true;
            std::cout << "\n^Success REV match" << std::endl; // Debug
          }
        }
      }
    }
    if (success) {
      int sz = info_vec.size();
      if (sz >= 4096) {
        info_vec.reserve(sz*1.2); // grow slowly in capacity
      }
      info_vec.push_back(&info);
    }
  } else {
    throw std::runtime_error("String search corrupted; discrepancy with hash map");
  }
}

void AhoCorasick::search(const char* corpus, size_t len, std::vector<ContigInfo*>& info_vec) {
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
            processWord(corpus, len, word, position, info_vec);
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
  info.s = contig;
  info.rule = (contig.size() % 2 != 0);
  insert(word, dictionary_index + 1); // Adjust index by 1 to avoid multiplication by 0
  dictionary_index++;
  std::string word_r = revcomp(word);
  bool lex = (word < word_r); // We want to index the canonical (i.e. lexicographically smaller) word
  info.fwd = lex; // If we're storing the sequence as-is (i.e. the non-reverse-complemented sequence is the canonical one)
  infomap[lex ? word : word_r] = info; // Index into hash map
  insert(word_r, dictionary_index + 1); // Insert the reverse complement
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

void AhoCorasick::searchInCorpus(const char* corpus, size_t len, std::vector<ContigInfo*>& info_vec) {
  if (!initialized) {
    throw std::runtime_error("String search not initialized; cannot do search");
    return;
  }
  search(corpus, len, info_vec);
}
