#ifndef EXPRESSION_PARSER_H
#define EXPRESSION_PARSER_H

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

enum TokenType {
    VALUE,
    UNION,
    INTERSECT,
    DIFFERENCE,
    NAND,
    XOR,
    OPEN_PAR,
    CLOSE_PAR
};

struct Token {
    TokenType type;
    char value;
    Token(TokenType t, char v) : type(t), value(v) {}
};

struct Node {
    char value;
    Node* left;
    Node* right;

    Node(char v) : value(v), left(nullptr), right(nullptr) {}
};

class ExpressionParser {
public:
    ExpressionParser(const std::string& input);
    Node* parse();
    std::vector<Token> tokenize(const std::string& input);

private:
    std::string input;
    std::vector<Token> tokens;

    Node* parsePrimary(size_t& index);
    Node* parseExpression(size_t& index);
};

#endif // EXPRESSION_PARSER_H
