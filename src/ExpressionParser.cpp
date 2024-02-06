#include "ExpressionParser.h"

ExpressionParser::ExpressionParser(const std::string& input) : input(input) { tokenize(this->input); }

Node* ExpressionParser::parse() {
    size_t index = 0;
    return parseExpression(index);
}

std::vector<Token> ExpressionParser::tokenize(const std::string& input) {
    tokens.clear();
    int open_pars = 0; // Track open parentheses
    for (char c : input) {
        switch (c) {
        case 'U':
            tokens.push_back(Token(UNION, c));
            break;
        case 'I':
            tokens.push_back(Token(INTERSECT, c));
            break;
        case '\\':
            tokens.push_back(Token(DIFFERENCE, c));
            break;
        case '(':
            tokens.push_back(Token(OPEN_PAR, c));
            open_pars++;
            break;
        case ')':
            tokens.push_back(Token(CLOSE_PAR, c));
            open_pars--;
            break;
        case 'N':
            tokens.push_back(Token(NAND, c));
            break;
        case 'X':
            tokens.push_back(Token(XOR, c));
            break;
        default:
            tokens.push_back(Token(VALUE, c));
            break;
        }
    }
    if (open_pars != 0) {
        throw std::runtime_error("[WARNING] Unmatched parentheses in input.");
    }
    return tokens;
}


Node* ExpressionParser::parsePrimary(size_t& index) {
    if (index >= tokens.size()) return nullptr;
    Node* node = nullptr;
    if (tokens[index].type == VALUE) {
        node = new Node(tokens[index].value);
        index++;
    }
    else if (tokens[index].type == OPEN_PAR) {
        index++; // consume '('
        node = parseExpression(index);
        if (index < tokens.size() && tokens[index].type == CLOSE_PAR) {
            index++; // consume ')'
        }
    }
    return node;
}

Node* ExpressionParser::parseExpression(size_t& index) {
    if (index >= tokens.size()) return nullptr;
    Node* left = parsePrimary(index);
    while (index < tokens.size() && (tokens[index].type == UNION || tokens[index].type == INTERSECT || tokens[index].type == DIFFERENCE || tokens[index].type == NAND || tokens[index].type == XOR)) {
        Node* operation = new Node(tokens[index].value);
        index++;
        Node* right = parsePrimary(index);
        operation->left = left;
        operation->right = right;
        left = operation;
    }
    return left;
}
