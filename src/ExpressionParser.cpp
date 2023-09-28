#include "ExpressionParser.h"

ExpressionParser::ExpressionParser(const std::string& input) : input(input) { tokenize(this->input); }

Node* ExpressionParser::parse() {
    size_t index = 0;
    return parseExpression(index);
}

std::vector<Token> ExpressionParser::tokenize(const std::string& input) {
    tokens.clear();
    for (char c : input) {
        switch (c) {
        case 'U':
            tokens.push_back(Token(UNION, c));
            break;
        case 'I':
            tokens.push_back(Token(INTERSECT, c));
            break;
        case '(':
            tokens.push_back(Token(OPEN_PAR, c));
            break;
        case ')':
            tokens.push_back(Token(CLOSE_PAR, c));
            break;
        default:  // A, B, C, etc.
            tokens.push_back(Token(VALUE, c));
            break;
        }
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

    while (index < tokens.size() && (tokens[index].type == UNION || tokens[index].type == INTERSECT)) {
        Node* operation = new Node(tokens[index].value);
        index++; // consume operation

        Node* right = parsePrimary(index);

        operation->left = left;
        operation->right = right;
        left = operation;  // the current operation becomes the new left operand for the next loop iteration
    }

    return left;
}
