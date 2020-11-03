#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include <sstream>
#include <stack>

// Type of argument - either fill or set bool by flag
namespace ArgParse {
enum class Type { DEFAULT, PARSE_TRUE, PARSE_FALSE, FLAG };
template<typename T, typename Validator>
bool validate(T& ref, Validator& validator) {
    bool v = validator(ref);
    if (!v) {
        std::cerr << std::endl << ref << " error: " << validator.ERR << std::endl;
    }
    return v;
}

struct DefaultValidator {
    template<typename T>
    bool operator()(T& ref) {
        return true;
    }
    const char* ERR = "Failed to validate!";
} default_validator;

// Simple struct for holding items from argv
struct Token {
    bool passed;
    bool required;
    const char* description;
    std::initializer_list<const char*> identifiers;
    Type type;
    Token(bool req, const char* descr, std::initializer_list<const char*> c, Type t)
            : passed(false), required(req), description(descr), identifiers(c), type(t) {}
    inline void print_identifiers(std::ostream& o) const {
        std::initializer_list<const char*>::const_iterator ft = identifiers.begin();
        o << *ft;
        ++ft;
        while (ft != identifiers.end()) {
            o << "/" << *ft;
            ++ft;
        }
        if (type == Type::FLAG)
            o << " "
              << "<arg>";
    }
    inline void print(std::ostream& o) {
        // Output non-required with surrounding brackets
        o << " ";
        if (!required) {
            o << "[";
            print_identifiers(o);
            o << "]";
        }
        // Output required
        else {
            print_identifiers(o);
        }
    }
};

class ArgParse {
public:
    // Create ArgParse and parse system arguments to vector
    ArgParse(int& argc, char** argv, const char* d = "ArgParse program!") : _description(d), help_requested(false) {
        // Parse all except for name of program
        std::stack<std::string> tmp;
        _name = argv[0];
        for (int i = 1; i < argc; i++) {
            tmp.push(std::string(argv[i]));
            // Check if help flag was passed
            if (tmp.top() == "-h" || tmp.top() == "--help") {
                help_requested = true;
            }
        }
        // Push to object storage in better order
        while (!tmp.empty()) {
            args.push(tmp.top());
            tmp.pop();
        }
        // Auto-set help if nothing is passed
        if (args.size() == 0)
            help_requested = true;
    }
    // Pass reference to variable to fill
    template<typename T, typename Validator = DefaultValidator>
    void
    refer(T& ref,
          std::initializer_list<const char*> flags,
          const char* description,
          Type ptype = Type::DEFAULT,
          bool required = true,
          Validator& validator = default_validator);
    // Call to parse results and display help if failes
    void parse_or_exit();

private:
    static const std::string ERR;
    // Vector of arguments from argv
    std::stack<std::string> args;
    // Parsed tokens with pass/fail status
    std::stack<Token> tokens;
    // Short description
    std::string _description;
    // Check if help was requested
    bool help_requested;
    // Name of program
    std::string _name;
    // Help screen
    void get_help(const std::string& err_msg);
    void usage_statement();
    // Parse positional argument
    template<typename T, typename Validator>
    void positional(T& ref, Token& data, Validator& val);
    // Parse flagged argument
    template<typename T, typename Validator>
    void flagged(T& ref, Token& data, Type& ptype, std::initializer_list<const char*>& it, Validator& val);
    // Exception on parsing
    void parse_fail(const std::string c = "Error parsing passed arguments") {
        get_help(c);
        exit(0);
    }
};

const std::string ArgParse::ERR = "Unable to properly parse arguments";

void ArgParse::parse_or_exit() {
    std::stack<Token> tmp;
    while (!tokens.empty()) {
        tmp.push(tokens.top());
        tokens.pop();
    }
    tokens = tmp;
    if (help_requested)
        parse_fail("");
    if (tokens.size() > 0 && !args.empty())
        parse_fail(ArgParse::ERR);
    Token top = tokens.top();
    tokens.pop();
    std::stack<Token> iter;
    while (!tokens.empty()) {
        iter.push(top);
        tokens.pop();
        if (top.required && !top.passed) {
            while (!iter.empty()) {
                tokens.push(iter.top());
                iter.pop();
            }
            parse_fail(std::string("Argument failed to parse: ") + *top.identifiers.begin());
        }
        top = tokens.top();
    }
}

template<typename T, typename Validator>
void ArgParse::refer(
        T& ref,
        std::initializer_list<const char*> flags,
        const char* description,
        Type ptype,
        bool required,
        Validator& validator) {
    // Create token to store
    Token data(required, description, flags, ptype);
    switch (ptype) {
    case Type::DEFAULT:
        positional(ref, data, validator);
        break;
    default:
        flagged(ref, data, ptype, flags, validator);
    }
}

template<typename T, typename Validator>
void ArgParse::positional(T& ref, Token& data, Validator& validator) {
    std::stack<std::string> tmp;
    // Check args to see if one of the flags is present
    while (!args.empty()) {
        // Start search at top of args
        std::string top = args.top();
        args.pop();
        tmp.push(top);
        // Check non-flagged arg
        if (top.at(0) != '-') {
            // Store value
            std::stringstream(top) >> ref;
            // Store as passed
            data.passed = validate(ref, validator);
            tmp.pop();
            break;
        }
    }
    // Add to tokens stack
    tokens.push(data);
    while (!tmp.empty()) {
        args.push(tmp.top());
        tmp.pop();
    }
}

template<typename T, typename Validator>
void ArgParse::flagged(
        T& ref, Token& data, Type& ptype, std::initializer_list<const char*>& flags, Validator& validator) {
    std::stack<std::string> tmp;
    std::initializer_list<const char*>::iterator it = flags.begin();
    while (it != flags.end()) {
        while (!args.empty()) {
            std::string top = args.top();
            args.pop();
            tmp.push(top);
            // Value found
            if (top == *it) {
                // Consume value
                tmp.pop();
                // Store as passed
                data.passed = validate(ref, validator);
                // Parse value
                switch (ptype) {
                // Store value following flag
                case Type::FLAG:
                    // Confirm next value is present and store
                    if (args.empty())
                        parse_fail(ArgParse::ERR);
                    // throw error if failed
                    if ((std::stringstream(args.top()) >> ref).fail())
                        parse_fail(top + " not provided");
                    // Pop off this next value, too
                    args.pop();
                    break;
                // Store flag to make value false
                case Type::PARSE_FALSE:
                    ref = false;
                    break;
                // Store flag to make value true
                case Type::PARSE_TRUE:
                    ref = true;
                    break;
                default:
                    parse_fail(ArgParse::ERR);
                }
                // Add to tokens stack
                tokens.push(data);
                // Return unused args back to stack
                while (!tmp.empty()) {
                    args.push(tmp.top());
                    tmp.pop();
                }
                return;
            }
        }
        ++it;
    }
    // Add to tokens stack
    tokens.push(data);
    while (!tmp.empty()) {
        args.push(tmp.top());
        tmp.pop();
    }
}

void ArgParse::usage_statement() {
    std::cerr << "Usage " << _name;
    Token& it = tokens.top();
    std::stack<Token> tmp;
    tmp.push(it);
    tokens.pop();
    while (!tokens.empty()) {
        it.print(std::cerr);
        // Move to next value
        it = tokens.top();
        tmp.push(tokens.top());
        tokens.pop();
    }
    it.print(std::cerr);
    while (!tmp.empty()) {
        tokens.push(tmp.top());
        tmp.pop();
    }
    std::cerr << std::endl;
}

void ArgParse::get_help(const std::string& err_msg) {
    // Usage statement
    usage_statement();
    std::cerr << std::endl;
    // Description of program
    if (_description.size() > 0)
        std::cerr << _description << "\n" << std::endl;

    Token& it = tokens.top();
    std::stack<Token> tmp;
    tmp.push(tokens.top());
    tokens.pop();
    while (!tokens.empty()) {
        it.print(std::cerr);
        std::cerr << std::endl;
        std::cerr << "\t\t\t" << it.description << std::endl;
        it = tokens.top();
        tmp.push(tokens.top());
        tokens.pop();
    }
    it.print(std::cerr);
    std::cerr << std::endl;
    std::cerr << "\t\t\t" << it.description << std::endl;
    if (err_msg.size() > 0)
        std::cerr << err_msg << std::endl;
}
}  // namespace ArgParse

#endif