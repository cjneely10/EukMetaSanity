#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>
#include <string>

class FileException : public std::runtime_error {
public:
    FileException(const std::string& message) : std::runtime_error(message) {}
    virtual ~FileException() throw() {}
};

class StrandComparisonException : public std::runtime_error {
public:
    StrandComparisonException(const std::string& message) : std::runtime_error(message) {}
    virtual ~StrandComparisonException() throw() {}
};

class RangeException : public std::runtime_error {
public:
    RangeException(const std::string& message) : std::runtime_error(message) {}
    virtual ~RangeException() throw() {}
};

#endif