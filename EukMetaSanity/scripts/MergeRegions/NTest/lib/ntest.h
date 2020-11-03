#ifndef NTEST_H
#define NTEST_H

#define RESET "\033[0m"
#define RED "\033[31m"   /* Red */
#define GREEN "\033[32m" /* Green */
#define CYAN "\033[36m"  /* Cyan */

#include <chrono>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

struct NTestExpectError : std::exception {};

class NTest {
public:
    // Instantiate with 0 tests passed
    // Determine if exception will end program at test failure
    NTest(bool _throw = false) : _passed(0), throw_on_fail(_throw){};
    // Return the number of tests that have passed
    int passed() const { return _passed; }
    // Return the total number of tests stored
    int total() const { return tests.size(); }
    // Run each test and update pass status
    virtual void test() {
        // Check if exception on error
        std::function<void(bool)> check_method = [](bool) { return; };
        if (throw_on_fail)
            check_method = [](bool _check) {
                if (!_check)
                    throw std::exception();
            };
        for (auto& test : tests) {
            // Name of test
            std::string& test_name = std::get<0>(test);
            std::cerr << CYAN << "***** Begin " << test_name << " *****" << RESET << std::endl;
            // Run test and output result
            bool result;
            double _time = NTest::experiment([&result, &test]() { result = std::get<1>(test)(); });
            if (!result) {
                std::cerr << RED << test_name << " failed in " << _time << "ms" << RESET << std::endl;
                // Raise error and exit
                check_method(false);
            } else {
                std::cerr << GREEN << test_name << " passed in " << _time << "ms" << RESET << std::endl;
                _passed++;
            }
            std::cerr << std::endl;
        }
    }
    // Generate a random string of size length
    static std::string random_string(size_t length) {
        auto randchar = []() -> char {
            const char charset[] = "0123456789"
                                   "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                   "abcdefghijklmnopqrstuvwxyz";
            const size_t max_index = (sizeof(charset) - 1);
            return charset[rand() % max_index];
        };
        std::string str(length, 0);
        std::generate_n(str.begin(), length, randchar);
        return str;
    }
    // Preferred method for adding a test
    void make_test(std::string name, std::function<bool()> test) { push(std::make_tuple(name, test)); }
    // Simple comparison method for test result, default is to
    // compare against `true`
    static bool expect(bool a, bool _throw = false, bool b = true) {
        bool _true = (a == b);
        if (!_true && _throw)
            throw NTestExpectError();
        return _true;
    }
    // Time the execution of a particular function
    static double experiment(std::function<void()> test) {
        auto start = std::chrono::high_resolution_clock::now();
        test();
        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = stop - start;
        return duration.count();
    }

protected:
    // Add tuple to list of tests to run
    void push(std::tuple<std::string, std::function<bool()>> test) { tests.push_back(test); }

private:
    std::vector<std::tuple<std::string, std::function<bool()>>> tests;
    int _passed;
    bool throw_on_fail;
};
#endif