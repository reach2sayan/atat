#ifndef __STRINGO_H__
#define __STRINGO_H__

#include "misc.h"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

inline int strcmp(const std::string &lhs, const char *rhs) {
  return lhs.compare(rhs == nullptr ? "" : rhs);
}

inline int strcmp(const char *lhs, const std::string &rhs) {
  return std::strcmp(lhs == nullptr ? "" : lhs, rhs.c_str());
}

inline int strcmp(const std::string &lhs, const std::string &rhs) {
  return lhs.compare(rhs);
}

inline size_t strlen(const std::string &value) { return value.size(); }

inline char *strcpy(char *dest, const std::string &src) {
  return std::strcpy(dest, src.c_str());
}

Real to_real(const std::string &s);

#include "binstream.h"

inline ostream &bin_ostream(ostream &file, const std::string &str) {
  bin_ostream(file, static_cast<int>(str.size()));
  for (char c : str) {
    bin_ostream(file, c);
  }
  return file;
}

inline istream &bin_istream(istream &file, std::string &str) {
  int len;
  bin_istream(file, len);
  str.assign(static_cast<std::string::size_type>(len), '\0');
  for (int i = 0; i < len; i++) {
    bin_istream(file, str[static_cast<std::string::size_type>(i)]);
  }
  return file;
}

#endif
