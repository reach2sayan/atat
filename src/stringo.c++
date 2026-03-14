#include "stringo.h"

/*
istream& operator >> (istream &file, std::string &str) {
        static char buf[1000];
        file.get(buf,1000);
        return file;
}
*/

Real to_real(const std::string &s) {
  istringstream str(s);
  Real r;
  str >> r;
  return r;
}
