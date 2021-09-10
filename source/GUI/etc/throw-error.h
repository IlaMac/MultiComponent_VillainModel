#ifndef THROW_ERROR_H
#define THROW_ERROR_H

#include <string>
#include <sstream>

const std::string enterBlue   = "\033[38;5;12m",
                  enterYellow = "\033[38;5;11m",
                  enterGreen  = "\033[38;5;10m",
                  enterRed    = "\033[38;5;196m",
                  enterOrange = "\033[38;5;202m",
                  enterBold   = "\033[1m",
                  resetStyle  = "\033[0m";

////
//// use these macros rather than the class
////
#define throwError()  _ThrowError(__FILE__, __LINE__, true)
#define _throwError() _ThrowError(__FILE__, __LINE__, false)


////
//// pretty much a copy of AtomicWriter
////
class _ThrowError {
  protected:

    const char * file;
    const int    line;
    const bool   _exit;

    std::ostringstream st;
    std::ostream &     stream;

  public:
    _ThrowError (
      const char *   file,
      const int      line,
      const bool     _exit = true,
      std::ostream & s     = std::cout
    ) : file{file},
        line{line} ,
        _exit{_exit},
        stream(s)
    { }

    template <typename T>
    _ThrowError & operator<< (T const& t) {
      st << t;
      return *this;
    }

    _ThrowError & operator<< (
      std::ostream&(*f)(std::ostream&)
    ) {
      st << f;
      return *this;
    }

    ~_ThrowError () {
      // #pragma omp critical   // apparently it should also work without this, but it does not...
      stream << enterRed
             << "-----------------[ ERROR ]-----------------" << std::endl
             << file << ":" << line << std::endl
             << st.str()
             << resetStyle;

      if (_exit) exit(EXIT_FAILURE);
    }
};


#endif