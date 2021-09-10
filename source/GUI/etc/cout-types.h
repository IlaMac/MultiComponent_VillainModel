#ifndef COUT_TYPES_H
#define COUT_TYPES_H

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <map>
#include <set>
#include <utility>
#include <tuple>
#include <string>
#include <optional>
#include <variant>

namespace std {

  template<typename T>
  ostream& operator<< (ostream & out, const vector<T> & v) {
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
      out << v[i];
      if (i != last)
        out << ", ";
    }
    out << "]";
    return out;
  }

  template<typename T, size_t size>
  ostream& operator<< (ostream & out, const array<T, size> & v) {
    out << "[";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
      out << v[i];
      if (i != last)
        out << ", ";
    }
    out << "]";
    return out;
  }

  template<typename T1, typename T2>
  ostream& operator<< (ostream & out, const map<T1, T2> & m) {
    out << "[" << endl;
    for (auto it = m.begin(); it != m.end(); it++) {
      out << "\t" << it->first << ": " << it->second;
      if (next(it, 1) != m.end())
        out << ",";
      out << endl;
    }
    out << "]";
    return out;
  }

  template<typename T1>
  ostream& operator<< (ostream & out, const set<T1> & s) {
    out << "[";
    for (auto it = s.begin(); it != s.end(); it++) {
      out << *it;
      if (next(it, 1) != s.end()) out << ", ";
    }
    out << "]";
    return out;
  }

  template<typename T1, typename T2>
  ostream& operator<< (ostream & out, const pair<T1, T2> & p) {
    out << "(" << p.first << ", " << p.second << ")";
    return out;
  }

  template<typename T>
  ostream& operator<< (ostream & out, const optional<T> & opt) {
    if (opt.has_value()) out << opt.value();
    else                 out << "-";
    return out;
  }

  template <typename T0, typename ... Ts>
  ostream& operator<< (ostream & out, variant<T0, Ts...> const & v) {
    visit([&] (auto && arg) { out << arg;}, v);
    return out;
  }

  namespace aux {
    template<size_t...>
    struct seq{ };


    ////
    //// gen_seq<N>  ->  gen_seq<N-1, N-1>  ->  gen_seq<N-2, N-2, N-1>  ->  gen_seq<N-3, N-3, N-2, N-1>  ->  ...
    ////
    template<size_t N, size_t... Is>                 // one or more size_t's
    struct gen_seq : gen_seq<N-1, N-1, Is...> { };   // "struct xxx : yyy {};" makes xxx inherit yyy


    ////
    //// gen_seq<0, 0, 1, 2, ..., N-1>  ->  seq<0, 1, 2, ..., N-1>
    ////
    template<size_t... Is>   // zero or more size_t's
    struct gen_seq<0, Is...> : seq<Is...> { };


    template<class Tuple, size_t... Is>
    void print_tuple(
      ostream & out,
      Tuple const& t,
      seq<Is...>
    ) {

      ////
      //// int a[] = {get1(), get2()}           will execute get1 before executing get2 (using x = int[] in order to avoid "-Wunused-variable")
      ////
      //// int a[] = {( void( ... ) , 0 )}      use the comma operator to support operations which do not return a proper value
      ////

      using swallow = int[];
      (void) swallow{(void(out << (Is == 0 ? "" : ", ") << get<Is>(t)), 0)...};

    }
  }

  template<typename... Ts>
  ostream& operator<< (
    ostream & out,
    tuple<Ts...> const& t
  ) {
    out << "(";
    aux::print_tuple(out, t, aux::gen_seq<sizeof...(Ts)>());
    return out << ")";
  }



  // static inline
  // ostream& operator<< (ostream & out, const HoppingAmplitude & ampl) {
  //   out << "{ t=" << ampl << " }";
  //   return out;
  // }

  static inline
  ostream& operator<< (ostream & out, const Hop & hop) {
    out << "{ i=" << hop.i << ", ampl=" << hop.ampl << ", dist=" << hop.dist << " }";
    return out;
  }


}


#endif