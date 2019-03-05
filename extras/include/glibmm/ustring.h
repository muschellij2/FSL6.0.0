#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>

#ifndef ENABLE_USTRING

#if defined (__GNUC__)
#  define G_STRFUNC     ((const char*) (__PRETTY_FUNCTION__))
#elif defined (G_HAVE_ISO_VARARGS)
#  define G_STRFUNC     ((const char*) (__func__))
#else
#  define G_STRFUNC     ((const char*) ("???"))
#endif

namespace Glib
{
  //typedef std::string ustring;

  class ustring : public std::string {
  public:
    ustring(const char* c) : std::string(c) {}
    ustring() : std::string() {}
    ustring(const std::string& str) : std::string(str) {} 
    template <class InputIterator>
    ustring(InputIterator first, InputIterator last) : std::string(first,last) {}
    size_t bytes() const { return this->size(); }
    template <typename T>
    static std::string format(const T num) { return boost::lexical_cast<std::string>(num); }
    template <typename T, typename P>
      static std::string format(std::ios_base& (& base)(std::ios_base&), P precision,  const T& num ) {
      std::stringstream ss;
      ss << base << precision << num;
      return ss.str();
    }
  };

  
  inline std::string locale_from_utf8(const std::string& mystr) { return mystr; }
  inline std::string convert(const std::string in, const std::string& par1, const std::string& par2) { return in; }
}

//#define g_unichar_isspace(str) std::isspace(str, std::locale("C"))   
inline bool g_unichar_isspace(const char& ch) { return std::isspace(ch, std::locale("C")); }
#define ENABLE_USTRING
#endif
