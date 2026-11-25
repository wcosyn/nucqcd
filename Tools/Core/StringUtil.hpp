/*! \file StringUtil.hpp
 * \brief Small std::string helpers that replace the TString methods used in the kinematics classes and
 * Deuteron (roadmap step 9, phase 3)
 * \author Wim Cosyn
 * \date 22/11/2026
 *
 * \addtogroup MePhys
 * @{
 */
#ifndef STRINGUTIL_HPP
#define STRINGUTIL_HPP

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace nucqcd {

/*! \param s a string \return s in lower case (ASCII), as TString::ToLower */
inline std::string to_lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

/*! \param a first string \param b second string \return whether a and b are equal ignoring case, as
 * !a.CompareTo(b, TString::kIgnoreCase) */
inline bool iequals(const std::string &a, const std::string &b) { return to_lower(a) == to_lower(b); }

/*! \param s a string \param part the substring \return whether s contains part ignoring case, as
 * s.Contains(part, TString::kIgnoreCase) */
inline bool icontains(const std::string &s, const std::string &part) {
  return to_lower(s).find(to_lower(part)) != std::string::npos;
}

/*! Splits a string at a separator, leaving out empty tokens, as TString::Tokenize
 * \param s the string \param sep the separator \return the non-empty tokens */
inline std::vector<std::string> split(const std::string &s, char sep) {
  std::vector<std::string> tokens;
  std::string::size_type start = 0;
  while (start <= s.size()) {
    const std::string::size_type end = std::min(s.find(sep, start), s.size());
    if (end > start) tokens.push_back(s.substr(start, end - start));
    start = end + 1;
  }
  return tokens;
}

}  // namespace nucqcd
/*! @} */
#endif
