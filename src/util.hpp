#ifndef UTIL_HPP
#define UTIL_HPP

#include <algorithm>
#include <set>
#include <string_view>
#include <string>
#include <vector>
#include <ROOT/RDataFrame.hxx>

#include <iostream>

namespace util
{
[[gnu::const]] inline ROOT::Detail::RDF::ColumnNames_t apply_blacklist(ROOT::Detail::RDF::ColumnNames_t branches)
{
    // These bracnhes have been found not to exist in every ROOT file in a
    // given dataset, and can cause the program to crash when a snapshot is
    // made
    static const std::set<std::string_view> blacklist{{"HLT_HcalIsolatedbunch",
                                                       "HLT_L1FatEvents"}};

    branches.erase(std::remove_if(branches.begin(),
                                  branches.end(),
                                  [](const std::string_view& s) { return blacklist.find(s) != blacklist.end(); }),
                   branches.end());
    return branches;
}
} // namespace util


#endif /* UTIL_HPP */
