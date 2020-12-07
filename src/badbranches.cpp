#include "badbranches.hpp"

#include <ROOT/RDataFrame.hxx>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>
#include <set>

using namespace std;

namespace
{
std::vector<std::string> get_branches(const std::string_view& f)
{
    ROOT::RDataFrame d{"Events", f};
    return d.GetColumnNames();
}
} // namespace

void badbranches()//int argc, char* argv[])
{
    const std::vector<std::string> files{
        "/home/data/nanoAOD_2017/mumuRun2017B/1CEF9BDB-BC44-E811-996F-801844DF001C.root",
        "/home/data/nanoAOD_2017/mumuRun2017B/324D0DE2-BC44-E811-BA69-FA163EFC9F83.root",
        "/home/data/nanoAOD_2017/mumuRun2017B/6CC2B88E-C444-E811-87B6-1866DAEA7E28.root",
        "/home/data/nanoAOD_2017/mumuRun2017B/7E3D164F-BF44-E811-88CE-FA163E96C69B.root",
        "/home/data/nanoAOD_2017/mumuRun2017B/9EB58CB8-1C47-E811-9382-FA163E67A014.root",
        "/home/data/nanoAOD_2017/mumuRun2017B/A41543AB-C444-E811-9599-FA163E9592B6.root",
        "/home/data/nanoAOD_2017/mumuRun2017B/B405DFEE-BC44-E811-A346-FA163EA1976B.root",
        "/home/data/nanoAOD_2017/mumuRun2017B/BC7C8A2C-BE44-E811-93C8-FA163E55A87B.root",
        "/home/data/nanoAOD_2017/muRun2017B/0415ADA8-C144-E811-B740-0CC47A78A408.root",
        "/home/data/nanoAOD_2017/muRun2017B/08409896-EB44-E811-95DC-0025905B8572.root",
        "/home/data/nanoAOD_2017/muRun2017B/0E6B1434-C444-E811-AEF3-0025905A60F8.root",
        "/home/data/nanoAOD_2017/muRun2017B/1254CE0A-CE44-E811-8886-0025905A6060.root",
        "/home/data/nanoAOD_2017/muRun2017B/1656E658-C644-E811-9BFE-0025905B85FC.root",
        "/home/data/nanoAOD_2017/muRun2017B/16E2A292-E444-E811-9D25-0CC47A4D75EC.root",
        "/home/data/nanoAOD_2017/muRun2017B/18C1E7FA-D844-E811-A30A-0025905B858A.root",
        "/home/data/nanoAOD_2017/muRun2017B/1A50C5CC-D944-E811-89B5-0CC47A4D75F8.root",
        "/home/data/nanoAOD_2017/muRun2017B/1A6B1E36-0645-E811-AE6A-0025905B85CA.root",
        "/home/data/nanoAOD_2017/muRun2017B/1E5D58E3-C944-E811-AD70-0CC47A7C347A.root",
        "/home/data/nanoAOD_2017/muRun2017B/2A78DDCA-D144-E811-9195-0CC47A78A3B4.root",
        "/home/data/nanoAOD_2017/muRun2017B/309DFAB5-DA44-E811-B994-0025905B8604.root",
        "/home/data/nanoAOD_2017/muRun2017B/387C2232-E144-E811-878D-0025905A48D6.root",
        "/home/data/nanoAOD_2017/muRun2017B/3A878CDD-D244-E811-8726-0025905A612A.root",
        "/home/data/nanoAOD_2017/muRun2017B/48538D43-D844-E811-A947-0CC47A4D7632.root",
        "/home/data/nanoAOD_2017/muRun2017B/4CDAADA7-C844-E811-B132-0CC47A4D769A.root",
        "/home/data/nanoAOD_2017/muRun2017B/5241FE1A-E146-E811-A220-0CC47A4D75F0.root",
        "/home/data/nanoAOD_2017/muRun2017B/580DD3AD-D644-E811-826C-003048FFD7AA.root",
        "/home/data/nanoAOD_2017/muRun2017B/5AB1AECE-CF44-E811-8CA1-0CC47A7C3572.root",
        "/home/data/nanoAOD_2017/muRun2017B/5AED7CD2-D544-E811-87C6-003048FFD7A4.root",
        "/home/data/nanoAOD_2017/muRun2017B/5CF71049-D344-E811-8C01-0CC47A4D7606.root",
        "/home/data/nanoAOD_2017/muRun2017B/5E712F35-D444-E811-B51F-0CC47A78A4BA.root",
        "/home/data/nanoAOD_2017/muRun2017B/60BDBEAD-CE44-E811-8158-0CC47A4C8E56.root",
        "/home/data/nanoAOD_2017/muRun2017B/6CFB1049-D344-E811-8D28-0CC47A4D7606.root",
        "/home/data/nanoAOD_2017/muRun2017B/725870E1-D244-E811-A932-0025905A6092.root",
        "/home/data/nanoAOD_2017/muRun2017B/74FE61A3-D644-E811-A84C-0CC47A4C8F30.root",
        "/home/data/nanoAOD_2017/muRun2017B/78C0F886-D744-E811-9191-0025905B85D0.root",
        "/home/data/nanoAOD_2017/muRun2017B/78C6A5D4-DB44-E811-8209-0CC47A4D7630.root",
        "/home/data/nanoAOD_2017/muRun2017B/8E402B35-D444-E811-B415-0CC47A78A4BA.root",
        "/home/data/nanoAOD_2017/muRun2017B/982289F1-1345-E811-BD4B-001F29082E68.root",
        "/home/data/nanoAOD_2017/muRun2017B/98FF1773-D144-E811-A7FC-0CC47A745282.root",
        "/home/data/nanoAOD_2017/muRun2017B/9C593350-CB44-E811-A980-0025905B85DA.root",
        "/home/data/nanoAOD_2017/muRun2017B/A00B80C0-D044-E811-8E00-0025905A60B6.root",
        "/home/data/nanoAOD_2017/muRun2017B/AE5CBCF5-1345-E811-A638-00269E95ACFC.root",
        "/home/data/nanoAOD_2017/muRun2017B/B416A542-1445-E811-9392-A0369FC5B844.root",
        "/home/data/nanoAOD_2017/muRun2017B/B68812B8-D644-E811-909B-0025905A48D6.root",
        "/home/data/nanoAOD_2017/muRun2017B/B6B99584-BF44-E811-BD70-0025905A609A.root",
        "/home/data/nanoAOD_2017/muRun2017B/B8675ED2-D544-E811-83AB-0025905A613C.root",
        "/home/data/nanoAOD_2017/muRun2017B/BA1FD7CE-E546-E811-8F21-0025905B861C.root",
        "/home/data/nanoAOD_2017/muRun2017B/BA27E502-BC44-E811-AB03-FA163E4E91B6.root",
        "/home/data/nanoAOD_2017/muRun2017B/C8E85408-CD44-E811-9B17-0025905B8600.root",
        "/home/data/nanoAOD_2017/muRun2017B/CC51608A-D744-E811-B6C7-0025905B8612.root",
        "/home/data/nanoAOD_2017/muRun2017B/CC712F35-D444-E811-9352-0CC47A78A4BA.root",
        "/home/data/nanoAOD_2017/muRun2017B/CCF9C7E3-1345-E811-B6A6-FA163EBDA738.root",
        "/home/data/nanoAOD_2017/muRun2017B/CE48AC39-DE44-E811-B638-0025905A60DE.root",
        "/home/data/nanoAOD_2017/muRun2017B/D2234164-BD44-E811-9741-0CC47A7C34EE.root",
        "/home/data/nanoAOD_2017/muRun2017B/DA269143-D844-E811-A9ED-0CC47A4D7632.root",
        "/home/data/nanoAOD_2017/muRun2017B/E216A8DF-1345-E811-8EA3-0025905A48F0.root",
        "/home/data/nanoAOD_2017/muRun2017B/EE345672-D144-E811-ABF5-0CC47A7C35A8.root",
        "/home/data/nanoAOD_2017/muRun2017B/F205F211-D544-E811-94DD-0CC47A4D7630.root",
        "/home/data/nanoAOD_2017/muRun2017B/F284C01E-D544-E811-BD7E-0025905A6122.root"};

    std::vector<std::vector<std::string>> file_branches{};

    std::transform(files.begin(), files.end(), std::back_inserter(file_branches), get_branches);

    // Make union of all brach names
    std::set<std::string_view> u{};
    for (const auto& strings: file_branches)
    {
        std::copy(strings.begin(), strings.end(), std::inserter(u, u.end()));
    }

    std::cout << "Branches which do not exist in all files:" << '\n';
    for (const auto& branch: u)
    {
        if (std::any_of(file_branches.cbegin(), file_branches.cend(), [&](const std::vector<std::string>& bs) { return std::find(bs.begin(), bs.end(), branch) == bs.end(); }))
        {
            std::cout << '\t' << branch << '\n';
        }
    }
    std::cout << std::flush;
}
