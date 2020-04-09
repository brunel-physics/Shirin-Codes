#include "analyse.hpp"
#include <algorithm>
#include "TLorentzVector.h"
#include "sf.hpp"
#include <TCanvas.h>
#include <TText.h>
#include <THStack.h>
#include <TH1D.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include "eval_complex.hpp"
#include <TStyle.h>
#include "tdrstyle.C"
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include <ROOT/RResultPtr.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vdt/vdtMath.h>
#include <numeric>



using namespace std;

using doubles = ROOT::VecOps::RVec<double>;
using floats = ROOT::VecOps::RVec<float>;
using ints = ROOT::VecOps::RVec<int>;
using bools = ROOT::VecOps::RVec<bool>;
using chars = ROOT::VecOps::RVec<UChar_t>; // aka 1 byte ints
using strings = ROOT::VecOps::RVec<string>;



/*class CSVReader
{
	std::string fileName;
	std::string delimeter;

public:
	CSVReader(std::string filename, std::string delm = ",") :
			fileName(filename), delimeter(delm)
	{ }

	// Function to fetch data from a CSV File
	std::vector<std::vector<std::string> > getData();
};

std::vector<std::vector<std::string> > CSVReader::getData()
{
	std::ifstream file("CSVv2_94XSF_V2_B_F.csv");

	std::vector<std::vector<std::string> > dataList;
	std::string line = "";
	// Iterate through each line and split the content using delimeter
	while (getline(file, line))
	{
		cout<<"file is opened in CSVreader"<<endl;
		std::vector<std::string> vec;
		boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
		dataList.push_back(vec);
	}
	// Close the File
	file.close();

	return dataList;
}
*/
/*class CSVRow
{
    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream, cell, ','))
            {
                m_data.push_back(cell);
            }
            // This checks for a trailing comma with no data after it.
            if (!lineStream && cell.empty())
            {
                // If there was a trailing comma then add an empty element.
                m_data.push_back("");
            }
        }
    private:
        std::vector<std::string>    m_data;
};
std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
} 

class CSVIterator
{
    public:
        typedef std::input_iterator_tag     iterator_category;
        typedef CSVRow                      value_type;
        typedef std::size_t                 difference_type;
        typedef CSVRow*                     pointer;
        typedef CSVRow&                     reference;

        CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) { ++(*this); }
        CSVIterator()                   :m_str(NULL) {}

        // Pre Increment
        CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
        // Post increment
        CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
        CSVRow const& operator*()   const       {return m_row;}
        CSVRow const* operator->()  const       {return &m_row;}

        bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
        bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
    private:
        std::istream*       m_str;
        CSVRow              m_row;
};

*/
namespace
{

constexpr double MAX_ELE_NUM{1};
constexpr double MIN_ELE_PT{45}; //{15}//min 12, AP 45, 
constexpr float MIN_ELE_LEADING_PT{35.f};
constexpr double MAX_ELE_ETA{2.5};
constexpr double ENDCAP_MIN_ETA{1.566};
constexpr double BARREL_MAX_ETA{1.4442};

constexpr double MAX_MU_NUM{1};
constexpr double MIN_MU_PT{40};//min 33,40 AP
constexpr float MIN_MU_LEADING_PT{26.f};
constexpr double MAX_MU_ETA{2.4};
constexpr float MU_LOOSE_ISO{0.15f};
constexpr float MU_TIGHT_ISO{0.25f};

constexpr double MIN_MET_PT{40};

constexpr float Z_MASS{91.1876f};
constexpr float Z_MASS_CUT{20.f};

constexpr float MAX_JET_ETA{4.7f};
constexpr float MIN_JET_PT{30.f};//{30.f};
constexpr float JET_ISO{0.4f};
constexpr unsigned MIN_JETS{4};
constexpr unsigned MAX_JETS{6};

constexpr float MAX_BJET_ETA{2.4f};
constexpr float MIN_BTAG_DISC{0.8838f};
constexpr unsigned MIN_BJETS{1};
constexpr unsigned MAX_BJETS{3};

constexpr float W_MASS{80.385f};
constexpr float W_MASS_CUT{20.f};

constexpr float TOP_MASS{172.5f};
constexpr float TOP_MASS_CUT{20.f};

constexpr float DELTA_R_ZL{1.6f};
constexpr float DELTA_PHI_ZW{2};
constexpr float DELTA_PHI_ZMET{2};

constexpr double PI{3.14};
constexpr double TZQ_W{0.0128};
constexpr double WWLNQQ_W{2.1740};
constexpr double WZLNQQ_W{0.2335};
constexpr double TTZQQ_W{0.0237};
constexpr double ZZLLQQ_W{0.0485};



enum class channels
{
    enu,
    mnu
};

[[gnu::const]] auto delta_phi(const float phi1, const float phi2)
{
    return vdt::fast_atan2f(vdt::fast_sinf(phi1 - phi2), vdt::fast_cosf(phi1 - phi2));
}

[[gnu::const]] auto deltaR(const float eta1, const float phi1, const float eta2, const float phi2)
{
    return std::sqrt(std::pow(eta1 - eta2, 2) + std::pow(delta_phi(phi1, phi2), 2));
}

template<typename T, typename U>
[[gnu::const]] bool all_equal(const T& t, const U& u)
{
    return t == u;
}

template<typename T, typename U, typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
{
    return t == u && all_equal(u, args...);
}

[[gnu::const]] auto inv_mass(const floats& pts, const floats& etas, const floats& phis, const floats& ms)
{
    if (!all_equal(pts.size(), etas.size(), phis.size(), ms.size()))
    {
        throw std::logic_error("Collections must be the same size");
    }
    else if (pts.empty())
    {
        throw std::logic_error("Collections must not be empty");
    }

    TLorentzVector vec{};
    for (size_t i{0}; i < pts.size(); i++)
    {
        TLorentzVector p{};
        p.SetPtEtaPhiM(pts[i], etas[i], phis[i], ms[i]);
        vec += p;
    }
	return boost::numeric_cast<float>(vec.M());
}

template<typename T>
[[gnu::const]] T select(const T& a, const ints& mask)
{
    return a[mask];
}
} // namespace

void analyse(int argc, char* argv[])
{
       //ROOT::EnableImplicitMT();

        std::cout << "I am starting"<<std::endl;
	setTDRStyle();

   	ROOT::RDataFrame dc{"Events", "/data/disk3/nanoAOD_2017/tZqlvqq/*.root"};

	//TChain
	/*TChain MCBG("Events");
	MCBG.Add("/data/disk0/nanoAOD_2017/WWToLNuQQ/*.root");
	MCBG.Add("/data/disk0/nanoAOD_2017/WZTo1L1Nu2Q/*.root");
	MCBG.Add("/data/disk0/nanoAOD_2017/ttZToQQ/*.root");
	MCBG.Add("/data/disk0/nanoAOD_2017/ZZTo2L2Q/*.root");

	ROOT::RDataFrame bgc(MCBG);*/

	ROOT::RDataFrame wwc{"Events", "/data/disk0/nanoAOD_2017/WWToLNuQQ/*.root"};
	ROOT::RDataFrame wzc{"Events", "/data/disk0/nanoAOD_2017/WZTo1L1Nu2Q/*.root"};
	ROOT::RDataFrame ttZc{"Events", "/data/disk0/nanoAOD_2017/ttZToQQ/*.root"};
	ROOT::RDataFrame zzc{"Events", "/data/disk0/nanoAOD_2017/ZZTo2L2Q/*.root"};

//need to add the chain for real data...

	//ROOT::RDataFrame Se{"Events","/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_Run*/*.root"};
	//ROOT::RDataFrame Sm{"Events","/data/disk3/nanoAOD_2017/SingleMuon_*/*.root"};
	//ROOT::RDataFrame met{"Events","/data/disk0/nanoAOD_2017/MET*/*.root"};

	auto d = dc.Range(0, 100);
	//auto bg = bgc.Range(0, 10000);

	auto ww = wwc.Range(0, 10000);
	auto wz = wzc.Range(0, 10000);
	auto ttZ = ttZc.Range(0, 10000);
	auto zz = zzc.Range(0, 10000);
	//auto se = Sec.Range(0, 10000);
	//auto sm = Smc.Range(0, 10000);
	//auto met = metc.Range(0, 10000);

	std::cout << "I have looked up the dataset"<<std::endl;
/////////////////////////////////////////////////////////////////////////// Number of Particles Per Event /////////////////////////////////////////////////////////////////////
/*	std::cout << " gonna do particle statistics"<<std::endl;
	auto countZ = [](const ints ids) -> int {return std::count(ids.begin(),ids.end(),23);};//  writing a function with constant integer input as id and function called countZ and i am saying read from start of ids to ends of ids and find all the ids with number 23 and count them all
        auto Znum = d.Define("numZ",countZ,{"GenPart_pdgId"}) // making pointer Znum which its column is filled with numZ that is the GenPArt_pdgId and it uses the countZ function
			.Filter([](int numZ){return numZ >= 1;},{"numZ"}) //setting filters for the numZ to exist in my events
                        .Count(); // count them all
//////////////////////////////////////////////////////////////// Mean Number of Particles //////////////////////////////////////////////////////////////////////////
	 std::cout << "a quick mean calculation"<<std::endl;
         auto bnumMean = d.Define("numb",countb,{"GenPart_pdgId"})
                         .Filter([](int numb){return numb >= 1;},{"numb"})
                         .Mean("numb"); // in here instead of counting sum of all events with numb (number of b per event), I am taking a mean of how many bs r per event

////////////////////////////////////////////////////////////////// Histograms of cross sections  //////////////////////////////////////////////////////////////////////////////////////
	 std::cout << "historgram plotting"<<std::endl;
	 auto bnumHist = d.Define("numb",countb,{"GenPart_pdgId"})
                         .Filter([](int numb){return numb >= 1;},{"numb"})
                         .Histo1D({"numb","Number of b quarks per event",2,0,10},"numb"); // This is how a histogtam is done, saying show me the dist. of numb
	 auto bHist = new TCanvas("Number of b quarks per event", "Number of b quarks per event",10,10,700,700);
	 bnumHist->DrawClone();
	 bHist->SaveAs("bnumHist.root"); // this is how u make a canvas out of it for bnumHist defined earlier :D

/////////////////////////////////////////////////////////////////////////// Pt Histograms ///////////////////////////////////////////////////////
 	std::cout << "pt histogram"<<std::endl;
 	auto bptfunc = [](const ints id, const floats pts) {return pts[id ==5];};
 	auto bptHist = d.Define("bpts",bptfunc,{"GenPart_pdgId","GenPart_pt"})
                        .Filter([](floats bpts){return std::all_of(bpts.cbegin(), bpts.cend(), [](float pt){return pt >= 20.;});},{"bpts"})
                         .Histo1D({"bpt","Pt of b quarks per event",15,0,150},"bpts"); // This is how a histogtam is done, sayin$
         auto bptcanvas = new TCanvas("Pt of b quarks per event", "Pt of b quarks per event",10,10,700,700);
         bptHist->SetLineColor(kBlue);
 	bptHist->DrawClone();
         bptcanvas->SaveAs("bptHist.root");


// /////////////////////////////////////////////////////////////////// Eta Histograms //////////////////////////////////////////////////////////////////////////////////
 	std::cout << "eta histogram"<<std::endl;
         auto betafunc = [](const ints id, const floats eta) {return eta[id ==5];};
         auto betaHist = d.Define("betas",betafunc,{"GenPart_pdgId","GenPart_eta"})
                        .Filter([](floats betas){return std::all_of(betas.cbegin(), betas.cend(), [](float eta){return abs(eta) <= 2.5;});},{"betas"})
                         .Histo1D({"beta","Eta of b quarks per event",20,-5,5},"betas"); // This is how a histogtam is done, sayin$
         auto betacanvas = new TCanvas("Eta of b quarks per event", "Eta of b quarks per event",10,10,700,700);
         betaHist->SetLineColor(kBlue);
         betaHist->DrawClone();
         betacanvas->SaveAs("betaHist.root");


///////////////////////////////////////////////////////////////// Phi Histogram ////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Phi Histogram"<<std::endl;
	auto dphifunc = [](const ints &id, const floats &pt, const floats &eta, const floats &phi)
	{
		ROOT::VecOps::RVec<int> index(id.size()); std::iota(index.begin(), index.end(),0);
		double phi1;
		for (const auto i : index[id == 1])
                {
                        if(pt.at(i)>= 30 && abs(eta.at(i))<=2.4)
                        {
				phi1 = phi.at(i);
                        }
                }

		return phi1;
	};



	auto dphihist = d.Define("dphis",dphifunc, {"GenPart_pdgId","GenPart_pt","GenPart_eta","GenPart_phi"})
			 .Histo1D({"dphi","Phi dist. of d quarks",50,-5,5},"dphis");

	auto dphicanvas = new TCanvas("d phi dist.","d phi dist.",10,10,700,700);

	dphicanvas ->cd();
	dphihist->SetLineColor(kRed);
	dphihist->DrawClone();

*/
////////////////////////////////////////////////////////////////////////Working Histo 2D Z Reconsturcion/////////////////////////////////////////////////////////////////////////////
/*        auto dZReconMassCut =[](const ints &id, const floats &pt, const floats &eta, const floats &phi, const floats &mass)
        {
          	ROOT::VecOps::RVec<int> index(id.size()); std::iota(index.begin(), index.end(),0); //making a RVec to find the index i need for genPArt_Pd
		double ZmassDs;
		if(index[id == 1].size() > 0 && index[id == -1].size() > 0)
		{
	             	for (const auto i : index[id ==1])
                	{
				if(index[id == 5].size()>0 || index[id == -5].size()>0 && pt.at(i)>= 70 && abs(eta.at(i))<=2.4)
				{
                       			ZmassD.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), mass.at(i));// zmassd will fill only for q
                		}
				for (const auto y : index[id == -1])
                		{
					if(index[id == 5].size()>0 || index[id == -5].size()>0 && pt.at(y)>= 70 && abs(eta.at(y))<=2.4)
					{
                       				ZmassDbar.SetPtEtaPhiM(pt.at(y), eta.at(y), phi.at(y), mass.at(y));//Zmassdbar will fill only for q bar
						ZmassDs = (ZmassD+ZmassDbar).M();
					}
					if(ZmassDs>=76 && ZmassDs <= 106)
					{
						index1= i;
						index2= y;
						return ZmassDs;
					}
					else
					{
						ZmassDs = 0;
					}
				}
			}
			if(ZmassDs == 0)
			{
				return 0.00;
			}
		}
		else
		{
			return 0.00;
		}
        }; // function calculating  reconstructed mass of q and q bar
       std::cout << "Z Recon is done "<<std::endl;
       std::cout << "Z Histogram"<<std::endl;
       auto dZReconMassPhiCut =[](const ints &id, const floats &pt, const floats &eta, const floats &phi, const floats &mass)
       {
                ROOT::VecOps::RVec<int> index(id.size()); std::iota(index.begin(), index.end(),0); //making a RVec to find the index i need for genPArt_Pd
		double phi1;
                double phi2;
		double deltaphi;
                if(index[id == 1].size()> 0 && index[id == -1].size()> 0)
		{
			for(const auto i: index[id==1])
                	{
				if(i == index1)
				{
                        			phi1 = phi.at(i);
				}
               			for(const auto y:index[id==-1])
                		{
					if(y == index2)
					{
	                       			phi2 = phi.at(y);
					}
				}
				deltaphi = phi1 - phi2;
				if(deltaphi > pi)
				{
					deltaphi = deltaphi - 2*pi;
				}
				else if(deltaphi < -pi)
				{
					deltaphi = deltaphi + 2*pi;
				}
				else
				{
					deltaphi = deltaphi;
				}
			}
			return abs(deltaphi);
         	}
		else
		{
		return 0.00;
		}
	};
       auto ZReconD2DHist = d.Define("dZReconMass",dZReconMassCut,{"GenPart_pdgId","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass"})
                           .Filter([](double dZReconMass){return dZReconMass > 0;},{"dZReconMass"})
	                   .Define("dZReconPhi",dZReconMassPhiCut,{"GenPart_pdgId","GenPart_pt","GenPart_eta","GenPart_phi","GenPart_mass"})
			   .Filter([](double dZReconPhi){return dZReconPhi > 2.26;},{"dZReconPhi"})
       			   .Histo2D({"ZReconDMassPhiHist","Recon. Z mass of d quarks Vs. delta phi",50,0,10,20,0,200},"dZReconPhi","dZReconMass");
       	auto dZReconMPhicanvas = new TCanvas("Recon Z  mass of d quarks per event", "Recon Z mass of d quarks per event",10,10,700,700);
       	ZReconD2DHist->DrawClone("colz");
       	dZReconMPhicanvas->SaveAs("ZReconDMassPhiHist_bPtEtaPhiCut.root");
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	auto is_good_ele{[](const int target_id, const bools& isPFs, const floats& pts, const floats& etas, const ints& ids){
		const auto abs_etas{abs(etas)};
        	return (isPFs && pts > MIN_ELE_PT
                && ((abs_etas < MAX_ELE_ETA && abs_etas > ENDCAP_MIN_ETA) || (abs_etas < BARREL_MAX_ETA))
                && ids >= target_id);
		}};

	auto is_good_tight_ele{[&is_good_ele](const bools& isPFs, const floats& pts, const floats& etas, const ints& ids) {
        	return is_good_ele(4, isPFs, pts, etas, ids);
	}};

    	auto is_good_loose_ele{[&is_good_ele](bools& isPFs, floats& pts, floats& etas, ints& ids) {
        	return is_good_ele(1, isPFs, pts, etas, ids);
   	}};

  	auto is_good_mu{[](const float target_iso, const bools& isPFs, const floats& pts, const floats& etas, const bools& ids, const floats& isos) {
        	auto abs_etas{abs(etas)};
        	return (isPFs && pts > MIN_MU_PT && abs_etas < MAX_MU_ETA && ids && isos <= target_iso);
   	}};
	auto is_good_tight_mu{[is_good_mu](const bools& isPFs, const floats& pts, const floats& etas, const bools& ids, const floats& isos) {
        	return is_good_mu(MU_TIGHT_ISO, isPFs, pts, etas, ids, isos);
    	}};
    	auto is_good_loose_mu{[is_good_mu](const bools& isPFs, const floats& pts, const floats& etas, const bools& ids, const floats& isos) {
        	return is_good_mu(MU_LOOSE_ISO, isPFs, pts, etas, ids, isos);
    	}};


	auto e_cut{[](const floats& tight_ele_pts, const floats& loose_ele_pts) {
        	const bool ele_cut{tight_ele_pts.size() == 1 && tight_ele_pts.size() == loose_ele_pts.size()};
        	//bool lead_pt_cut{false};
        	//lead_pt_cut = tight_ele_pts.empty() ? false : *std::max_element(tight_ele_pts.begin(), tight_ele_pts.end()) > MIN_ELE_PT;
		return ele_cut;
        	//return lead_pt_cut && ele_cut;
    	}};

	auto mu_cut{[](const floats& tight_mu_pts, const floats& loose_mu_pts) {
        	const bool mu_cut{tight_mu_pts.size() == 1 && tight_mu_pts.size() == loose_mu_pts.size()};
        	//bool lead_pt_cut{false};
        	//lead_pt_cut = tight_mu_pts.empty() ? false : *std::max_element(tight_mu_pts.begin(), tight_mu_pts.end()) > MIN_MU_PT;
        	//return lead_pt_cut && mu_cut;
		return mu_cut;
        }};

	auto get_w_e_quantity_selector{[](const std::string& s){
                return "Electron_" + s + "[tight_eles]";
        }};

        auto get_w_mu_quantity_selector{[](const std::string& s){
        	return "Muon_" + s + "[tight_mus]";
        }};

	auto ele_met_selection_function{[](const float& MET_electron_pt_Selection) ->bool{
		return MET_electron_pt_Selection > 80;
	}};

        auto mu_met_selection_function{[](const float& MET_muon_pt_Selection) ->bool { 
        	return MET_muon_pt_Selection > 40;
        }};

	auto transvers_W_mass{[](floats lep_pt,floats lep_phi,float met_pt,float met_phi){
  		//float w_reco_mass{std::numeric_limits<float>::infinity()};
                //size_t l_index_1{std::numeric_limits<size_t>::max()};
		floats w_mass_vec;
		for(int i{0}; i< lep_pt.size();i++)
		{
			const float  reco_mass = sqrt( 2 * lep_pt.at(i) * met_pt * (1 - cos(delta_phi(lep_phi.at(i), met_phi))) );
			//if (std::abs(W_MASS - reco_mass) < std::abs(W_MASS - w_reco_mass))
                        //{
                        	//w_reco_mass = reco_mass;
                                //l_index_1 = i;
                        //}
			w_mass_vec.push_back(reco_mass);
		}
		return w_mass_vec;
	}};


	/*auto w_mass_cut{[](const floats& w_mass) {;
        	return std::abs(w_mass - W_MASS) < W_MASS_CUT;
   	}};*/


	auto jet_lep_min_deltaR{[](const floats& jet_etas, const floats& jet_phis, const floats& lep_etas, const floats& lep_phis) {
        	floats min_dRs{};
        	std::transform(jet_etas.begin(), jet_etas.end(), jet_phis.begin(), std::back_inserter(min_dRs), [&](float jet_eta, float jet_phi) { return deltaR(jet_eta, jet_phi, lep_etas.at(0), lep_phis.at(0));});
		return min_dRs;
    	}};


	auto tight_jet_id{[](const floats& jet_lep_min_dRs, const floats& pts, const floats& etas, const ints& ids) {
		return (pts > MIN_JET_PT && etas < MAX_JET_ETA && jet_lep_min_dRs > JET_ISO && ids >= 2);
	}};

        auto jet_cut{[](const ints& tight_jets) {
        	auto njet{std::count_if(tight_jets.begin(), tight_jets.end(), [](int i) { return i; })};
        	return (njet >= MIN_JETS) && (njet <= MAX_JETS);
	}};


	auto bjet_id{[](const ints& tight_jets, const floats& btags, const floats& etas){
		return  tight_jets && (btags > 0.8838f) && (abs(etas) < 2.4f);
	}};

	auto bjet_cut{[](const ints& bjets) {const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        	return nbjet >= 1 && nbjet <= 3;
	}};

	auto bjet_variable{[](const floats& Jet_variable, const unsigned int& nJet, const ints& lead_bjet){
		floats vec{};
		for(int i = 0; i < nJet; i++)
		{
 			if(lead_bjet.at(i) == 1)
			{
				vec.push_back(Jet_variable.at(i));
			}
		}
		return vec;
	}};

	auto numberofbjets{[](const ints& bjets) {const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        	return nbjet;
	}};

	auto find_lead_mask{[](const ints& mask, const floats& vals){
		const auto masked_vals{mask * vals};
		const auto max_idx{boost::numeric_cast<size_t>(std::distance(masked_vals.begin(), max_element(masked_vals.begin(), masked_vals.end())))};
		ints lead_mask(masked_vals.size(), 0);
		lead_mask.at(max_idx) = 1;
		return lead_mask;
	}};

	auto jet_deltaphi_func{[](const floats& phis){
		float jet_deltaphi;
		floats deltaphis;
		for(int i{0};i < phis.size(); i++)
		{
			for(int j{i+1}; j < phis.size(); j++)
			{

                		jet_deltaphi = abs(delta_phi(phis.at(i),phis.at(j)));
			}
			deltaphis.push_back(jet_deltaphi);
		}
		return deltaphis;
	}};


	auto find_z_pair{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet){
		double z_reco_mass{std::numeric_limits<double>::infinity()};
		size_t jet_index_1{std::numeric_limits<size_t>::max()};
		size_t jet_index_2{std::numeric_limits<size_t>::max()};
		const size_t njets{pts.size()};
		for (size_t i{0}; i < njets; ++i)
		{
			for (size_t j{i + 1}; j < njets; ++j)
 			{
				if(tight_jets[i] != 0 && tight_jets[j] != 0 && lead_bjet[i] != 1 && lead_bjet[j] != 1)
                		{
                    			continue;
	               		}
                                auto jet1{TLorentzVector{}};
                		auto jet2{TLorentzVector{}};
                		jet1.SetPtEtaPhiM(pts.at(i), etas.at(i), phis.at(i), ms.at(i));
                		jet2.SetPtEtaPhiM(pts.at(j), etas.at(j), phis.at(j), ms.at(j));

               			if (const double reco_mass{(jet1 + jet2).M()}; std::abs(Z_MASS - reco_mass) < std::abs(Z_MASS - z_reco_mass))
                		{
					z_reco_mass = reco_mass;
                    			jet_index_1 = i;
                    			jet_index_2 = j;
       				}
        		}
        	}

		ints z_pair(njets, 0);
        	z_pair.at(jet_index_1) = 1;
        	z_pair.at(jet_index_2) = 1;
   		return z_pair;
	}};

	auto z_mass_cut{[](const float& z_mass){
   		return abs(z_mass - Z_MASS) < Z_MASS_CUT;
	}};


	auto deltaR_z_l{[](const floats& deltaRzl){
		return std::any_of(deltaRzl.cbegin(), deltaRzl.cend(), [](float delta){return delta >= DELTA_R_ZL;});
        }};

	auto ZW_deltaphi_func{[](const floats& phis1, const floats& phis2){
		float phi1;
		float phi2;
		float deltaphi;
		floats deltaphi_vec;
                for(int i; i < phis1.size(); i++)
                {
                        phi1 = phis1.at(i);
                        for(int j; j < phis2.size(); j++)
                        {
				phi2 = phis2.at(j);
				deltaphi = abs(delta_phi(phi1,phi2));
				deltaphi_vec.push_back(deltaphi);
			}
		}
		return deltaphi_vec;
	}};

        auto ZW_deltaphi_cut{[](const floats deltaphi){
		return std::any_of(deltaphi.cbegin(), deltaphi.cend(), [](float delta){return delta >= DELTA_PHI_ZW;});
        }};

	auto ZMet_deltaphi_func{[](const floats& z_phi, const float met_pt){
		float deltaphi;
		floats deltaphi_vec;
			for(int i{0}; i < z_phi.size(); i++)
			{
				deltaphi = abs(delta_phi(z_phi.at(i), met_pt));
				deltaphi_vec.push_back(deltaphi);
			}
		return deltaphi_vec;
	}};


        auto ZMet_deltaphi_cut{[](const floats deltaPhi){
		return std::any_of(deltaPhi.cbegin(), deltaPhi.cend(), [](float delta){return delta >= DELTA_PHI_ZMET;});
        }};

	auto TLorentzVectorMass{[](const TLorentzVector& object){
		const float mass{object.M()};
		return mass;
	}};

	auto TLorentzVectorPt{[](const TLorentzVector& object){
		return object.Pt();
	}};

	auto TLorentzVectorPhi{[](const TLorentzVector& object){
		return object.Phi();
	}};

	auto TLorentzVectorEta{[](const TLorentzVector& object){
		return object.Eta();
	}};


	auto BLorentzVector{[](const floats& bjet_pt, const floats& bjet_eta, const floats& bjet_phi, const floats& bjet_mass, const long& nbjets){
		auto BJets = TLorentzVector{};
		for(int i = 0; i < bjet_pt.size(); i++)
		{
			BJets.SetPtEtaPhiM(bjet_pt.at(i), bjet_eta.at(i), bjet_phi.at(i), bjet_mass.at(i));
		}
		return BJets;
	}};

	auto top_reconstruction_function{[](const floats& bjets_pt, const floats& bjets_eta, const floats& bjets_phi, const floats& bjets_mass,
		const floats& w_pair_pt, const floats& w_pair_eta, const floats& w_pair_phi, const floats& w_mass){

		float t_reco_mass{std::numeric_limits<float>::infinity()};
		const size_t nbjets{bjets_pt.size()};
		const size_t nWs{w_pair_pt.size()};
		size_t bjet_index{std::numeric_limits<size_t>::max()};
                size_t W_index{std::numeric_limits<size_t>::max()};

		auto BJets{TLorentzVector{}};
                auto RecoW{TLorentzVector{}};
                auto reco_top{TLorentzVector{}};

		for(int i{0}; i < nbjets; i++)
		{
			for(int j{0}; j < nWs; j++)
			{
                               	BJets.SetPtEtaPhiM(bjets_pt.at(i), bjets_eta.at(i), bjets_phi.at(i), bjets_mass.at(i));
                               	RecoW.SetPtEtaPhiM(w_pair_pt.at(j), w_pair_eta.at(j), w_pair_phi.at(j), w_mass.at(j));
				if(abs(RecoW.M() - W_MASS) < W_MASS_CUT)
				{
		 			if (const float reco_mass{(RecoW + BJets).M()}; std::abs(TOP_MASS - reco_mass) < std::abs(TOP_MASS - t_reco_mass))
                                        {
						reco_top = RecoW + BJets;
						//t_reco_mass = reco_mass;
					}
				}
			}
		}
		return reco_top;
	}};

        auto top_mass_cut{[](const float& top_mass){
                return abs(top_mass - TOP_MASS) < TOP_MASS_CUT;
        }};


	vector<string> bjet_mass_strings = {"Jet_mass", "nJet", "lead_bjet"};
	vector<string> bjet_eta_strings = {"Jet_eta", "nJet", "lead_bjet"};
	vector<string> bjet_pt_strings = {"Jet_pt", "nJet", "lead_bjet"};
	vector<string> bjet_phi_strings = {"Jet_phi", "nJet", "lead_bjet"};

// Try to open the file 
	auto open_CSV_file{[](const floats& dummy){
		int value;
		ifstream ip;
		ip.open("/home/eepgssg/Shirin-Codes/CSVv2_94XSF_V2_B_F.csv", ios::in);
		if(ip.is_open())
		{
			cout<<"CSV is open"<<endl;
			value = 0;
		}
		else
		{
			cout<<"CSV can not open"<<endl;
			value = 1;
		}
		ip.close();
		return value;
		/*ofstream tf;
		tf.open("test.txt");
		if(tf.is_open())
		{
			tf << "open\n";
			tf.close();
		}
		return 1;
		*/
/*		cout<<"i am in open csv"<<endl;
		CSVReader reader("/home/eepgssg/Shirin-Codes/CSVv2_94XSF_V2_B_F.csv");
		cout<<"passed giving the address"<<endl;
		// Get the data from CSV File
		vector<vector<string> > dataList = reader.getData();
		cout<<"got passed reader.getData()"<<endl;
		// Print the content of row by row on screen
		for(std::vector<std::string> vec : dataList)
		{
			cout<<"i am inside the first loop of reading file"<<endl;
			for(string data : vec)
			{
				cout<<"I am reading the file"<<endl;
				cout<<data << " , ";
			}
		}
		return 0;*/
		/*ifstream file("CSVv2_94XSF_V2_B_F.csv");
		cout<<"file is declared"<<endl;
    		CSVRow row;
    		while(file >> row)
    		{
        		std::cout << "4th Element(" << row[3] << endl;;
    		}
		*/
/*		ifstream file("/home/eepgssg/Shirin-Codes/CSVv2_94XSF_V2_B_F.csv");
		cout<<"file is declared"<<endl;
    		for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
    		{
			cout<< "file is iterating"<<endl;
        		//std::cout << "4th Element(" << (*loop)[3] << ")\n";
    		}
		return 1;
*/
		}};

// B Tag Efficiency centre
	auto btag_CSVv2_formula{[](const floats& btag, const floats& pt, const floats& eta){
		strings formula;
		string x ("x");
		floats result;

	  	ifstream ip("/home/eepgssg/Shirin-Codes/CSVv2_94XSF_V2_B_F.csv");
		//ip.open("CSVv2_94XSF_V2_B_F.csv");
  		if(!ip.is_open()) std::cout << "ERROR: File Open" << '\n';
		if(ip.is_open()) cout<< "CSVV IS OPEN"<< endl;
        	string CSVv2;
        	string measure_type;
        	string sys_type;
        	string jet_flav;
        	string eta_min;
        	string eta_max;
        	string pt_min;
        	string pt_max;
        	string CSV_min;
        	string CSV_max;
        	string formular;

	  	while(ip.good())
		{
			cout<<"File is read"<<endl;
               		getline(ip, CSVv2, ',');
               		getline(ip, measure_type, ',');
               		getline(ip, sys_type, ',');
               		getline(ip, jet_flav, ',');
               		getline(ip, eta_min, ',');
               		getline(ip, eta_max, ',');
               		getline(ip, pt_min, ',');
               		getline(ip, pt_max, ',');
               		getline(ip, CSV_min, ',');
               		getline(ip, CSV_max, ',');
               		getline(ip, formular, '\n');

			cout<< "This is CSVvs operator "<<CSVv2<<endl;
			cout<< "This is measure_type "<<measure_type<<endl;
	      		cout<< "This is sys_type "<<sys_type<<endl;
			cout<< "This is jet_flav "<<jet_flav<<endl;
			cout<< "eta min "<<eta_min<<endl;
			cout<< "eta max "<<eta_max<<endl;
			cout<< "This is pt min "<<pt_min<<endl;
			cout<< "This is pt max "<<pt_max<<endl;
			cout<< "This is CSV min "<<CSV_min<<endl;
			cout<< "This is CSV max "<<CSV_max<<endl;
			cout<< "this is formular "<<formular<<endl;
			cout<<"size of pt is "<<pt.size()<<" pt is "<<pt.at(0)<<" pt_min is "<<pt_min<<" pt max is "<<pt_max<<endl;
			cout<<"before loop"<<endl;
			for(int i{0}; i < pt.size(); i++)
                        {
				float CSVv2_v = boost::lexical_cast<float>(CSVv2); cout<<"CSV value is "<<CSVv2_v<<endl;
				float jet_flav_v = boost::lexical_cast<float>(jet_flav); cout<<"jet flav value is "<<jet_flav_v<<endl;
				float eta_min_v = boost::lexical_cast<float>(eta_min); cout<<"eta min value is "<<eta_min_v<<endl;
				float eta_max_v = boost::lexical_cast<float>(eta_max); cout<<"eta max value is "<<eta_max_v<<endl;
				float pt_min_v = boost::lexical_cast<float>(pt_min); cout<<"pt min value is "<<pt_min_v<<endl;
				float pt_max_v = boost::lexical_cast<float>(pt_max); cout<<"pt max value is "<<pt_max_v<<endl;
				float CSV_min_v = boost::lexical_cast<float>(CSV_min); cout<<"CSV min value is "<<CSV_min_v<<endl;
				float CSV_max_v = boost::lexical_cast<float>(CSV_max); cout<<"CSV max value is "<<CSV_max_v<<endl;

                        	cout << "CSV works"<<endl;
				cout<<"pt is "<<pt.at(i)<<" pt_min is "<<pt_min<<"pt max is "<<pt_max<<endl;
                                if(CSVv2_v >= 0.8838 && measure_type == "comb" && jet_flav_v  == 0 && eta.at(i) > eta_min_v && eta.at(i) < eta_max_v && pt.at(i) > pt_min_v && pt.at(i) < pt_max_v && btag.at(i) > CSV_min_v && btag.at(i) < CSV_max_v)
                                {
                                	cout <<formular<<"  This is the formula"<<endl;
                                        formula.push_back(formular);
                                }
                                else
                                {
					cout<<"combination for the formula was not found"<<endl;
                                	formula.push_back("0");
                                }
			}
			cout<<"after loop"<<endl;
		}
		ip.close();
		cout<<"size of formula is "<<formula.size()<<endl;
	      	for(int j; j< formula.size(); j++)
	        {
			cout<<"i am starting to parse the formula   "<<"formula is"<<formula.at(j)<<endl;
			// convert bpt to string
			string bpt_string = boost::lexical_cast<string>(pt.at(j));
			//replace all x with bpts values and "" with space
			string form;
			form = formula.at(j);
			form.replace(form.begin(), form.end(), 'x', 'bpt_string');
			form.replace(form.begin(), form.end(), '"', ' ');
			//use the parser
			Eval ev;
			complex<float> res;
			res = ev.eval((char*)form.c_str());
    			result.push_back(res.real());

		}
	        cout<< "returning evaluated result in btag"<<endl;
		return result;

	}};

	auto non_btag_CSVv2_formula{[](const floats& btag, const floats& pt, const floats& eta){
                strings formula;
                string x ("x");
                floats result;

                ifstream ip;//("CSVv2_94XSF_V2_B_F.csv");
		ip.open("/home/eepgssg/Shirin-Codes/CSVv2_94XSF_V2_B_F.csv");
                if(!ip.is_open()) std::cout << "ERROR: File Open" << '\n';
                string CSVv2;
                string measure_type;
                string sys_type;
                string jet_flav;
                string eta_min;
                string eta_max;
                string pt_min;
                string pt_max;
                string CSV_min;
                string CSV_max;
                string formular;

                while(ip.good())
                {
                        getline(ip, CSVv2, ',');
                        getline(ip, measure_type, ',');
                        getline(ip, sys_type, ',');
                        getline(ip, jet_flav, ',');
                        getline(ip, eta_min, ',');
                        getline(ip, eta_max, ',');
                        getline(ip, pt_min, ',');
                        getline(ip, pt_max, ',');
                        getline(ip, CSV_min, ',');
                        getline(ip, CSV_max, ',');
                        getline(ip, formular, '\n');

                        for(int i; i <pt.size(); i++)
                        {
                                cout << "CSV works"<<endl;
                                if(measure_type == "comb" && jet_flav  == "0" && eta.at(i) > eta_min && eta.at(i) < eta_max && pt.at(i) > pt_min && pt.at(i) < pt_max) //&& btag.at(i) > CSV_min && btag.at(i) < CSV_max)

                                {
                                        cout <<formular<<"  This is the formula"<<endl;
                                        formula.push_back(formular);
				}
                                else
                                {
                                        formula.push_back("0");
                                }
			}
                        for(int j; j< formula.size(); j++)
                        {
                               // convert bpt to string
                               string pt_string = boost::lexical_cast<string>(pt.at(j));
                               //replace all x with bpts values and "" with space
                               string form;
                               form = formula.at(j);
                               form.replace(form.begin(), form.end(), 'x', 'pt_string');
                               form.replace(form.begin(), form.end(), '"', ' ');
                               //use the parser
                               Eval ev;
                               complex<float> res;
                               res = ev.eval((char*)form.c_str());
                               result.push_back(res.real());

                        }
                        ip.close();

                }
		ip.close();
                cout<< "returning evaluated result"<<endl;
		cout<< " size of pt is "<<pt.size()<<endl;
                return result;

        }};


	// e calculation for b tagging efficiency (b tagged),ei
	// numerator

	auto bjet_id_numer{[](const ints& tight_jets, const floats& btags, const floats& etas,const ints& Gen_id){
		ints CountVec{};
		for(int i; i< btags.size(); i++) // this loop could be incorrectly defined , due to not selecting the jets from tight_jets but selecting the number of $
                {
			if(Gen_id.at(i) == 5 && btags.at(i) > 0.8838f && (abs(etas.at(i)) < 2.4f) && tight_jets.at(i) > 0)
			{
				CountVec.push_back(1.0);
			}
			else
			{

				CountVec.push_back(0.0);
			}
		}
		return  CountVec;

	}};
	// denominator
	auto bjet_id_denom{[](const ints& tight_jets, const floats& etas, const ints& Gen_id){
        	ints CountVec;
		for(int i; i< etas.size(); i++)
		{

			if(Gen_id.at(i) == 5 && (abs(etas.at(i)) < 2.4f) && tight_jets.at(i) >0)
                        {
                                CountVec.push_back(1.0);
                        }
                        else
                        {
                                CountVec.push_back(0.0);
                        }
                }
               return  CountVec;

        }};

        auto non_bjet_id_numer{[](const ints& tight_jets, const floats& btags, const floats& etas,const ints& Gen_id){
                ints CountVec;
                for(int i; i< etas.size(); i++)
                {
                        if(tight_jets.at(i) > 0 && btags.at(i) > 0.8838f && (abs(etas.at(i)) < 2.4f) && ((Gen_id.at(i) > 0 && Gen_id.at(i) <= 4) || Gen_id.at(i) == 21))
                        {
                                CountVec.push_back(1.0);
                        }
			else
			{
				CountVec.push_back(0.0);
			}
                }
		return CountVec;

        }};
        // denominator
   	auto non_bjet_id_denom{[](const ints& tight_jets, const floats& etas, const ints& Gen_id){
		ints CountVec;
                for(int i; i< etas.size(); i++)
                {
                        if(tight_jets.at(i)>0 && (abs(etas.at(i)) < 2.4f) && ((Gen_id.at(i) > 0 && Gen_id.at(i) <= 4) || Gen_id.at(i) == 21))
                        {
                                CountVec.push_back(1.0);
                        }
			else
			{
				CountVec.push_back(0.0);
			}
                }
		return CountVec;
        }};

	//Product formula for MC
	auto EffBTaggedProduct{[](const floats& EffBTagged){
		cout << "inside EffBTaggedProduct" << endl;
		float initial = 1;
		//floats Vec;
		for(int i = 0; i < EffBTagged.size(); i++)
		{
			initial = EffBTagged.at(i) * initial;
			//Vec.push_back(initial);
		}
		//cout<<"in effbtagged size of vec is "<<Vec.size()<<endl;
		return initial;
	}};
	auto EffNonBTaggedProduct{[](const floats& EffNonBTagged){
		float initial = 1;
		//floats Vec;
		for(int i = 0; i < EffNonBTagged.size(); i++ )
		{
			initial = (1 - EffNonBTagged.at(i)) * initial;
			//Vec.push_back(initial);
		}
		//cout<<"in effnonbtagged size of vec is "<<Vec.size()<<endl;
		return initial;
	}};
	auto P_MC_func{[](const float pi_ei, const float pi_ej){
		cout<<"I am inside P(MC)"<<endl;
		/*floats Vec;
		int size = (pi_ei.size() < pi_ej.size()) ? pi_ei.size() : pi_ej.size();
		for(int i = 0; i< size; i++)
		{
			Vec.push_back(pi_ei.at(i) * pi_ej.at(i));
		}
		return Vec;*/
		return pi_ei * pi_ej;
	}};
        auto Sfi_EffBTaggedProduct{[](const floats& EffBTagged, const floats sfi){
                cout << "inside Si_EffBTaggedProduct" << endl;
                float initial = 1;
                int size = (EffBTagged.size() < sfi.size()) ? EffBTagged.size() : sfi.size();
                for(int i = 0; i < size; i++)
                {
                        initial = sfi.at(i) * EffBTagged.at(i) * initial;
                }
                //cout<<"in effbtagged size of vec is "<<Vec.size()<<endl;
                return initial;
        }};

        auto Sfj_EffNonBTaggedProduct{[](const floats& EffNonBTagged, const floats sfj){
                cout << "inside Si_EffNonBTaggedProduct" << endl;
                float initial = 1;
                int size = (EffNonBTagged.size() < sfj.size()) ? EffNonBTagged.size() : sfj.size();
                for(int i = 0; i < size; i++)
                {
                        initial = (1 - EffNonBTagged.at(i) * sfj.at(i)) * initial;
                }
                //cout<<"in effbtagged size of vec is "<<Vec.size()<<endl;
                return initial;
        }};
	auto P_data_func{[](const float pi_ei, const float pi_ej){
                cout<<"I am inside P(Data)"<<endl;
                /*floats Vec;
                int size = (pi_ei.size() < pi_ej.size()) ? pi_ei.size() : pi_ej.size();
                for(int i = 0; i< size; i++)
                {
                        Vec.push_back(pi_ei.at(i) * pi_ej.at(i));
                }
                return Vec;*/
                return pi_ei * pi_ej;
        }};


// Variable Luminosity scale factor for all
	auto VarFact_func_double{[](const double& i){// this function make the variable which is equal to one and will be used for all scale factors
		return 1.0f;
	}};
	auto VarFact_func{[](const float& i){// this function make the variable which is equal to one and will be used for all scale factors
                return 1.0f;
        }};

//Signal Luminosity normalization
	auto NormScaleFact_func_double{[&VarFact_func_double](const doubles& i){// this function calculates the weight scale factor
		doubles weight_vec;
		for(int w; w < i.size(); w++)
		{
			double n_w;
			n_w = VarFact_func_double(i.at(w))*TZQ_W;
		        weight_vec.push_back(n_w);
		}
		return weight_vec;
	}};

        auto NormScaleFact_func{[&VarFact_func](const floats& i){// this function calculates the weight scale factor
                floats weight_vec;
                for(int w; w < i.size(); w++)
                {
			float n_w;
                        n_w = VarFact_func(i.at(w))*TZQ_W;
			weight_vec.push_back(n_w);
                }
                return weight_vec;
        }};

        auto NormScaleFact_func_novec{[&VarFact_func](const float& i){// this function calculates the weight scale factor for no Rvector variables
                float weight;
                weight = VarFact_func(i)*TZQ_W;
		return weight;
        }};
// WW scale factors

        auto WW_NormScaleFact_func_double{[&VarFact_func_double](const doubles& i){// this function calculates the weight scale factor
                doubles weight_vec;
                for(int w; w < i.size(); w++)
                {
                        double n_w;
                        n_w =  VarFact_func_double(i.at(w))*WWLNQQ_W;
                        weight_vec.push_back(n_w);
                }
                return weight_vec;
        }};

        auto WW_NormScaleFact_func{[&VarFact_func](const floats& i){// this function calculates the weight scale factor
                floats weight_vec;
                for(int w; w < i.size(); w++)
                {
                        float n_w;
                        n_w = VarFact_func(i.at(w))*WWLNQQ_W;
                        weight_vec.push_back(n_w);
                }
                return weight_vec;
        }};

        auto WW_NormScaleFact_func_novec{[&VarFact_func](const float& i){// this function calculates the weight scale factor for no Rvector variables
                float weight;
                weight = VarFact_func(i)*WWLNQQ_W;
                return weight;
        }};


//WZ
	auto WZ_NormScaleFact_func_double{[&VarFact_func_double](const doubles& i){// this function calculates the weight scale factor
                doubles weight_vec;
                for(int w; w < i.size(); w++)
                {
                        double n_w;
                        n_w = VarFact_func_double(i.at(w))*WZLNQQ_W;
                        weight_vec.push_back(n_w);
                }
                return weight_vec;
        }};

	auto WZ_NormScaleFact_func{[&VarFact_func](const floats& i){// this function calculates the weight scale factor
                floats weight_vec;
                for(int w; w < i.size(); w++)
                {
                        float n_w;
                        n_w = VarFact_func(i.at(w))*WZLNQQ_W;
                        weight_vec.push_back(n_w);
                }
                return weight_vec;
     	}};

	auto WZ_NormScaleFact_func_novec{[&VarFact_func](const float& i){// this function calculates the weight scale factor for no Rvector variables
                float weight;
                weight = VarFact_func(i)*WZLNQQ_W;
                return weight;
        }};


//ZZ
	auto ZZ_NormScaleFact_func_double{[&VarFact_func_double](const doubles& i){// this function calculates the weight scale factor
                doubles weight_vec;
                for(int w; w < i.size(); w++)
                {
                        double n_w;
                        n_w =  VarFact_func_double(i.at(w))*ZZLLQQ_W;
                        weight_vec.push_back(n_w);
                }
                return weight_vec;
        }};

        auto ZZ_NormScaleFact_func{[&VarFact_func](const floats& i){// this function calculates the weight scale factor
                floats weight_vec;
                for(int w; w < i.size(); w++)
                {
                        float n_w;
                        n_w = VarFact_func(i.at(w))*ZZLLQQ_W;
                        weight_vec.push_back(n_w);
                }
                return weight_vec;
        }};

        auto ZZ_NormScaleFact_func_novec{[&VarFact_func](const float& i){// this function calculates the weight scale factor for no Rvector variables
                float weight;
                weight = VarFact_func(i)*ZZLLQQ_W;
                return weight;
        }};

//TTZ
        auto TTZ_NormScaleFact_func_double{[&VarFact_func_double](const doubles& i){// this function calculates the weight scale factor
                doubles weight_vec;
                for(int w; w < i.size(); w++)
                {
                        double n_w;
                        n_w =  VarFact_func_double(i.at(w))*TTZQQ_W;
                        weight_vec.push_back(n_w);
                }
                return weight_vec;
        }};

	auto TTZ_NormScaleFact_func{[&VarFact_func](const floats& i){// this function calculates the weight scale factor
                floats weight_vec;
                for(int w; w < i.size(); w++)
                {
                        float n_w;
                        n_w = VarFact_func(i.at(w))*TTZQQ_W;
                        weight_vec.push_back(n_w);
                }
                return weight_vec;
        }};

        auto TTZ_NormScaleFact_func_novec{[&VarFact_func](const float& i){// this function calculates the weight scale factor for no Rvector variables
                float weight;
                weight = VarFact_func(i)*TTZQQ_W;
                return weight;
        }};
//////////////////////////////////////////////////////////////////////////////Signal ////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////Electron Channel/////////////////////////////////////////////////////////////////////////
	auto d_enu_event_selection = d.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("tight_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("tight_ele_eta", select<floats>, {"Electron_eta", "tight_eles"})
					.Define("tight_ele_phi", select<floats>, {"Electron_phi", "tight_eles"})
                   			.Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("loose_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_electron_pt_Selection",{"MET_pt"})
					.Filter(ele_met_selection_function, {"MET_electron_pt_Selection"}, "MET PT CUT")
					.Filter(e_cut, {"tight_ele_pt", "loose_ele_pt"}, "lepton cut");

	auto d_enu_w_selection = d_enu_event_selection.Define("w_e_eta", get_w_e_quantity_selector("eta"))
                 					.Define("w_e_phi", get_w_e_quantity_selector("phi"))
                 					.Define("w_e_pt", get_w_e_quantity_selector("pt"))
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"});
                 					//.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto d_enu_jets_selection = d_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
							.Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
							.Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
							.Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
							.Define("tight_jets_Genid", select<ints>, {"GenPart_pdgId", "tight_jets"})
							.Define("tight_jets_btagCSVV2", select<floats>, {"Jet_btagCSVV2", "tight_jets"})
							.Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                   					.Filter(jet_cut, {"tight_jets"}, "Jet cut");

	auto d_enu_jets_bjets_selection = d_enu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
								.Define("btag_numer", bjet_id_numer, {"tight_jets", "Jet_btagCSVV2", "Jet_eta", "GenPart_pdgId"})
								.Define("btag_denom", bjet_id_denom, {"tight_jets", "Jet_eta", "GenPart_pdgId"})
								.Define("btag_numer_pt", select<floats>, {"Jet_pt", "btag_numer"})
								.Define("btag_numer_eta",select<floats>, {"Jet_eta", "btag_numer"})
								.Define("btag_denom_pt", select<floats>, {"Jet_pt", "btag_denom"})
								.Define("btag_denom_eta",select<floats>, {"Jet_eta", "btag_denom"})
								.Define("non_btag_numer", non_bjet_id_numer,{"tight_jets", "Jet_btagCSVV2", "Jet_eta", "GenPart_pdgId"})
								.Define("non_btag_denom", non_bjet_id_denom, {"tight_jets", "Jet_eta", "GenPart_pdgId"})
                                                                .Define("non_btag_numer_pt", select<floats>, {"Jet_pt", "non_btag_numer"})
                                                                .Define("non_btag_numer_eta",select<floats>, {"Jet_eta", "non_btag_numer"})
                                                                .Define("non_btag_denom_pt", select<floats>, {"Jet_pt", "btag_denom"})
                                                                .Define("non_btag_denom_eta",select<floats>, {"Jet_eta", "btag_denom"})
								.Define("sfi", btag_CSVv2_formula, {"Jet_btagCSVV2","tight_jets_pt","tight_jets_eta"})
								.Define("sfj", non_btag_CSVv2_formula, {"Jet_btagCSVV2","tight_jets_pt","tight_jets_eta"});



	auto d_enu_z_rec_selection = d_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                                                                .Define("z_e_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_ele_eta", "tight_ele_phi"})
							        .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_e_phi"})
								.Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_electron_pt_Selection"});
								//.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
								//.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
								//.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                 						//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


	auto d_enu_brec_selection = d_enu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
							.Define("bjetpt", bjet_variable, bjet_pt_strings)
							.Define("bjeteta", bjet_variable, bjet_eta_strings)
							.Define("bjetphi", bjet_variable, bjet_phi_strings)
							.Define("nbjets", numberofbjets, {"bjets"})
							.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

	auto d_enu_top_selection = d_enu_brec_selection	.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_e_pt", "w_e_eta", "w_e_phi", "w_e_mass"})
							.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
							.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
							.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
							.Define("nw_tight_ele_pt", NormScaleFact_func, {"tight_ele_pt"})
							.Define("nw_tight_ele_eta", NormScaleFact_func, {"tight_ele_eta"})
							.Define("nw_tight_jets_eta", NormScaleFact_func, {"tight_jets_eta"})
							.Define("nw_tight_jets_pt", NormScaleFact_func, {"tight_jets_pt"})
							.Define("nw_jet_e_min_dR", NormScaleFact_func, {"jet_e_min_dR"})
							.Define("nw_z_e_min_dR", NormScaleFact_func, {"z_e_min_dR"})
							.Define("nw_w_e_mass", NormScaleFact_func, {"w_e_mass"})
							.Define("nw_z_mass", NormScaleFact_func_novec, {"z_mass"})
							.Define("nw_tight_jets_deltaphi", NormScaleFact_func, {"tight_jets_deltaphi"})
							.Define("nw_ZMet_deltaphi", NormScaleFact_func, {"ZMet_deltaphi"})
							.Define("nw_ZW_deltaphi", NormScaleFact_func, {"ZW_deltaphi"});

        auto h_d_enu_events_btag_numer_PtVsEta = d_enu_top_selection.Histo2D({"MC btag_Pt_vs_eta_enu_Channel","MC btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3},"btag_numer_pt", "btag_numer_eta");
        auto h_d_enu_events_non_btag_numer_PtVsEta = d_enu_top_selection.Histo2D({"MC non btag_Pt_vs_eta_enu_Channel","MC non btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3},"non_btag_numer_pt","non_btag_numer_eta");
        auto h_d_enu_events_btag_denom_PtVsEta = d_enu_top_selection.Histo2D({"MC btag_Pt_vs_eta_enu_Channel","MC btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3},"btag_denom_pt","btag_denom_eta");
        auto h_d_enu_events_non_btag_denom_PtVsEta = d_enu_top_selection.Histo2D({"MC non btag_Pt_vs_eta_enu_Channel","MC non btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3},"non_btag_denom_pt","non_btag_denom_eta");

//Write histogram to a root file:
	h_d_enu_events_btag_numer_PtVsEta->Write();
        h_d_enu_events_non_btag_numer_PtVsEta->Write();
        h_d_enu_events_btag_denom_PtVsEta->Write();
        h_d_enu_events_non_btag_denom_PtVsEta->Write();


	auto BTaggedBinFunction{[&h_d_enu_events_btag_numer_PtVsEta, &h_d_enu_events_btag_denom_PtVsEta](const floats& pts, const floats& etas){
		floats BTaggedEff{};
                for(int i = 0; i < pts.size(); i++)
                {
			int PtNum;
			int EtaNum;
			int PtDenom;
			int EtaDenom;
			PtNum = h_d_enu_events_btag_numer_PtVsEta->GetXaxis()->FindBin(pts.at(i));
			EtaNum = h_d_enu_events_btag_numer_PtVsEta->GetYaxis()->FindBin(etas.at(i));
			PtDenom = h_d_enu_events_btag_denom_PtVsEta->GetXaxis()->FindBin(pts.at(i));
			EtaDenom = h_d_enu_events_btag_denom_PtVsEta->GetYaxis()->FindBin(etas.at(i));

			float Numerator = h_d_enu_events_btag_numer_PtVsEta->GetBinContent(PtNum, EtaNum);
			float Denominator = h_d_enu_events_btag_denom_PtVsEta->GetBinContent(PtDenom, EtaDenom);
			float eff = Numerator / Denominator;
			BTaggedEff.push_back(eff);
		}
		return BTaggedEff;
	}};
        auto NonBTaggedBinFunction{[&h_d_enu_events_non_btag_numer_PtVsEta, &h_d_enu_events_non_btag_denom_PtVsEta](const floats& pts, const floats& etas){
                floats NonBTaggedEff{};
		int PtNum;
		int EtaNum;
		int PtDenom;
		int EtaDenom;
                for(int i = 0; i < pts.size(); i++)
                {

                        PtNum = h_d_enu_events_non_btag_numer_PtVsEta->GetXaxis()->FindBin(pts.at(i));
                        EtaNum = h_d_enu_events_non_btag_numer_PtVsEta->GetYaxis()->FindBin(etas.at(i));
                        PtDenom = h_d_enu_events_non_btag_denom_PtVsEta->GetXaxis()->FindBin(pts.at(i));
                        EtaDenom = h_d_enu_events_non_btag_denom_PtVsEta->GetYaxis()->FindBin(etas.at(i));

			float Numerator = h_d_enu_events_non_btag_numer_PtVsEta->GetBinContent(PtNum, EtaNum);
                       	float Denominator = h_d_enu_events_non_btag_denom_PtVsEta->GetBinContent(PtDenom, EtaDenom);
                       	float eff = Numerator / Denominator;
                      	NonBTaggedEff.push_back(eff);

                }
                return NonBTaggedEff;
        }};


	auto d_enu_btag_eff = d_enu_top_selection.Define("EffBTagged", BTaggedBinFunction, {"tight_jets_pt", "tight_jets_eta"})
						.Define("NonEffBTagged", NonBTaggedBinFunction, {"tight_jets_pt", "tight_jets_eta"});
	auto h_d_enu_events_btag_eff = d_enu_btag_eff.Histo1D({"MC btag EFF","MC btag EFF electro and neutrino channel",50,0,400},"EffBTagged");
        auto h_d_enu_events_non_btag_eff = d_enu_btag_eff.Histo1D({"MC non btag EFF","MC non btag EFF electro and neutrino channel",50,0,400},"NonEffBTagged");

	auto d_enu_P_btag = d_enu_btag_eff.Define("Pi_ei",EffBTaggedProduct, {"EffBTagged"})
					.Define("Pi_ej", EffNonBTaggedProduct,{"NonEffBTagged"})
					.Define("Pi_sfei",Sfi_EffBTaggedProduct, {"EffBTagged", "sfi"})
					.Define("Pi_sfej", Sfj_EffNonBTaggedProduct, {"NonEffBTagged", "sfj"})
					.Define("P_MC", P_MC_func, {"Pi_ei", "Pi_ej"});
					//.Define("P_Data", P_data_func, {"Pi_sfei", "Pi_sfej"});
					//.Define("dummy", open_CSV_file, {"Jet_pt"});
/*
	auto h_dummy_csv = d_enu_P_btag.Histo1D({"dummy","dummy", 50,0,3}, "dummy");
	h_dummy_csv->Write();
*/
	/*auto h_d_enu_Pi_ei = d_enu_P_btag.Histo1D({"Pi ei histogram","Pi ei histogram",50,0,400},"Pi_ei"); h_d_enu_Pi_ei->Write();
	auto h_d_enu_Pi_ej = d_enu_P_btag.Histo1D({"Pi ej histogram","Pi ej histogram",50,0,400},"Pi_ej"); h_d_enu_Pi_ej->Write();
	*/ auto h_d_enu_sfi = d_enu_P_btag.Histo1D({"sfi histogram","sfi histogram",50,0,400},"sfi"); h_d_enu_sfi->Write();
	//auto h_d_enu_sfj = d_enu_P_btag.Histo1D({"sfj histogram","sfj histogram",50,0,400},"sfj"); h_d_enu_sfj->Write();
	/* auto h_d_enu_Pi_sfei = d_enu_P_btag.Histo1D({"Pi sfei histogram","Pi sfei histogram",50,0,400},"Pi_sfei"); h_d_enu_Pi_sfei->Write();
	auto h_d_enu_Pi_sfej = d_enu_P_btag.Histo1D({"Pi sfej histogram","Pi sfej histogram",50,0,400},"Pi_sfej"); h_d_enu_Pi_sfej->Write();
	auto h_d_enu_P_MC = d_enu_P_btag.Histo1D({"P(MC) histogram","P(MC) histogram",50,0,400},"P_MC"); h_d_enu_P_MC->Write();
	auto h_d_enu_P_Data = d_enu_P_btag.Histo1D({"P(Data) histogram","P(Data) histogram",50,0,400},"P_Data"); //h_d_enu_P_Data->Write();
	*/
//////////////////////////////////////////////////////////////////////////// Muon Channel/////////////////////////////////////////////////////////////////////////////
        auto d_munu_event_selection = d.Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                      .Define("tight_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
				      .Define("tight_mu_eta", select<floats>, {"Muon_eta", "tight_mus"})
				      .Define("tight_mu_phi", select<floats>, {"Muon_phi", "tight_mus"})
                                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("loose_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
                                      .Define("MET_phi_Selection",{"MET_phi"})
                                      .Define("MET_mu_pt_Selection",{"MET_pt"})
                                      .Filter(mu_met_selection_function, {"MET_mu_pt_Selection"}, "MET PT CUT")
                                      .Filter(e_cut, {"tight_mu_pt", "loose_mu_pt"}, "lepton cut");

        auto d_munu_w_selection = d_munu_event_selection.Define("w_mu_eta", get_w_mu_quantity_selector("eta"))
                                                        .Define("w_mu_phi", get_w_mu_quantity_selector("phi"))
                                                        .Define("w_mu_pt", get_w_mu_quantity_selector("pt"))
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"});
                                                        //.Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto d_munu_jets_selection = d_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                                                        .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

        auto d_munu_jets_bjets_selection = d_munu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto d_munu_z_rec_selection = d_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                                                                .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                                                                .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                                                                .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                                                                .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                                                                .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
								.Define("z_mu_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_mu_eta", "tight_mu_phi"})
	                                                        .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_mu_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_mu_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                                                                //.Filter(z_mass_cut, {"z_mass"}, "z mass cut");

        auto d_munu_brec_selection = d_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto d_munu_top_selection = d_munu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_mu_pt", NormScaleFact_func, {"tight_mu_pt"})
                                                        .Define("nw_tight_mu_eta", NormScaleFact_func, {"tight_mu_eta"})
                                                        .Define("nw_tight_jets_eta", NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_mu_min_dR", NormScaleFact_func, {"jet_mu_min_dR"})
                                                        .Define("nw_z_mu_min_dR", NormScaleFact_func, {"z_mu_min_dR"})
                                                        .Define("nw_w_mu_mass", NormScaleFact_func, {"w_mu_mass"})
                                                        .Define("nw_z_mass", NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", NormScaleFact_func, {"ZW_deltaphi"});


////////////////////////////////////////////////////////////////////////// WW///////////////////////////////////////////////////////////////////////////////////////

	auto ww_enu_event_selection = ww.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("tight_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("tight_ele_eta", select<floats>, {"Electron_eta", "tight_eles"})
					.Define("tight_ele_phi", select<floats>, {"Electron_phi", "tight_eles"})
                   			.Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("loose_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_electron_pt_Selection",{"MET_pt"})
					.Filter(ele_met_selection_function, {"MET_electron_pt_Selection"}, "MET PT CUT")
					.Filter(e_cut, {"tight_ele_pt", "loose_ele_pt"}, "lepton cut");

	auto ww_enu_w_selection = ww_enu_event_selection.Define("w_e_eta", get_w_e_quantity_selector("eta"))
                 					.Define("w_e_phi", get_w_e_quantity_selector("phi"))
                 					.Define("w_e_pt", get_w_e_quantity_selector("pt"))
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"});
                 					//.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto ww_enu_jets_selection = ww_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
							.Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

	auto ww_enu_jets_bjets_selection = ww_enu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto ww_enu_z_rec_selection = ww_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                                                                .Define("z_e_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_ele_eta", "tight_ele_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_e_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_electron_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut");
                 						//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


	auto ww_enu_brec_selection = ww_enu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
							.Define("bjetpt", bjet_variable, bjet_pt_strings)
							.Define("bjeteta", bjet_variable, bjet_eta_strings)
							.Define("bjetphi", bjet_variable, bjet_phi_strings)
							.Define("nbjets", numberofbjets, {"bjets"})
							.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

	auto ww_enu_top_selection = ww_enu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_e_pt", "w_e_eta", "w_e_phi", "w_e_mass"})
							.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
							.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
							.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_ele_pt", WW_NormScaleFact_func, {"tight_ele_pt"})
                                                        .Define("nw_tight_ele_eta", WW_NormScaleFact_func, {"tight_ele_eta"})
                                                        .Define("nw_tight_jets_eta", WW_NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", WW_NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_e_min_dR", WW_NormScaleFact_func, {"jet_e_min_dR"})
                                                        .Define("nw_z_e_min_dR", WW_NormScaleFact_func, {"z_e_min_dR"})
                                                        .Define("nw_w_e_mass", WW_NormScaleFact_func, {"w_e_mass"})
                                                        .Define("nw_z_mass", WW_NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", WW_NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", WW_NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", WW_NormScaleFact_func, {"ZW_deltaphi"});

///////////////////////////////////////////////////////////////////////// Muon Channel/////////////////////////////////////////////////////////////////////////////
        auto ww_munu_event_selection = ww.Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                      .Define("tight_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
				      .Define("tight_mu_eta", select<floats>, {"Muon_eta", "tight_mus"})
				      .Define("tight_mu_phi", select<floats>, {"Muon_phi", "tight_mus"})
                                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("loose_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
                                      .Define("MET_phi_Selection",{"MET_phi"})
                                      .Define("MET_mu_pt_Selection",{"MET_pt"})
                                      .Filter(mu_met_selection_function, {"MET_mu_pt_Selection"}, "MET PT CUT")
                                      .Filter(e_cut, {"tight_mu_pt", "loose_mu_pt"}, "lepton cut");

        auto ww_munu_w_selection = ww_munu_event_selection.Define("w_mu_eta", get_w_mu_quantity_selector("eta"))
                                                        .Define("w_mu_phi", get_w_mu_quantity_selector("phi"))
                                                        .Define("w_mu_pt", get_w_mu_quantity_selector("pt"))
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"});
                                                        //.Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto ww_munu_jets_selection = ww_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                                                        .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

        auto ww_munu_jets_bjets_selection = ww_munu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto ww_munu_z_rec_selection = ww_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                                                                .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                                                                .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                                                                .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                                                                .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                                                                .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
								.Define("z_mu_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_mu_eta", "tight_mu_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_mu_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_mu_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                                                                //.Filter(z_mass_cut, {"z_mass"}, "z mass cut");

        auto ww_munu_brec_selection = ww_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto ww_munu_top_selection = ww_munu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_mu_pt", WW_NormScaleFact_func, {"tight_mu_pt"})
                                                        .Define("nw_tight_mu_eta", WW_NormScaleFact_func, {"tight_mu_eta"})
                                                        .Define("nw_tight_jets_eta", WW_NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", WW_NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_mu_min_dR", WW_NormScaleFact_func, {"jet_mu_min_dR"})
                                                        .Define("nw_z_mu_min_dR", WW_NormScaleFact_func, {"z_mu_min_dR"})
                                                        .Define("nw_w_mu_mass", WW_NormScaleFact_func, {"w_mu_mass"})
                                                        .Define("nw_z_mass", WW_NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", WW_NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", WW_NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", WW_NormScaleFact_func, {"ZW_deltaphi"});

/////////////////////////////////////////////////////////////////////////// wz Electron Channel/////////////////////////////////////////////////////////////////////////
	auto wz_enu_event_selection = wz.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("tight_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("tight_ele_eta", select<floats>, {"Electron_eta", "tight_eles"})
					.Define("tight_ele_phi", select<floats>, {"Electron_phi", "tight_eles"})
                   			.Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("loose_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_electron_pt_Selection",{"MET_pt"})
					.Filter(ele_met_selection_function, {"MET_electron_pt_Selection"}, "MET PT CUT")
					.Filter(e_cut, {"tight_ele_pt", "loose_ele_pt"}, "lepton cut");

	auto wz_enu_w_selection = wz_enu_event_selection.Define("w_e_eta", get_w_e_quantity_selector("eta"))
                 					.Define("w_e_phi", get_w_e_quantity_selector("phi"))
                 					.Define("w_e_pt", get_w_e_quantity_selector("pt"))
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"});
                 					//.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto wz_enu_jets_selection = wz_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                   					.Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

	auto wz_enu_jets_bjets_selection = wz_enu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto wz_enu_z_rec_selection = wz_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                                                                .Define("z_e_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_ele_eta", "tight_ele_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_e_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_electron_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                 						//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


	auto wz_enu_brec_selection = wz_enu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
							.Define("bjetpt", bjet_variable, bjet_pt_strings)
							.Define("bjeteta", bjet_variable, bjet_eta_strings)
							.Define("bjetphi", bjet_variable, bjet_phi_strings)
							.Define("nbjets", numberofbjets, {"bjets"})
							.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

	auto wz_enu_top_selection = wz_enu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_e_pt", "w_e_eta", "w_e_phi", "w_e_mass"})
							.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
							.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
							.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
							.Define("nw_tight_ele_pt", WZ_NormScaleFact_func, {"tight_ele_pt"})
                                                        .Define("nw_tight_ele_eta", WZ_NormScaleFact_func, {"tight_ele_eta"})
                                                        .Define("nw_tight_jets_eta", WZ_NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", WZ_NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_e_min_dR", WZ_NormScaleFact_func, {"jet_e_min_dR"})
                                                        .Define("nw_z_e_min_dR", WZ_NormScaleFact_func, {"z_e_min_dR"})
                                                        .Define("nw_w_e_mass", WZ_NormScaleFact_func, {"w_e_mass"})
                                                        .Define("nw_z_mass", WZ_NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", WZ_NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", WZ_NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", WZ_NormScaleFact_func, {"ZW_deltaphi"});

/////////////////////////////////////////////////////////////////////////wz Muon Channel/////////////////////////////////////////////////////////////////////////////
        auto wz_munu_event_selection = wz.Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                      .Define("tight_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
				      .Define("tight_mu_eta", select<floats>, {"Muon_eta", "tight_mus"})
				      .Define("tight_mu_phi", select<floats>, {"Muon_phi", "tight_mus"})
                                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("loose_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
                                      .Define("MET_phi_Selection",{"MET_phi"})
                                      .Define("MET_mu_pt_Selection",{"MET_pt"})
                                      .Filter(mu_met_selection_function, {"MET_mu_pt_Selection"}, "MET PT CUT")
                                      .Filter(e_cut, {"tight_mu_pt", "loose_mu_pt"}, "lepton cut");

        auto wz_munu_w_selection = wz_munu_event_selection.Define("w_mu_eta", get_w_mu_quantity_selector("eta"))
                                                        .Define("w_mu_phi", get_w_mu_quantity_selector("phi"))
                                                        .Define("w_mu_pt", get_w_mu_quantity_selector("pt"))
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"});
                                                        //.Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto wz_munu_jets_selection = wz_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                                                        .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

        auto wz_munu_jets_bjets_selection = wz_munu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto wz_munu_z_rec_selection = wz_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                                                                .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                                                                .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                                                                .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                                                                .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                                                                .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
								.Define("z_mu_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_mu_eta", "tight_mu_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_mu_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_mu_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
								//.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                                                                //.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


        auto wz_munu_brec_selection = wz_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto wz_munu_top_selection = wz_munu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_mu_pt", WZ_NormScaleFact_func, {"tight_mu_pt"})
                                                        .Define("nw_tight_mu_eta", WZ_NormScaleFact_func, {"tight_mu_eta"})
                                                        .Define("nw_tight_jets_eta", WZ_NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", WZ_NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_mu_min_dR", WZ_NormScaleFact_func, {"jet_mu_min_dR"})
                                                        .Define("nw_z_mu_min_dR", WZ_NormScaleFact_func, {"z_mu_min_dR"})
                                                        .Define("nw_w_mu_mass", WZ_NormScaleFact_func, {"w_mu_mass"})
                                                        .Define("nw_z_mass", WZ_NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", WZ_NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", WZ_NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", WZ_NormScaleFact_func, {"ZW_deltaphi"});

///////////////////////////////////////////////////////////////////////////ttz Electron Channel/////////////////////////////////////////////////////////////////////////
	auto ttZ_enu_event_selection = ttZ.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("tight_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("tight_ele_eta", select<floats>, {"Electron_eta", "tight_eles"})
					.Define("tight_ele_phi", select<floats>, {"Electron_phi", "tight_eles"})
                   			.Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("loose_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_electron_pt_Selection",{"MET_pt"})
					.Filter(ele_met_selection_function, {"MET_electron_pt_Selection"}, "MET PT CUT")
					.Filter(e_cut, {"tight_ele_pt", "loose_ele_pt"}, "lepton cut");

	auto ttZ_enu_w_selection = ttZ_enu_event_selection.Define("w_e_eta", get_w_e_quantity_selector("eta"))
                 					.Define("w_e_phi", get_w_e_quantity_selector("phi"))
                 					.Define("w_e_pt", get_w_e_quantity_selector("pt"))
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"});
                 					//.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto ttZ_enu_jets_selection = ttZ_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                   					.Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

	auto ttZ_enu_jets_bjets_selection = ttZ_enu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto ttZ_enu_z_rec_selection = ttZ_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                                                                .Define("z_e_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_ele_eta", "tight_ele_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_e_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_electron_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                 						//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


	auto ttZ_enu_brec_selection = ttZ_enu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
							.Define("bjetpt", bjet_variable, bjet_pt_strings)
							.Define("bjeteta", bjet_variable, bjet_eta_strings)
							.Define("bjetphi", bjet_variable, bjet_phi_strings)
							.Define("nbjets", numberofbjets, {"bjets"})
							.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

	auto ttZ_enu_top_selection = ttZ_enu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_e_pt", "w_e_eta", "w_e_phi", "w_e_mass"})
							.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
							.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
							.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_ele_pt", TTZ_NormScaleFact_func, {"tight_ele_pt"})
                                                        .Define("nw_tight_ele_eta", TTZ_NormScaleFact_func, {"tight_ele_eta"})
                                                        .Define("nw_tight_jets_eta", TTZ_NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", TTZ_NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_e_min_dR", TTZ_NormScaleFact_func, {"jet_e_min_dR"})
                                                        .Define("nw_z_e_min_dR", TTZ_NormScaleFact_func, {"z_e_min_dR"})
                                                        .Define("nw_w_e_mass", TTZ_NormScaleFact_func, {"w_e_mass"})
                                                        .Define("nw_z_mass", TTZ_NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", TTZ_NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", TTZ_NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", TTZ_NormScaleFact_func, {"ZW_deltaphi"});

/////////////////////////////////////////////////////////////////////////ttz Muon Channel/////////////////////////////////////////////////////////////////////////////
        auto ttZ_munu_event_selection = ttZ.Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                      .Define("tight_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
				      .Define("tight_mu_eta", select<floats>, {"Muon_eta", "tight_mus"})
				      .Define("tight_mu_phi", select<floats>, {"Muon_phi", "tight_mus"})
                                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("loose_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
                                      .Define("MET_phi_Selection",{"MET_phi"})
                                      .Define("MET_mu_pt_Selection",{"MET_pt"})
                                      .Filter(mu_met_selection_function, {"MET_mu_pt_Selection"}, "MET PT CUT")
                                      .Filter(e_cut, {"tight_mu_pt", "loose_mu_pt"}, "lepton cut");

        auto ttZ_munu_w_selection = ttZ_munu_event_selection.Define("w_mu_eta", get_w_mu_quantity_selector("eta"))
                                                        .Define("w_mu_phi", get_w_mu_quantity_selector("phi"))
                                                        .Define("w_mu_pt", get_w_mu_quantity_selector("pt"))
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"});
                                                        //.Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto ttZ_munu_jets_selection = ttZ_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                                                        .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

        auto ttZ_munu_jets_bjets_selection = ttZ_munu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto ttZ_munu_z_rec_selection = ttZ_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                                                                .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                                                                .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                                                                .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                                                                .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                                                                .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
								.Define("z_mu_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_mu_eta", "tight_mu_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_mu_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_mu_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                                                                //.Filter(z_mass_cut, {"z_mass"}, "z mass cut");

        auto ttZ_munu_brec_selection = ttZ_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto ttZ_munu_top_selection = ttZ_munu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_mu_pt", TTZ_NormScaleFact_func, {"tight_mu_pt"})
                                                        .Define("nw_tight_mu_eta", TTZ_NormScaleFact_func, {"tight_mu_eta"})
                                                        .Define("nw_tight_jets_eta", TTZ_NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", TTZ_NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_mu_min_dR", TTZ_NormScaleFact_func, {"jet_mu_min_dR"})
                                                        .Define("nw_z_mu_min_dR", TTZ_NormScaleFact_func, {"z_mu_min_dR"})
                                                        .Define("nw_w_mu_mass", TTZ_NormScaleFact_func, {"w_mu_mass"})
                                                        .Define("nw_z_mass", TTZ_NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", TTZ_NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", TTZ_NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", TTZ_NormScaleFact_func, {"ZW_deltaphi"});



///////////////////////////////////////////////////////////////// zz Electron channel ////////////////////////////////////////////////////////////////////////////
	auto zz_enu_event_selection = zz.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("tight_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("tight_ele_eta", select<floats>, {"Electron_eta", "tight_eles"})
					.Define("tight_ele_phi", select<floats>, {"Electron_phi", "tight_eles"})
                   			.Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("loose_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_electron_pt_Selection",{"MET_pt"})
					.Filter(ele_met_selection_function, {"MET_electron_pt_Selection"}, "MET PT CUT")
					.Filter(e_cut, {"tight_ele_pt", "loose_ele_pt"}, "lepton cut");

	auto zz_enu_w_selection = zz_enu_event_selection.Define("w_e_eta", get_w_e_quantity_selector("eta"))
                 					.Define("w_e_phi", get_w_e_quantity_selector("phi"))
                 					.Define("w_e_pt", get_w_e_quantity_selector("pt"))
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"});
                 					//.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto zz_enu_jets_selection = zz_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
							.Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

	auto zz_enu_jets_bjets_selection = zz_enu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto zz_enu_z_rec_selection = zz_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                                                                .Define("z_e_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_ele_eta", "tight_ele_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_e_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_electron_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
		       						//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


	auto zz_enu_brec_selection = zz_enu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
							.Define("bjetpt", bjet_variable, bjet_pt_strings)
							.Define("bjeteta", bjet_variable, bjet_eta_strings)
							.Define("bjetphi", bjet_variable, bjet_phi_strings)
							.Define("nbjets", numberofbjets, {"bjets"})
							.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

	auto zz_enu_top_selection = zz_enu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_e_pt", "w_e_eta", "w_e_phi", "w_e_mass"})
							.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
							.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
							.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_ele_pt", ZZ_NormScaleFact_func, {"tight_ele_pt"})
                                                        .Define("nw_tight_ele_eta", ZZ_NormScaleFact_func, {"tight_ele_eta"})
                                                        .Define("nw_tight_jets_eta", ZZ_NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", ZZ_NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_e_min_dR", ZZ_NormScaleFact_func, {"jet_e_min_dR"})
                                                        .Define("nw_z_e_min_dR", ZZ_NormScaleFact_func, {"z_e_min_dR"})
                                                        .Define("nw_w_e_mass", ZZ_NormScaleFact_func, {"w_e_mass"})
                                                        .Define("nw_z_mass", ZZ_NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", ZZ_NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", ZZ_NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", ZZ_NormScaleFact_func, {"ZW_deltaphi"});

/////////////////////////////////////////////////////////////////////////zz Muon Channel/////////////////////////////////////////////////////////////////////////////
        auto zz_munu_event_selection = zz.Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                      .Define("tight_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
				      .Define("tight_mu_eta", select<floats>, {"Muon_eta", "tight_mus"})
				      .Define("tight_mu_phi", select<floats>, {"Muon_phi", "tight_mus"})
                                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("loose_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
                                      .Define("MET_phi_Selection",{"MET_phi"})
                                      .Define("MET_mu_pt_Selection",{"MET_pt"})
                                      .Filter(mu_met_selection_function, {"MET_mu_pt_Selection"}, "MET PT CUT")
                                      .Filter(e_cut, {"tight_mu_pt", "loose_mu_pt"}, "lepton cut");

        auto zz_munu_w_selection = zz_munu_event_selection.Define("w_mu_eta", get_w_mu_quantity_selector("eta"))
                                                        .Define("w_mu_phi", get_w_mu_quantity_selector("phi"))
                                                        .Define("w_mu_pt", get_w_mu_quantity_selector("pt"))
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"});
                                                        //.Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto zz_munu_jets_selection = zz_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                                                        .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

        auto zz_munu_jets_bjets_selection = zz_munu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto zz_munu_z_rec_selection = zz_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                                                                .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                                                                .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                                                                .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                                                                .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                                                                .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
								.Define("z_mu_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_mu_eta", "tight_mu_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_mu_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_mu_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                                                                //.Filter(z_mass_cut, {"z_mass"}, "z mass cut");

        auto zz_munu_brec_selection = zz_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto zz_munu_top_selection = zz_munu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_mu_pt", ZZ_NormScaleFact_func, {"tight_mu_pt"})
                                                        .Define("nw_tight_mu_eta", ZZ_NormScaleFact_func, {"tight_mu_eta"})
                                                        .Define("nw_tight_jets_eta", ZZ_NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", ZZ_NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_mu_min_dR", ZZ_NormScaleFact_func, {"jet_mu_min_dR"})
                                                        .Define("nw_z_mu_min_dR", ZZ_NormScaleFact_func, {"z_mu_min_dR"})
                                                        .Define("nw_w_mu_mass", ZZ_NormScaleFact_func, {"w_mu_mass"})
                                                        .Define("nw_z_mass", ZZ_NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", ZZ_NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", ZZ_NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", ZZ_NormScaleFact_func, {"ZW_deltaphi"});


///////////////////////////////////////////////////////////////////////////Electron Channel/////////////////////////////////////////////////////////////////////////
/*	auto bg_enu_event_selection = bg.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("tight_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("tight_ele_eta", select<floats>, {"Electron_eta", "tight_eles"})
					.Define("tight_ele_phi", select<floats>, {"Electron_phi", "tight_eles"})
                   			.Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("loose_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_electron_pt_Selection",{"MET_pt"})
					.Filter(ele_met_selection_function, {"MET_electron_pt_Selection"}, "MET PT CUT")
					.Filter(e_cut, {"tight_ele_pt", "loose_ele_pt"}, "lepton cut");

	auto bg_enu_w_selection = bg_enu_event_selection.Define("w_e_eta", get_w_e_quantity_selector("eta"))
                 					.Define("w_e_phi", get_w_e_quantity_selector("phi"))
                 					.Define("w_e_pt", get_w_e_quantity_selector("pt"))
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"});
                 					//.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto bg_enu_jets_selection = bg_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
							.Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
							.Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
							.Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
							.Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                   					.Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

	auto bg_enu_jets_bjets_selection = bg_enu_jets_selection.Define("bjets", bjet_id, {"Jet_btagCSVV2", "Jet_eta"})
                    					        .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto bg_enu_z_rec_selection = bg_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                                                                .Define("z_e_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_ele_eta", "tight_ele_phi"})
							        .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_e_phi"})
								.Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_electron_pt_Selection"})
								.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
								.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
								.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                 						//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


	auto bg_enu_brec_selection = bg_enu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
							.Define("bjetpt", bjet_variable, bjet_pt_strings)
							.Define("bjeteta", bjet_variable, bjet_eta_strings)
							.Define("bjetphi", bjet_variable, bjet_phi_strings)
							.Define("nbjets", numberofbjets, {"bjets"})
							.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

	auto bg_enu_top_selection = bg_enu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_e_pt", "w_e_eta", "w_e_phi", "w_e_mass"})
							.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
							.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
							.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_ele_pt", NormScaleFact_func, {"tight_ele_pt"})
                                                        .Define("nw_tight_ele_eta", NormScaleFact_func, {"tight_ele_eta"})
                                                        .Define("nw_tight_jets_eta", NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_e_min_dR", NormScaleFact_func, {"jet_e_min_dR"})
                                                        .Define("nw_z_e_min_dR", NormScaleFact_func, {"z_e_min_dR"})
                                                        .Define("nw_w_e_mass", NormScaleFact_func, {"w_e_mass"})
                                                        .Define("nw_z_mass", NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", NormScaleFact_func, {"ZW_deltaphi"});



///////////////////////////////////////////////////////////////////////// Muon Channel/////////////////////////////////////////////////////////////////////////////
        auto bg_munu_event_selection = bg.Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                      .Define("tight_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
				      .Define("tight_mu_eta", select<floats>, {"Muon_eta", "tight_mus"})
				      .Define("tight_mu_phi", select<floats>, {"Muon_phi", "tight_mus"})
                                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("loose_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
                                      .Define("MET_phi_Selection",{"MET_phi"})
                                      .Define("MET_mu_pt_Selection",{"MET_pt"})
                                      .Filter(mu_met_selection_function, {"MET_mu_pt_Selection"}, "MET PT CUT")
                                      .Filter(e_cut, {"tight_mu_pt", "loose_mu_pt"}, "lepton cut");

        auto bg_munu_w_selection = bg_munu_event_selection.Define("w_mu_eta", get_w_mu_quantity_selector("eta"))
                                                        .Define("w_mu_phi", get_w_mu_quantity_selector("phi"))
                                                        .Define("w_mu_pt", get_w_mu_quantity_selector("pt"))
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"});
                                                        //.Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto bg_munu_jets_selection = bg_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                                                        .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

        auto bg_munu_jets_bjets_selection = bg_munu_jets_selection.Define("bjets", bjet_id, {"Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto bg_munu_z_rec_selection = bg_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                                                                .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                                                                .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                                                                .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                                                                .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                                                                .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
								.Define("z_mu_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_mu_eta", "tight_mu_phi"})
	                                                        .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_mu_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_mu_pt_Selection"})
                                                                .Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
                                                                .Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                .Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                                                                //.Filter(z_mass_cut, {"z_mass"}, "z mass cut");

        auto bg_munu_brec_selection = bg_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto bg_munu_top_selection = bg_munu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
                                                        .Define("nw_tight_mu_pt", NormScaleFact_func, {"tight_mu_pt"})
                                                        .Define("nw_tight_mu_eta", NormScaleFact_func, {"tight_mu_eta"})
                                                        .Define("nw_tight_jets_eta", NormScaleFact_func, {"tight_jets_eta"})
                                                        .Define("nw_tight_jets_pt", NormScaleFact_func, {"tight_jets_pt"})
                                                        .Define("nw_jet_mu_min_dR", NormScaleFact_func, {"jet_mu_min_dR"})
                                                        .Define("nw_z_mu_min_dR", NormScaleFact_func, {"z_mu_min_dR"})
                                                        .Define("nw_w_mu_mass", NormScaleFact_func, {"w_mu_mass"})
                                                        .Define("nw_z_mass", NormScaleFact_func_novec, {"z_mass"})
                                                        .Define("nw_tight_jets_deltaphi", NormScaleFact_func, {"tight_jets_deltaphi"})
                                                        .Define("nw_ZMet_deltaphi", NormScaleFact_func, {"ZMet_deltaphi"})
                                                        .Define("nw_ZW_deltaphi", NormScaleFact_func, {"ZW_deltaphi"});


*/
//////////////////////////////////////////////////////////////////////  Single Electron /////////////////////////////////////////////////////////////////////////
/*	auto se_enu_event_selection = Se.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("tight_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("tight_ele_eta", select<floats>, {"Electron_eta", "tight_eles"})
					.Define("tight_ele_phi", select<floats>, {"Electron_phi", "tight_eles"})
                   			.Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("loose_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_electron_pt_Selection",{"MET_pt"})
					.Filter(ele_met_selection_function, {"MET_electron_pt_Selection"}, "MET PT CUT")
					.Filter(e_cut, {"tight_ele_pt", "loose_ele_pt"}, "lepton cut");

	auto se_enu_w_selection = se_enu_event_selection.Define("w_e_eta", get_w_e_quantity_selector("eta"))
                 					.Define("w_e_phi", get_w_e_quantity_selector("phi"))
                 					.Define("w_e_pt", get_w_e_quantity_selector("pt"))
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"});
                 					//.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto se_enu_jets_selection = se_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                   					.Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

	auto se_enu_jets_bjets_selection = se_enu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto se_enu_z_rec_selection = se_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                                                                .Define("z_e_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_ele_eta", "tight_ele_phi"})
							        .Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_e_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_electron_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                 						//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


	auto se_enu_brec_selection = se_enu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
							.Define("bjetpt", bjet_variable, bjet_pt_strings)
							.Define("bjeteta", bjet_variable, bjet_eta_strings)
							.Define("bjetphi", bjet_variable, bjet_phi_strings)
							.Define("nbjets", numberofbjets, {"bjets"})
							.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

	auto se_enu_top_selection = se_enu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_e_pt", "w_e_eta", "w_e_phi", "w_e_mass"})
							.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
							.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
							.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});
*/

//////////////////////////////////////////////////////////////////////////////Data Single Muon///////////////////////////////////////////////////////////////////////////
/*	auto sm_munu_event_selection = Sm.Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                      .Define("tight_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
				      .Define("tight_mu_eta", select<floats>, {"Muon_eta", "tight_mus"})
				      .Define("tight_mu_phi", select<floats>, {"Muon_phi", "tight_mus"})
                                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("loose_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
                                      .Define("MET_phi_Selection",{"MET_phi"})
                                      .Define("MET_mu_pt_Selection",{"MET_pt"})
                                      .Filter(mu_met_selection_function, {"MET_mu_pt_Selection"}, "MET PT CUT")
                                      .Filter(e_cut, {"tight_mu_pt", "loose_mu_pt"}, "lepton cut");

        auto sm_munu_w_selection = sm_munu_event_selection.Define("w_mu_eta", get_w_mu_quantity_selector("eta"))
                                                        .Define("w_mu_phi", get_w_mu_quantity_selector("phi"))
                                                        .Define("w_mu_pt", get_w_mu_quantity_selector("pt"))
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"});
                                                        //.Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto sm_munu_jets_selection = sm_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                                                        .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

        auto sm_munu_jets_bjets_selection = sm_munu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto sm_munu_z_rec_selection = sm_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                                                                .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                                                                .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                                                                .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                                                                .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                                                                .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
								.Define("z_mu_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_mu_eta", "tight_mu_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_mu_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_mu_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
								//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");

        auto sm_munu_brec_selection = sm_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto sm_munu_top_selection = sm_munu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});

////////////////////////////////////////////////////////////////////////// Data MET Electron /////////////////////////////////////////////////////////////////////////
	auto met_enu_event_selection = met.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("tight_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("tight_ele_eta", select<floats>, {"Electron_eta", "tight_eles"})
					.Define("tight_ele_phi", select<floats>, {"Electron_phi", "tight_eles"})
                   			.Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                   			.Define("loose_ele_pt", select<floats>, {"Electron_pt", "tight_eles"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_electron_pt_Selection",{"MET_pt"})
					.Filter(ele_met_selection_function, {"MET_electron_pt_Selection"}, "MET PT CUT")
					.Filter(e_cut, {"tight_ele_pt", "loose_ele_pt"}, "lepton cut");

	auto met_enu_w_selection = met_enu_event_selection.Define("w_e_eta", get_w_e_quantity_selector("eta"))
                 					.Define("w_e_phi", get_w_e_quantity_selector("phi"))
                 					.Define("w_e_pt", get_w_e_quantity_selector("pt"))
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"});
                 					//.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto met_enu_jets_selection = met_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                   					.Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

	auto met_enu_jets_bjets_selection = met_enu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto met_enu_z_rec_selection = met_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                                                                .Define("z_e_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_ele_eta", "tight_ele_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_e_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_electron_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
								//.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


	auto met_enu_brec_selection = met_enu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
							.Define("bjetpt", bjet_variable, bjet_pt_strings)
							.Define("bjeteta", bjet_variable, bjet_eta_strings)
							.Define("bjetphi", bjet_variable, bjet_phi_strings)
							.Define("nbjets", numberofbjets, {"bjets"})
							.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

	auto met_enu_top_selection = met_enu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_e_pt", "w_e_eta", "w_e_phi", "w_e_mass"})
							.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
							.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
							.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});



//////////////////////////////////////////////////////////////////////////Data MET Muon//////////////////////////////////////////////////////////////////////////////////


	auto met_munu_event_selection = met.Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                      .Define("tight_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
				      .Define("tight_mu_eta", select<floats>, {"Muon_eta", "tight_mus"})
				      .Define("tight_mu_phi", select<floats>, {"Muon_phi", "tight_mus"})
                                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("loose_mu_pt", select<floats>, {"Muon_pt", "tight_mus"})
                                      .Define("MET_phi_Selection",{"MET_phi"})
                                      .Define("MET_mu_pt_Selection",{"MET_pt"})
                                      .Filter(mu_met_selection_function, {"MET_mu_pt_Selection"}, "MET PT CUT")
                                      .Filter(e_cut, {"tight_mu_pt", "loose_mu_pt"}, "lepton cut");

        auto met_munu_w_selection = met_munu_event_selection.Define("w_mu_eta", get_w_mu_quantity_selector("eta"))
                                                        .Define("w_mu_phi", get_w_mu_quantity_selector("phi"))
                                                        .Define("w_mu_pt", get_w_mu_quantity_selector("pt"))
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"});
                                                        //.Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto met_munu_jets_selection = met_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                                                        .Define("tight_jets_pt", select<floats>, {"Jet_pt", "tight_jets"})
                                                        .Define("tight_jets_eta", select<floats>, {"Jet_eta", "tight_jets"})
                                                        .Define("tight_jets_phi", select<floats>, {"Jet_phi", "tight_jets"})
                                                        .Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                                                        .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

        auto met_munu_jets_bjets_selection = met_munu_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto met_munu_z_rec_selection = met_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                                                                .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                                                                .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                                                                .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                                                                .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                                                                .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
								.Define("z_mu_min_dR", jet_lep_min_deltaR, {"z_pair_eta", "z_pair_phi", "tight_mu_eta", "tight_mu_phi"})
                                                                .Define("ZW_deltaphi", ZW_deltaphi_func, {"z_pair_phi", "w_mu_phi"})
                                                                .Define("ZMet_deltaphi", ZMet_deltaphi_func, {"z_pair_phi", "MET_mu_pt_Selection"});
                                                                //.Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
                                                                //.Filter(ZW_deltaphi_cut, {"ZW_deltaphi"}, "delta phi ZW cut")
                                                                //.Filter(ZMet_deltaphi_cut, {"ZMet_deltaphi"}, "Z met cut ");
                                                                //.Filter(z_mass_cut, {"z_mass"}, "z mass cut");

        auto met_munu_brec_selection = met_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto met_munu_top_selection = met_munu_brec_selection.Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});

*/

///////////////////////////////////////////////////////////////////// THSTACKS ///////////////////////////////////////////////////////////////////////////////

        TLegend legend_ed= TLegend(0.7,0.7,0.94,0.94);
        legend_ed.SetFillStyle(1001);
        legend_ed.SetBorderSize(1);
        legend_ed.SetFillColor(kWhite);
//        gStyle->SetOptStat(1111111);

        TLegend legend_ebg= TLegend(0.7,0.7,0.94,0.94);
        legend_ebg.SetFillStyle(1001);
        legend_ebg.SetBorderSize(1);
        legend_ebg.SetFillColor(kWhite);
//        gStyle->SetOptStat(1111111);

        TLegend legend_md= TLegend(0.7,0.7,0.94,0.94);
        legend_md.SetFillStyle(1001);
        legend_md.SetBorderSize(1);
        legend_md.SetFillColor(kWhite);
//        gStyle->SetOptStat(1111111);

        TLegend legend_mbg= TLegend(0.7,0.7,0.94,0.94);
        legend_mbg.SetFillStyle(1001);
        legend_mbg.SetBorderSize(1);
        legend_mbg.SetFillColor(kWhite);
//        gStyle->SetOptStat(1111111);


        TLegend legend_sed= TLegend(0.7,0.7,0.94,0.94);
        legend_sed.SetFillStyle(1001);
        legend_sed.SetBorderSize(1);
        legend_sed.SetFillColor(kWhite);
//        gStyle->SetOptStat(1111111);

        TLegend legend_emet= TLegend(0.7,0.7,0.94,0.94);
        legend_emet.SetFillStyle(1001);
        legend_emet.SetBorderSize(1);
        legend_emet.SetFillColor(kWhite);
//        gStyle->SetOptStat(1111111);


        TLegend legend_smd= TLegend(0.7,0.7,0.94,0.94);
        legend_smd.SetFillStyle(1001);
        legend_smd.SetBorderSize(1);
        legend_smd.SetFillColor(kWhite);
        gStyle->SetOptStat(1111111);

        TLegend legend_mmet= TLegend(0.7,0.7,0.94,0.94);
        legend_mmet.SetFillStyle(1001);
        legend_mmet.SetBorderSize(1);
        legend_mmet.SetFillColor(kWhite);
//	gStyle->SetOptStat(1111111);


	auto h_d_enu_events_ept = d_enu_top_selection.Histo1D({"MC electron_pt_enu_Channel","MC electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt", "nw_tight_ele_pt");
//	auto h_bg_enu_events_ept = bg_enu_top_selection.Histo1D({"bg electron_pt_enu_Channel","Back ground electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt", "nw_tight_ele_pt");
        auto h_d_munu_events_mupt = d_munu_top_selection.Histo1D({"MC muon_pt_munu_Channel","MC muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt","nw_tight_mu_pt");
//        auto h_bg_munu_events_mupt = bg_munu_top_selection.Histo1D({"bg muon_pt_munu_Channel","Background muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt","nw_tight_mu_pt");
//      auto h_se_enu_events_ept = se_enu_top_selection.Histo1D({"Single Electron electron_pt_enu_Channel","Single Electron electron pt in electron-neutrino channel",50,0,250}, "tight_ele_pt","nw_tight_ele_pt");
//      auto h_met_enu_events_ept = met_enu_top_selection.Histo1D({"MET electron_pt_enu_Channel","MET electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
//      auto h_smu_munu_events_mupt = sm_munu_top_selection.Histo1D({"Single Muon muon_pt_munu_Channel","Single muon muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt","nw_tight_mu_pt");
//      auto h_met_munu_events_mupt = met_munu_top_selection.Histo1D({"MET muon_pt_munu_Channel","MET muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");
	auto h_ww_enu_events_ept = ww_enu_top_selection.Histo1D({"WW electron_pt_enu_Channel","WW electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
        auto h_wz_enu_events_ept = wz_enu_top_selection.Histo1D({"WZ electron_pt_enu_Channel","WZ electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
        auto h_zz_enu_events_ept = zz_enu_top_selection.Histo1D({"ZZ electron_pt_enu_Channel","ZZ electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");
        auto h_ttZ_enu_events_ept = ttZ_enu_top_selection.Histo1D({"ttZ electron_pt_enu_Channel","ttZ electron pt in electron-neutrino channel",50,0,250},"tight_ele_pt","nw_tight_ele_pt");

	auto h_ww_munu_events_mupt = ww_munu_top_selection.Histo1D({"WW muon_pt_munu_Channel","WW muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");
        auto h_wz_munu_events_mupt = wz_munu_top_selection.Histo1D({"WZ muon_pt_munu_Channel","WZ muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");
        auto h_zz_munu_events_mupt = zz_munu_top_selection.Histo1D({"ZZ muon_pt_munu_Channel","ZZ muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");
        auto h_ttZ_munu_events_mupt = ttZ_munu_top_selection.Histo1D({"ttZ muon_pt_munu_Channel","ttZ muon pt in muon-neutrino channel",50,0,250},"tight_mu_pt", "nw_tight_mu_pt");

	THStack *lep_pt_Stack = new THStack("MC_Stack","Lepton Transverse Momentum ");
	h_d_enu_events_ept->SetLineColor(kBlack);
//	h_bg_enu_events_ept->SetLineColor(kRed);
	h_d_munu_events_mupt->SetLineColor(kGreen);
//	h_bg_munu_events_mupt->SetLineColor(kBlue);
//        h_se_enu_events_ept->SetLineColor(kPink);
//        h_met_enu_events_ept->SetLineColor(kCherry);
//        h_sm_munu_events_mupt->SetLineColor(kViolet);
//        h_met_munu_events_mupt->SetLineColor(kRose);
	h_ww_enu_events_ept->SetLineColor(kGray);
	h_wz_enu_events_ept->SetLineColor(kRed);
	h_zz_enu_events_ept->SetLineColor(kBlue);
	h_ttZ_enu_events_ept->SetLineColor(kCyan);
        h_ww_munu_events_mupt->SetLineColor(kOrange);
        h_wz_munu_events_mupt->SetLineColor(kSpring);
        h_zz_munu_events_mupt->SetLineColor(kAzure);
        h_ttZ_munu_events_mupt->SetLineColor(kTeal);



	lep_pt_Stack->Add((TH1*)&h_d_enu_events_ept.GetValue());
	//lep_pt_Stack->Add((TH1*)&h_bg_enu_events_ept.GetValue());
	lep_pt_Stack->Add((TH1*)&h_d_munu_events_mupt.GetValue());
	//lep_pt_Stack->Add((TH1*)&h_bg_munu_events_mupt.GetValue());
        //lep_pt_Stack->Add((TH1*)&h_se_enu_events_ept.GetValue());
        //lep_pt_Stack->Add((TH1*)&h_met_enu_events_ept.GetValue());
        //lep_pt_Stack->Add((TH1*)&h_smu_munu_events_mupt.GetValue());
        //lep_pt_Stack->Add((TH1*)&h_met_munu_events_mupt.GetValue());
	lep_pt_Stack->Add((TH1*)&h_ww_enu_events_ept.GetValue());
	lep_pt_Stack->Add((TH1*)&h_wz_enu_events_ept.GetValue());
	lep_pt_Stack->Add((TH1*)&h_zz_enu_events_ept.GetValue());
	lep_pt_Stack->Add((TH1*)&h_ttZ_enu_events_ept.GetValue());
        lep_pt_Stack->Add((TH1*)&h_ww_munu_events_mupt.GetValue());
        lep_pt_Stack->Add((TH1*)&h_wz_munu_events_mupt.GetValue());
        lep_pt_Stack->Add((TH1*)&h_zz_munu_events_mupt.GetValue());
        lep_pt_Stack->Add((TH1*)&h_ttZ_munu_events_mupt.GetValue());


	auto h_events_lep_pt_canvas = new TCanvas("electron pt", "electron pt",10,10,900,900);

        h_d_enu_events_ept->GetXaxis()->SetTitle("Pt/GeV");
        h_d_enu_events_ept->GetYaxis()->SetTitle("Events");
	legend_ed.AddEntry(h_d_enu_events_ept.GetPtr(),"tZq MC,electron pt","l");
	//legend_ebg.AddEntry(h_bg_enu_events_ept.GetPtr(),"tZq Background,electron pt","l");
        legend_md.AddEntry(h_d_munu_events_mupt.GetPtr(),"tZq MC,muon pt","l");
        //legend_mbg.AddEntry(h_bg_munu_events_mupt.GetPtr(),"tZq Background,muon pt","l");
        //legend_sed.AddEntry(h_se_enu_events_ept.GetPtr(),"Single electron pt","l");
        //legend_emet.AddEntry(h_se_enu_events_ept.GetPtr(),"MET data,electron pt","l");
        //legend_smd.AddEntry(h_smu_munu_events_mupt.GetPtr(),"Signal muon pt","l");
        //legend_mmet.AddEntry(h_met_munu_events_mupt.GetPtr(),"MET data,muon pt","l");


	h_d_enu_events_ept->Draw();
        //h_bg_enu_events_ept->Draw("SAME");
        h_d_munu_events_mupt->Draw("SAME");
        //h_bg_munu_events_mupt->Draw("SAME");
        //h_se_enu_events_ept->Draw("SAME");
        //h_met_enu_events_ept->Draw("SAME");
        //h_smu_munu_events_mupt->Draw("SAME");
        //h_met_munu_events_mupt->Draw("SAME");
	h_ww_enu_events_ept->Draw("SAME");
	h_wz_enu_events_ept->Draw("SAME");
	h_zz_enu_events_ept->Draw("SAME");
	h_ttZ_enu_events_ept->Draw("SAME");
        h_ww_munu_events_mupt->Draw("SAME");
        h_wz_munu_events_mupt->Draw("SAME");
        h_zz_munu_events_mupt->Draw("SAME");
        h_ttZ_munu_events_mupt->Draw("SAME");

        h_events_lep_pt_canvas->cd(2);
        lep_pt_Stack->Draw("HIST");


	h_events_lep_pt_canvas->BuildLegend();
        h_events_lep_pt_canvas->SaveAs("hist_lep_pt.root");
        h_events_lep_pt_canvas->SaveAs("hist_lep_pt.pdf");
	legend_ed.Clear();
//	legend_ebg.Clear();
       legend_md.Clear();
        //legend_mbg.Clear();
  //      legend_sed.Clear();
  //      legend_emet.Clear();
  //      legend_smd.Clear();
  //      legend_mmet.Clear();


	auto h_d_enu_events_jpt = d_enu_top_selection.Histo1D({"MC jet_pt_enu_Channel","MC jet pt in electron-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        //auto h_bg_enu_events_jpt = bg_enu_top_selection.Histo1D({"bg jet_pt_enu_Channel","Back ground jet pt in electron-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_d_munu_events_jpt = d_munu_top_selection.Histo1D({"MC jet_pt_munu_Channel","MC jet pt in muon-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        //auto h_bg_munu_events_jpt = bg_munu_top_selection.Histo1D({"Bg jet_pt_munu_Channel","Background jet pt in muon-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
//      auto h_se_enu_events_jpt = se_enu_top_selection.Histo1D({"Single Electron jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
//      auto h_met_enu_events_jpt = met_enu_top_selection.Histo1D({"MET jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
//      auto h_smu_munu_events_jpt = sm_munu_top_selection.Histo1D({"Single muon jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
//      auto h_met_munu_events_jpt = met_munu_top_selection.Histo1D({"MET jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_ww_enu_events_jpt = ww_enu_top_selection.Histo1D({"WW jet_pt_enu_Channel","WW jet pt in electron-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_wz_enu_events_jpt = wz_enu_top_selection.Histo1D({"WZ jet_pt_enu_Channel","WZ jet pt in electron-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_zz_enu_events_jpt = zz_enu_top_selection.Histo1D({"ZZ jet_pt_enu_Channel","ZZ jet pt in electron-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_ttZ_enu_events_jpt = ttZ_enu_top_selection.Histo1D({"ttZ jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_ww_munu_events_jpt = ww_munu_top_selection.Histo1D({"WW jet_pt_munu_Channel","WW jet pt in muon-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_wz_munu_events_jpt = wz_munu_top_selection.Histo1D({"WZ jet_pt_munu_Channel","WZ jet pt in muon-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_zz_munu_events_jpt = zz_munu_top_selection.Histo1D({"ZZ jet_pt_munu_Channel","ZZ jet pt in muon-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");
        auto h_ttZ_munu_events_jpt = ttZ_munu_top_selection.Histo1D({"ttZ jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,250},"tight_jets_pt","nw_tight_jets_pt");


        THStack *jet_pt_Stack = new THStack("MC_Stack","Jets Transverse Momentum ");
        h_d_enu_events_jpt->SetLineColor(kBlack);
	//h_bg_enu_events_jpt->SetLineColor(kRed);
        h_d_munu_events_jpt->SetLineColor(kGreen);
        //h_bg_munu_events_jpt->SetLineColor(kBlue);
//        h_se_enu_events_jpt->SetLineColor(kPink);
//        h_met_enu_events_jpt->SetLineColor(kCherry);
//        h_sm_munu_events_jpt->SetLineColor(kViolet);
//        h_met_munu_events_jpt->SetLineColor(kRose);
        h_ww_enu_events_jpt->SetLineColor(kGray);
        h_wz_enu_events_jpt->SetLineColor(kRed);
        h_zz_enu_events_jpt->SetLineColor(kBlue);
        h_ttZ_enu_events_jpt->SetLineColor(kCyan);
        h_ww_munu_events_jpt->SetLineColor(kOrange);
        h_wz_munu_events_jpt->SetLineColor(kSpring);
        h_zz_munu_events_jpt->SetLineColor(kAzure);
        h_ttZ_munu_events_jpt->SetLineColor(kTeal);


        jet_pt_Stack->Add((TH1*)&h_d_enu_events_jpt.GetValue());
        //jet_pt_Stack->Add((TH1*)&h_bg_enu_events_jpt.GetValue());
        jet_pt_Stack->Add((TH1*)&h_d_munu_events_jpt.GetValue());
        //jet_pt_Stack->Add((TH1*)&h_bg_munu_events_jpt.GetValue());
        //jet_pt_Stack->Add((TH1*)&h_se_enu_events_jpt.GetValue());
        //jet_pt_Stack->Add((TH1*)&h_met_enu_events_jpt.GetValue());
        //jet_pt_Stack->Add((TH1*)&h_smu_munu_events_jpt.GetValue());
        //jet_pt_Stack->Add((TH1*)&h_met_munu_events_jpt.GetValue());
	jet_pt_Stack->Add((TH1*)&h_ww_enu_events_jpt.GetValue());
	jet_pt_Stack->Add((TH1*)&h_wz_enu_events_jpt.GetValue());
	jet_pt_Stack->Add((TH1*)&h_zz_enu_events_jpt.GetValue());
	jet_pt_Stack->Add((TH1*)&h_ttZ_enu_events_jpt.GetValue());
        jet_pt_Stack->Add((TH1*)&h_ww_munu_events_jpt.GetValue());
        jet_pt_Stack->Add((TH1*)&h_wz_munu_events_jpt.GetValue());
        jet_pt_Stack->Add((TH1*)&h_zz_munu_events_jpt.GetValue());
        jet_pt_Stack->Add((TH1*)&h_ttZ_munu_events_jpt.GetValue());

        auto h_events_jet_pt_canvas = new TCanvas("e nu jet pt ", "enu jet pt",10,10,900,900);

        h_d_enu_events_jpt->GetXaxis()->SetTitle("pt / GeV");
        h_d_enu_events_jpt->GetYaxis()->SetTitle("Events");
        legend_ed.AddEntry(h_d_enu_events_jpt.GetPtr(),"tZq MC,jet pt","l");
        //legend_ebg.AddEntry(h_bg_enu_events_jpt.GetPtr(),"tZq Background,jet pt","l");
        legend_md.AddEntry(h_d_munu_events_jpt.GetPtr(),"tZq MC,jet pt","l");
        //legend_mbg.AddEntry(h_bg_munu_events_jpt.GetPtr(),"tZq Background,jet pt","l");
        //legend_sed.AddEntry(h_se_enu_events_jpt.GetPtr(),"Single electron jet pt","l");
        //legend_emet.AddEntry(h_se_enu_events_jpt.GetPtr(),"MET data,electron channel jet pt","l");
        //legend_smd.AddEntry(h_smu_munu_events_jpt.GetPtr(),"Signal muon jet pt","l");
        //legend_mmet.AddEntry(h_met_munu_events_jpt.GetPtr(),"MET data,muon channel jet pt","l");

        h_d_enu_events_jpt->Draw();
        //h_bg_enu_events_jpt->Draw("SAME");
        h_d_munu_events_jpt->Draw("SAME");
        //h_bg_munu_events_jpt->Draw("SAME");
        //h_se_enu_events_jpt->Draw("SAME");
        //h_met_enu_events_jpt->Draw("SAME");
        //h_smu_munu_events_jpt->Draw("SAME");
        //h_met_munu_events_jpt->Draw("SAME");
        h_ww_enu_events_jpt->Draw("SAME");
        h_wz_enu_events_jpt->Draw("SAME");
        h_zz_enu_events_jpt->Draw("SAME");
        h_ttZ_enu_events_jpt->Draw("SAME");
        h_ww_munu_events_jpt->Draw("SAME");
        h_wz_munu_events_jpt->Draw("SAME");
        h_zz_munu_events_jpt->Draw("SAME");
        h_ttZ_munu_events_jpt->Draw("SAME");


        h_events_jet_pt_canvas->cd(2);
        jet_pt_Stack->Draw("HIST");

	h_events_jet_pt_canvas->BuildLegend();
        h_events_jet_pt_canvas->SaveAs("hist_jet_pt.root");
        h_events_jet_pt_canvas->SaveAs("hist_jet_pt.pdf");
        legend_ed.Clear();
 //       legend_ebg.Clear();
        legend_md.Clear();
 //       legend_mbg.Clear();
 //       legend_sed.Clear();
 //       legend_emet.Clear();
 //       legend_smd.Clear();
 //       legend_mmet.Clear();

        auto h_d_enu_events_eeta = d_enu_top_selection.Histo1D({"MC electron_eta_enu_Channel","MC electron eta in electron-neutrino channel",50,-3,3},"tight_ele_eta", "nw_tight_ele_eta");
        //auto h_bg_enu_events_eeta = bg_enu_top_selection.Histo1D({"bg electron_eta_enu_Channel","Back ground electron eta in electron-neutrino channel",50,-3,3},"tight_ele_eta", "nw_tight_ele_eta");
        auto h_d_munu_events_mueta = d_munu_top_selection.Histo1D({"MC muon_eta_enu_Channel","MC muon eta in muon-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta");
        //auto h_bg_munu_events_mueta = bg_munu_top_selection.Histo1D({"Bg muon_eta_enu_Channel","Background muon eta in muon-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta");
//        auto h_se_enu_events_eeta = se_enu_top_selection.Histo1D({"Single Electron electron_eta_enu_Channel","Single Electron electron eta in electron-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta");
//        auto h_met_enu_events_eeta = met_enu_top_selection.Histo1D({"MET electron_eta_enu_Channel","MET electron eta in electron-neutrino channel",50,-3,3},"tight_ele_eta", "nw_tight_ele_eta");
//        auto h_smu_munu_events_mueta = sm_munu_top_selection.Histo1D({"Single Muon muon_eta_Channel","Single Muon Muon eta in muon-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta")$
//        auto h_met_munu_events_mueta = met_munu_top_selection.Histo1D({"MET Muon_eta_Channel","MET muon eta in muon-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta");
        auto h_ww_enu_events_eeta = ww_enu_top_selection.Histo1D({"WW electron_eta_enu_Channel","WW electron eta in electron-neutrino channel",50,-3,3},"tight_ele_eta", "nw_tight_ele_eta");
        auto h_wz_enu_events_eeta = wz_enu_top_selection.Histo1D({"WZ electron_eta_enu_Channel","WZ electron eta in electron-neutrino channel",50,-3,3},"tight_ele_eta", "nw_tight_ele_eta");
        auto h_zz_enu_events_eeta = zz_enu_top_selection.Histo1D({"ZZ electron_eta_enu_Channel","ZZ electron eta in electron-neutrino channel",50,-3,3},"tight_ele_eta", "nw_tight_ele_eta");
        auto h_ttZ_enu_events_eeta = ttZ_enu_top_selection.Histo1D({"ttZ electron_eta_enu_Channel","ttZ electron eta in electron-neutrino channel",50,-3,3},"tight_ele_eta", "nw_tight_ele_eta");
        auto h_ww_munu_events_mueta = ww_munu_top_selection.Histo1D({"WW muon_eta_enu_Channel","WW muon eta in muon-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta");
        auto h_wz_munu_events_mueta = wz_munu_top_selection.Histo1D({"WZ muon_eta_enu_Channel","WZ muon eta in muon-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta");
        auto h_zz_munu_events_mueta = zz_munu_top_selection.Histo1D({"ZZ muon_eta_enu_Channel","ZZ muon eta in muon-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta");
        auto h_ttZ_munu_events_mueta = ttZ_munu_top_selection.Histo1D({"ttZ muon_eta_enu_Channel","ttZ muon eta in muon-neutrino channel",50,-3,3}, "tight_mu_eta","nw_tight_mu_eta");



        THStack *lep_eta_Stack = new THStack("MC_Stack","Leptons Eta ");
        h_d_enu_events_eeta->SetLineColor(kBlack);
        //h_bg_enu_events_eeta->SetLineColor(kRed);
        h_d_munu_events_mueta->SetLineColor(kGreen);
        //h_bg_munu_events_mueta->SetLineColor(kBlue);
//        h_se_enu_events_eeta->SetLineColor(kPink);
//        h_met_enu_events_eeta->SetLineColor(kCherry);
//        h_sm_munu_events_mueta->SetLineColor(kViolet);
//        h_met_munu_events_mueta->SetLineColor(kRose);
        h_ww_enu_events_eeta->SetLineColor(kGray);
        h_wz_enu_events_eeta->SetLineColor(kRed);
        h_zz_enu_events_eeta->SetLineColor(kBlue);
        h_ttZ_enu_events_eeta->SetLineColor(kCyan);
        h_ww_munu_events_mueta->SetLineColor(kOrange);
        h_wz_munu_events_mueta->SetLineColor(kSpring);
        h_zz_munu_events_mueta->SetLineColor(kAzure);
        h_ttZ_munu_events_mueta->SetLineColor(kTeal);


        lep_eta_Stack->Add((TH1*)&h_d_enu_events_eeta.GetValue());
        //lep_eta_Stack->Add((TH1*)&h_bg_enu_events_eeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_d_munu_events_mueta.GetValue());
        //lep_eta_Stack->Add((TH1*)&h_bg_munu_events_mueta.GetValue());
        //lep_eta_Stack->Add((TH1*)&h_se_enu_events_eeta.GetValue());
        //lep_eta_Stack->Add((TH1*)&h_met_enu_events_eeta.GetValue());
        //lep_eta_Stack->Add((TH1*)&h_smu_munu_events_mueta.GetValue());
        //lep_eta_Stack->Add((TH1*)&h_met_munu_events_mueta.GetValue());
	lep_eta_Stack->Add((TH1*)&h_ww_enu_events_eeta.GetValue());
	lep_eta_Stack->Add((TH1*)&h_wz_enu_events_eeta.GetValue());
	lep_eta_Stack->Add((TH1*)&h_zz_enu_events_eeta.GetValue());
	lep_eta_Stack->Add((TH1*)&h_ttZ_enu_events_eeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_ww_munu_events_mueta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_wz_munu_events_mueta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_zz_munu_events_mueta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_ttZ_munu_events_mueta.GetValue());


        auto h_events_lep_eta_canvas = new TCanvas("e nu jet eta ", "enu eta",10,10,900,900);

        h_d_enu_events_eeta->GetXaxis()->SetTitle("eta");
        h_d_enu_events_eeta->GetYaxis()->SetTitle("Events");
        legend_ed.AddEntry(h_d_enu_events_eeta.GetPtr(),"tZq MC,electron eta","l");
        //legend_ebg.AddEntry(h_bg_enu_events_eeta.GetPtr(),"tZq Background,electron eta","l");
        legend_md.AddEntry(h_d_munu_events_mueta.GetPtr(),"tZq MC,muon eta","l");
        //legend_mbg.AddEntry(h_bg_munu_events_mueta.GetPtr(),"tZq Background,muon eta","l");
        //legend_sed.AddEntry(h_se_enu_events_eeta.GetPtr(),"Single electron, electron eta","l");
        //legend_emet.AddEntry(h_se_enu_events_eeta.GetPtr(),"MET data,electron channel electron eta","l");
        //legend_smd.AddEntry(h_smu_munu_events_mueta.GetPtr(),"Signal muon, muon eta","l");
        //legend_mmet.AddEntry(h_met_munu_events_mueta.GetPtr(),"MET data,muon channel muon eta","l");


        h_d_enu_events_eeta->Draw();
        //h_bg_enu_events_eeta->Draw("SAME");
        h_d_munu_events_mueta->Draw("SAME");
        //h_bg_munu_events_mueta->Draw("SAME");
        //h_se_enu_events_eeta->Draw("SAME");
        //h_met_enu_events_eeta->Draw("SAME");
        //h_smu_munu_events_mueta->Draw("SAME");
        //h_met_munu_events_mueta->Draw("SAME");
        h_ww_enu_events_eeta->Draw("SAME");
        h_wz_enu_events_eeta->Draw("SAME");
        h_zz_enu_events_eeta->Draw("SAME");
        h_ttZ_enu_events_eeta->Draw("SAME");
        h_ww_munu_events_mueta->Draw("SAME");
        h_wz_munu_events_mueta->Draw("SAME");
        h_zz_munu_events_mueta->Draw("SAME");
        h_ttZ_munu_events_mueta->Draw("SAME");

        h_events_lep_eta_canvas->cd(2);
        lep_eta_Stack->Draw("HIST");

	h_events_lep_eta_canvas->BuildLegend();
        h_events_lep_eta_canvas->SaveAs("hist_lep_eta.root");
        h_events_lep_eta_canvas->SaveAs("hist_lep_eta.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
 //       legend_mbg.Clear();
 //       legend_sed.Clear();
 //       legend_emet.Clear();
 //       legend_smd.Clear();
 //       legend_mmet.Clear();


        auto h_d_enu_events_jeta = d_enu_top_selection.Histo1D({"MC jet_eta_enu_Channel","MC jet eta in electron-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
//        auto h_bg_enu_events_jeta = bg_enu_top_selection.Histo1D({"bg jet_eta_enu_Channel","Back ground jet eta in electron-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_d_munu_events_jeta = d_munu_top_selection.Histo1D({"MC jet_eta_enu_Channel","MC jet eta in muon-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
//        auto h_bg_munu_events_jeta = bg_munu_top_selection.Histo1D({"Bg jet_eta_enu_Channel","Background jet eta in muon-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
//        auto h_se_enu_events_jeta = se_enu_top_selection.Histo1D({"Single Electron jet_eta_enu_Channel","Single Electron jet eta in electron-neutrino channel",50,-3,3}, "nw_tights_jets_eta");
//        auto h_met_enu_events_jeta = met_enu_top_selection.Histo1D({"MET jet_eta_enu_Channel","MET jet eta in electron-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
//        auto h_smu_munu_events_jeta = sm_munu_top_selection.Histo1D({"Single Muon jet_eta_Channel","Single Muon jet eta in muon-neutrino channel",50,-3,3},"nw_tights_jets_eta");
//        auto h_met_munu_events_jeta = met_munu_top_selection.Histo1D({"MET jet_eta_munu_Channel","MET jet in muon-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_ww_enu_events_jeta = ww_enu_top_selection.Histo1D({"WW jet_eta_enu_Channel","WW jet eta in electron-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_wz_enu_events_jeta = wz_enu_top_selection.Histo1D({"WZ jet_eta_enu_Channel","WZ jet eta in electron-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_zz_enu_events_jeta = zz_enu_top_selection.Histo1D({"ZZ jet_eta_enu_Channel","ZZ jet eta in electron-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_ttZ_enu_events_jeta = ttZ_enu_top_selection.Histo1D({"ttZ jet_eta_enu_Channel","ttZ jet eta in electron-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_ww_munu_events_jeta = ww_munu_top_selection.Histo1D({"WW jet_eta_munu_Channel","WW jet eta in Muon-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_wz_munu_events_jeta = wz_munu_top_selection.Histo1D({"WZ jet_eta_munu_Channel","WZ jet eta in Muon-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_zz_munu_events_jeta = zz_munu_top_selection.Histo1D({"ZZ jet_eta_munu_Channel","ZZ jet eta in Muon-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");
        auto h_ttZ_munu_events_jeta = ttZ_munu_top_selection.Histo1D({"ttZ jet_eta_munu_Channel","ttZ jet eta in Muon-neutrino channel",50,-3,3}, "tight_jets_eta","nw_tight_jets_eta");

        THStack *jet_eta_Stack = new THStack("MC_Stack","Jets Eta ");
        h_d_enu_events_jeta->SetLineColor(kBlack);
        //h_bg_enu_events_jeta->SetLineColor(kRed);
        h_d_munu_events_jeta->SetLineColor(kGreen);
        //h_bg_munu_events_jeta->SetLineColor(kBlue);
//        h_se_enu_events_jeta->SetLineColor(kPink);
//        h_met_enu_events_jeta->SetLineColor(kCherry);
//        h_sm_munu_events_jeta->SetLineColor(kViolet);
//        h_met_munu_events_jeta->SetLineColor(kRose);
        h_ww_enu_events_jeta->SetLineColor(kGray);
        h_wz_enu_events_jeta->SetLineColor(kRed);
        h_zz_enu_events_jeta->SetLineColor(kBlue);
        h_ttZ_enu_events_jeta->SetLineColor(kCyan);
        h_ww_munu_events_jeta->SetLineColor(kOrange);
        h_wz_munu_events_jeta->SetLineColor(kSpring);
        h_zz_munu_events_jeta->SetLineColor(kAzure);
        h_ttZ_munu_events_jeta->SetLineColor(kTeal);


        jet_eta_Stack->Add((TH1*)&h_d_enu_events_jeta.GetValue());
        //jet_eta_Stack->Add((TH1*)&h_bg_enu_events_jeta.GetValue());
        jet_eta_Stack->Add((TH1*)&h_d_munu_events_jeta.GetValue());
        //jet_eta_Stack->Add((TH1*)&h_bg_munu_events_jeta.GetValue());
        //jet_eta_Stack->Add((TH1*)&h_se_enu_events_jeta.GetValue());
        //jet_eta_Stack->Add((TH1*)&h_met_enu_events_jeta.GetValue());
        //jet_eta_Stack->Add((TH1*)&h_smu_munu_events_jeta.GetValue());
        //jet_eta_Stack->Add((TH1*)&h_met_munu_events_jeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_ww_enu_events_jeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_wz_enu_events_jeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_zz_enu_events_jeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_ttZ_enu_events_jeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_ww_munu_events_jeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_wz_munu_events_jeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_zz_munu_events_jeta.GetValue());
        lep_eta_Stack->Add((TH1*)&h_ttZ_munu_events_jeta.GetValue());

        auto h_events_jet_eta_canvas = new TCanvas("e nu jet eta ", "enu eta",10,10,900,900);

        h_d_enu_events_jeta->GetXaxis()->SetTitle("eta");
        h_d_enu_events_jeta->GetYaxis()->SetTitle("Events");
        legend_ed.AddEntry(h_d_enu_events_jeta.GetPtr(),"tZq MC,electron channel jet eta","l");
        //legend_ebg.AddEntry(h_bg_enu_events_jeta.GetPtr(),"tZq Background,electron channel jet eta","l");
        legend_md.AddEntry(h_d_munu_events_jeta.GetPtr(),"tZq MC,muon channel jet eta","l");
        //legend_mbg.AddEntry(h_bg_munu_events_jeta.GetPtr(),"tZq Background,muon channel jet eta","l");
        //legend_sed.AddEntry(h_se_enu_events_jeta.GetPtr(),"Single electron, jet eta","l");
        //legend_emet.AddEntry(h_se_enu_events_jeta.GetPtr(),"MET data,electron channel jet eta","l");
        //legend_smd.AddEntry(h_smu_munu_events_jeta.GetPtr(),"Signal muon, jet eta","l");
        //legend_mmet.AddEntry(h_met_munu_events_jeta.GetPtr(),"MET data,muon channel jet eta","l");

        h_d_enu_events_jeta->Draw();
        //h_bg_enu_events_jeta->Draw("SAME");
        h_d_munu_events_jeta->Draw("SAME");
        //h_bg_munu_events_jeta->Draw("SAME");
        //h_se_enu_events_jeta->Draw("SAME");
        //h_met_enu_events_jeta->Draw("SAME");
        //h_smu_munu_events_jeta->Draw("SAME");
        //h_met_munu_events_jeta->Draw("SAME");
        h_ww_enu_events_jeta->Draw("SAME");
        h_wz_enu_events_jeta->Draw("SAME");
        h_zz_enu_events_jeta->Draw("SAME");
        h_ttZ_enu_events_jeta->Draw("SAME");
        h_ww_munu_events_jeta->Draw("SAME");
        h_wz_munu_events_jeta->Draw("SAME");
        h_zz_munu_events_jeta->Draw("SAME");
        h_ttZ_munu_events_jeta->Draw("SAME");

	h_events_jet_eta_canvas->cd(2);
        jet_eta_Stack->Draw("HIST");

	h_events_jet_eta_canvas->BuildLegend();
        h_events_jet_eta_canvas->SaveAs("hist_jet_eta.root");
        h_events_jet_eta_canvas->SaveAs("hist_jet_eta.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
        //legend_mbg.Clear();
//        legend_sed.Clear();
//        legend_emet.Clear();
//        legend_smd.Clear();
//        legend_mmet.Clear();


        auto h_d_enu_events_wmass = d_enu_top_selection.Histo1D({"MC enu_w_mass","MC electron-neutrino transverse w mass",50,0,400},"w_e_mass", "nw_w_e_mass");
        //auto h_bg_enu_events_wmass = bg_enu_top_selection.Histo1D({"BG enu_w_mass","Back ground electron-neutrino transverse w mass",50,0,400}, "w_e_mass", "nw_w_e_mass");
        auto h_d_munu_events_wmass = d_munu_top_selection.Histo1D({"MC munu_w_mass","MC Muon-neutrino transverse w mass",50,0,400},"w_mu_mass","nw_w_mu_mass");
        //auto h_bg_munu_events_wmass = bg_munu_top_selection.Histo1D({"BG munu_w_mass","Background Muon-neutrino transverse w mass",50,0,400},"w_mu_mass","nw_w_mu_mass");
//        auto h_se_enu_events_wmass = se_enu_top_selection.Histo1D({"Single Electron enu_w_mass","Single Electron electron-neutrino transverse w mass",50,0,400}, "w_e_mass", "nw_w_e_mass");
//        auto h_met_enu_events_wmass = met_enu_top_selection.Histo1D({"MET enu_w_mass","MET electron-neutrino transverse w mass",50,0,400}, "w_e_mass", "nw_w_e_mass");
//        auto h_smu_munu_events_wmass = sm_munu_top_selection.Histo1D({"Single Muon munu_w_mass","Single Muon Muon-neutrino transverse w mass",50,0,400},"w_mu_mass","nw_w_mu_mass");
//        auto h_met_munu_events_wmass = met_munu_top_selection.Histo1D({"MET munu_w_mass","MET Muon-neutrino transverse w mass",50,0,400},"w_mu_mass","nw_w_mu_mass");
        auto h_ww_enu_events_wmass = ww_enu_top_selection.Histo1D({"WW enu_w_mass","WW electron-neutrino transverse w mass",50,0,400}, "w_e_mass", "nw_w_e_mass");
        auto h_wz_enu_events_wmass = wz_enu_top_selection.Histo1D({"WZ enu_w_mass","WZ electron-neutrino transverse w mass",50,0,400},"w_e_mass", "nw_w_e_mass");
        auto h_ttZ_enu_events_wmass = ttZ_enu_top_selection.Histo1D({"ttZ enu_w_mass","ttZ electron-neutrino transverse w mass",50,0,400},"w_e_mass", "nw_w_e_mass");
        auto h_zz_enu_events_wmass = zz_enu_top_selection.Histo1D({"ZZ enu_w_mass","ZZ electron-neutrino transverse w mass",50,0,400},"w_e_mass", "nw_w_e_mass");
        auto h_ww_munu_events_wmass = ww_munu_top_selection.Histo1D({"WW munu_w_mass","WW Muon-neutrino transverse w mass",50,0,400},"w_mu_mass","nw_w_mu_mass");
        auto h_wz_munu_events_wmass = wz_munu_top_selection.Histo1D({"WZ munu_w_mass","WZ Muon-neutrino transverse w mass",50,0,400},"w_mu_mass","nw_w_mu_mass");
        auto h_ttZ_munu_events_wmass = ttZ_munu_top_selection.Histo1D({"ttZ munu_w_mass","ttZ Muon-neutrino transverse w mass",50,0,400},"w_mu_mass","nw_w_mu_mass");
        auto h_zz_munu_events_wmass = zz_munu_top_selection.Histo1D({"ZZ munu_w_mass","ZZ Muon-neutrino transverse w mass",50,0,400},"w_mu_mass","nw_w_mu_mass");

        THStack *trans_wmass_Stack = new THStack("MC_Stack","W Boson Transverse Mass ");
        h_d_enu_events_wmass->SetLineColor(kBlack);
        //h_bg_enu_events_wmass->SetLineColor(kRed);
        h_d_munu_events_wmass->SetLineColor(kGreen);
        //h_bg_munu_events_wmass->SetLineColor(kBlue);
//        h_se_enu_events_wmass->SetLineColor(kPink);
//        h_met_enu_events_wmass->SetLineColor(kCherry);
//        h_sm_munu_events_wmass->SetLineColor(kViolet);
//        h_met_munu_events_wmass->SetLineColor(kRose);
        h_ww_enu_events_wmass->SetLineColor(kGray);
        h_wz_enu_events_wmass->SetLineColor(kRed);
        h_zz_enu_events_wmass->SetLineColor(kBlue);
        h_ttZ_enu_events_wmass->SetLineColor(kCyan);
        h_ww_munu_events_wmass->SetLineColor(kOrange);
        h_wz_munu_events_wmass->SetLineColor(kSpring);
        h_zz_munu_events_wmass->SetLineColor(kAzure);
        h_ttZ_munu_events_wmass->SetLineColor(kTeal);

	trans_wmass_Stack->Add((TH1*)&h_d_enu_events_wmass.GetValue());
        //trans_wmass_Stack->Add((TH1*)&h_bg_enu_events_wmass.GetValue());
        trans_wmass_Stack->Add((TH1*)&h_d_munu_events_wmass.GetValue());
        //trans_wmass_Stack->Add((TH1*)&h_bg_munu_events_wmass.GetValue());
        //trans_wmass_Stack->Add((TH1*)&h_se_enu_events_wmass.GetValue());
        //trans_wmass_Stack->Add((TH1*)&h_met_enu_events_wmass.GetValue());
        //trans_wmass_Stack->Add((TH1*)&h_smu_munu_events_wmass.GetValue());
        //trans_wmass_Stack->Add((TH1*)&h_met_munu_events_wmass.GetValue());
	trans_wmass_Stack->Add((TH1*)&h_ww_enu_events_wmass.GetValue());
	trans_wmass_Stack->Add((TH1*)&h_wz_enu_events_wmass.GetValue());
	trans_wmass_Stack->Add((TH1*)&h_zz_enu_events_wmass.GetValue());
	trans_wmass_Stack->Add((TH1*)&h_ttZ_enu_events_wmass.GetValue());
        trans_wmass_Stack->Add((TH1*)&h_ww_munu_events_wmass.GetValue());
        trans_wmass_Stack->Add((TH1*)&h_wz_munu_events_wmass.GetValue());
        trans_wmass_Stack->Add((TH1*)&h_zz_munu_events_wmass.GetValue());
        trans_wmass_Stack->Add((TH1*)&h_ttZ_munu_events_wmass.GetValue());


        auto h_events_wmass_canvas = new TCanvas("enu w mass ", "enu w mass",10,10,900,900);

        h_d_enu_events_wmass->GetXaxis()->SetTitle("mass/GeV/C^2");
        h_d_enu_events_wmass->GetYaxis()->SetTitle("Events");
        legend_ed.AddEntry(h_d_enu_events_wmass.GetPtr(),"tZq MC,electron channel transverse w mass","l");
        //legend_ebg.AddEntry(h_bg_enu_events_wmass.GetPtr(),"tZq Background,electron channel transvese w mass","l");
        legend_md.AddEntry(h_d_munu_events_wmass.GetPtr(),"tZq MC,muon channel transverse w mass","l");
        //legend_mbg.AddEntry(h_bg_munu_events_wmass.GetPtr(),"tZq Background,muon channel transverse w mass","l");
        //legend_sed.AddEntry(h_se_enu_events_wmass.GetPtr(),"Single electron, transverse w mass","l");
        //legend_emet.AddEntry(h_se_enu_events_wmass.GetPtr(),"MET data,electron transverse w mass","l");
        //legend_smd.AddEntry(h_smu_munu_events_wmass.GetPtr(),"Signal muon, transverse w mass","l");
        //legend_mmet.AddEntry(h_met_munu_events_wmass.GetPtr(),"MET data,muon transverse w mass","l");


        h_d_enu_events_wmass->Draw();
        //h_bg_enu_events_wmass->Draw("SAME");
        h_d_munu_events_wmass->Draw("SAME");
        //h_bg_munu_events_wmass->Draw("SAME");
        //h_se_enu_events_wmass->Draw("SAME");
        //h_met_enu_events_wmass->Draw("SAME");
        //h_smu_munu_events_wmass->Draw("SAME");
        //h_met_munu_events_wmass->Draw("SAME");
	h_ww_enu_events_wmass->Draw("SAME");
	h_wz_enu_events_wmass->Draw("SAME");
	h_zz_enu_events_wmass->Draw("SAME");
	h_ttZ_enu_events_wmass->Draw("SAME");
        h_ww_munu_events_wmass->Draw("SAME");
        h_wz_munu_events_wmass->Draw("SAME");
        h_zz_munu_events_wmass->Draw("SAME");
        h_ttZ_munu_events_wmass->Draw("SAME");

        h_events_wmass_canvas->cd(2);
        trans_wmass_Stack->Draw("HIST");

	h_events_wmass_canvas->BuildLegend();
        h_events_wmass_canvas->SaveAs("hist_trans_wmass.root");
        h_events_wmass_canvas->SaveAs("hist_trans_wmass.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
        //legend_mbg.Clear();
//        legend_sed.Clear();
//        legend_emet.Clear();
//        legend_smd.Clear();
//        legend_mmet.Clear();

        auto h_d_enu_events_zmass = d_enu_top_selection.Histo1D({"MC Z_mass_enu_Channel","MC Z mass in electron-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        //auto h_bg_enu_events_zmass = bg_enu_top_selection.Histo1D({"BG Z_mass_enu_Channel","Background Z mass in electron-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_d_munu_events_zmass = d_munu_top_selection.Histo1D({"MC Z_mass_munu_Channel","MC Z mass in Muon-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        //auto h_bg_munu_events_zmass = bg_munu_top_selection.Histo1D({"bg Z_mass_munu_Channel","Background Z mass in Muon-neutrino channel",50,0,450},"z_mass","nw_z_mass");
//      auto h_se_enu_events_zmass = se_enu_top_selection.Histo1D({"Single Electron Z_mass_enu_Channel","Single Electron Z mass in electron-neutrino channel",50,0,450},"z_mass","nw_z_mass");
//      auto h_met_enu_events_zmass = met_enu_top_selection.Histo1D({"MET Z_mass_enu_Channel","MET Z mass in electron-neutrino channel",50,0,450},"z_mass","nw_z_mass");
//      auto h_smu_munu_events_zmass = sm_munu_top_selection.Histo1D({"Single Muon Z_mass_enu_Channel","Single Muon Z mass in Muon-neutrino channel",50,0,450},"nw_z_mass);
//      auto h_met_munu_events_zmass = met_munu_top_selection.Histo1D({"MET Z_mass_munu_Channel","MET Z mass in Muon-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_ww_enu_events_zmass = ww_enu_top_selection.Histo1D({"WW Z_mass_enu_Channel","WW Z mass in electron-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_wz_enu_events_zmass = wz_enu_top_selection.Histo1D({"WZ Z_mass_enu_Channel","WZ Z mass in electron-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_ttZ_enu_events_zmass = ttZ_enu_top_selection.Histo1D({"ttZ Z_mass_enu_Channel","ttZ Z mass in electron-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_zz_enu_events_zmass = zz_enu_top_selection.Histo1D({"zz Z_mass_enu_Channel","zz Z mass in electron-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_ww_munu_events_zmass = ww_munu_top_selection.Histo1D({"WW Z_mass_munu_Channel","WW Z mass in Muon-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_wz_munu_events_zmass = wz_munu_top_selection.Histo1D({"WZ Z_mass_munu_Channel","WZ Z mass in Muon-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_ttZ_munu_events_zmass = ttZ_munu_top_selection.Histo1D({"ttZ Z_mass_munu_Channel","ttZ Z mass in Muon-neutrino channel",50,0,450},"z_mass","nw_z_mass");
        auto h_zz_munu_events_zmass = zz_munu_top_selection.Histo1D({"zz Z_mass_munu_Channel","zz Z mass in Muon-neutrino channel",50,0,450},"z_mass","nw_z_mass");


        THStack *zmass_Stack = new THStack("MC_Stack","Z Boson Mass ");
        h_d_enu_events_zmass->SetLineColor(kBlack);
        //h_bg_enu_events_zmass->SetLineColor(kRed);
        h_d_munu_events_zmass->SetLineColor(kGreen);
        //h_bg_munu_events_zmass->SetLineColor(kBlue);
//        h_se_enu_events_zmass->SetLineColor(kPink);
//        h_met_enu_events_zmass->SetLineColor(kCherry);
//        h_sm_munu_events_zmass->SetLineColor(kViolet);
//        h_met_munu_events_zmass->SetLineColor(kRose);
        h_ww_enu_events_zmass->SetLineColor(kGray);
        h_wz_enu_events_zmass->SetLineColor(kRed);
        h_zz_enu_events_zmass->SetLineColor(kBlue);
        h_ttZ_enu_events_zmass->SetLineColor(kCyan);
        h_ww_munu_events_zmass->SetLineColor(kOrange);
        h_wz_munu_events_zmass->SetLineColor(kSpring);
        h_zz_munu_events_zmass->SetLineColor(kAzure);
        h_ttZ_munu_events_zmass->SetLineColor(kTeal);

        zmass_Stack->Add((TH1*)&h_d_enu_events_zmass.GetValue());
        //zmass_Stack->Add((TH1*)&h_bg_enu_events_zmass.GetValue());
        zmass_Stack->Add((TH1*)&h_d_munu_events_zmass.GetValue());
        //zmass_Stack->Add((TH1*)&h_bg_munu_events_zmass.GetValue());
        //zmass_Stack->Add((TH1*)&h_se_enu_events_zmass.GetValue());
        //zmass_Stack->Add((TH1*)&h_met_enu_events_zmass.GetValue());
        //zmass_Stack->Add((TH1*)&h_smu_munu_events_zmass.GetValue());
        //zmass_Stack->Add((TH1*)&h_met_munu_events_zmass.GetValue());
	zmass_Stack->Add((TH1*)&h_ww_enu_events_zmass.GetValue());
	zmass_Stack->Add((TH1*)&h_wz_enu_events_zmass.GetValue());
	zmass_Stack->Add((TH1*)&h_zz_enu_events_zmass.GetValue());
	zmass_Stack->Add((TH1*)&h_ttZ_enu_events_zmass.GetValue());
        zmass_Stack->Add((TH1*)&h_ww_munu_events_zmass.GetValue());
        zmass_Stack->Add((TH1*)&h_wz_munu_events_zmass.GetValue());
        zmass_Stack->Add((TH1*)&h_zz_munu_events_zmass.GetValue());
        zmass_Stack->Add((TH1*)&h_ttZ_munu_events_zmass.GetValue());


        auto h_events_zmass_canvas = new TCanvas("enu z mass ", "enu z mass",10,10,900,900);

        h_d_enu_events_zmass->GetXaxis()->SetTitle("mass/GeV/C^2");
        h_d_enu_events_zmass->GetYaxis()->SetTitle("Events");
        legend_ed.AddEntry(h_d_enu_events_zmass.GetPtr(),"tZq MC,electron channel z mass","l");
        //legend_ebg.AddEntry(h_bg_enu_events_zmass.GetPtr(),"tZq Background,electron channel z mass","l");
        legend_md.AddEntry(h_d_munu_events_zmass.GetPtr(),"tZq MC,muon channel z mass","l");
        //legend_mbg.AddEntry(h_bg_munu_events_zmass.GetPtr(),"tZq Background,muon channel z mass","l");
        //legend_sed.AddEntry(h_se_enu_events_zmass.GetPtr(),"Single electron, z mass","l");
        //legend_emet.AddEntry(h_se_enu_events_zmass.GetPtr(),"MET data,electron z mass","l");
        //legend_smd.AddEntry(h_smu_munu_events_zmass.GetPtr(),"Signal muon, z mass","l");
        //legend_mmet.AddEntry(h_met_munu_events_zmass.GetPtr(),"MET data,muon z mass","l");

        h_d_enu_events_zmass->Draw();
        //h_bg_enu_events_zmass->Draw("SAME");
        h_d_munu_events_zmass->Draw("SAME");
        //h_bg_munu_events_zmass->Draw("SAME");
        //h_se_enu_events_zmass->Draw("SAME");
        //h_met_enu_events_zmass->Draw("SAME");
        //h_smu_munu_events_zmass->Draw("SAME");
        //h_met_munu_events_zmass->Draw("SAME");
	h_ww_enu_events_zmass->Draw("SAME");
	h_wz_enu_events_zmass->Draw("SAME");
	h_zz_enu_events_zmass->Draw("SAME");
	h_ttZ_enu_events_zmass->Draw("SAME");
        h_ww_munu_events_zmass->Draw("SAME");
        h_wz_munu_events_zmass->Draw("SAME");
        h_zz_munu_events_zmass->Draw("SAME");
        h_ttZ_munu_events_zmass->Draw("SAME");

        h_events_zmass_canvas->cd(2);
        zmass_Stack->Draw("HIST");


	h_events_zmass_canvas->BuildLegend();
        h_events_zmass_canvas->SaveAs("hist_zmass.root");
        h_events_zmass_canvas->SaveAs("hist_zmass.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
        //legend_mbg.Clear();
//        legend_sed.Clear();
//        legend_emet.Clear();
//        legend_smd.Clear();
//        legend_mmet.Clear();


        auto h_d_enu_events_jdphi = d_enu_top_selection.Histo1D({"MC jdphi_enu_Channel","MC jet deltaphi in electron-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        //auto h_bg_enu_events_jdphi = bg_enu_top_selection.Histo1D({"Background jdphi_enu_Channel","Back ground jet deltaphi in electron-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_d_munu_events_jdphi = d_munu_top_selection.Histo1D({"MC jdphi_munu_Channel","MC jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        //auto h_bg_munu_events_jdphi = bg_munu_top_selection.Histo1D({"bg jdphi_munu_Channel","Backgrouund jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
//        auto h_smu_munu_events_jdphi = sm_munu_top_selection.Histo1D({"Single Muon jdphi_munu_Channel","Single Muon jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
//        auto h_met_munu_events_jdphi = met_munu_top_selection.Histo1D({"MET jdphi_munu_Channel","MET jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
//        auto h_smu_munu_events_jdphi = sm_munu_top_selection.Histo1D({"Single Muon jdphi_munu_Channel","Single Muon jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
//        auto h_met_munu_events_jdphi = met_munu_top_selection.Histo1D({"MET jdphi_munu_Channel","MET jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_ww_enu_events_jdphi = ww_enu_top_selection.Histo1D({"WW jdphi_enu_Channel","WW jet deltaphi in electron-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_wz_enu_events_jdphi = wz_enu_top_selection.Histo1D({"WZ jdphi_enu_Channel","WZ jet deltaphi in electron-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_ttZ_enu_events_jdphi = ttZ_enu_top_selection.Histo1D({"ttZ jdphi_enu_Channel","ttZ jet deltaphi in electron-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_zz_enu_events_jdphi = zz_enu_top_selection.Histo1D({"zz jdphi_enu_Channel","zz jet deltaphi in electron-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_ww_munu_events_jdphi = ww_munu_top_selection.Histo1D({"WW jdphi_munu_Channel","WW jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_wz_munu_events_jdphi = wz_munu_top_selection.Histo1D({"WZ jdphi_munu_Channel","WZ jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_ttZ_munu_events_jdphi = ttZ_munu_top_selection.Histo1D({"ttZ jdphi_munu_Channel","ttZ jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");
        auto h_zz_munu_events_jdphi = zz_munu_top_selection.Histo1D({"zz jdphi_munu_Channel","zz jet deltaphi in Muon-neutrino channel",50,0,2},"tight_jets_deltaphi","nw_tight_jets_deltaphi");


        THStack *jdphi_Stack = new THStack("MC_Stack","Jets Delta Phi ");
        h_d_enu_events_jdphi->SetLineColor(kBlack);
        //h_bg_enu_events_jdphi->SetLineColor(kRed);
        h_d_munu_events_jdphi->SetLineColor(kGreen);
        //h_bg_munu_events_jdphi->SetLineColor(kBlue);
//        h_se_enu_events_jdphi->SetLineColor(kPink);
//        h_met_enu_events_jdphi->SetLineColor(kCherry);
//        h_sm_munu_events_jdphi->SetLineColor(kViolet);
//        h_met_munu_events_jdphi->SetLineColor(kRose);
        h_ww_enu_events_jdphi->SetLineColor(kGray);
        h_wz_enu_events_jdphi->SetLineColor(kRed);
        h_zz_enu_events_jdphi->SetLineColor(kBlue);
        h_ttZ_enu_events_jdphi->SetLineColor(kCyan);
        h_ww_munu_events_jdphi->SetLineColor(kOrange);
        h_wz_munu_events_jdphi->SetLineColor(kSpring);
        h_zz_munu_events_jdphi->SetLineColor(kAzure);
        h_ttZ_munu_events_jdphi->SetLineColor(kTeal);

        jdphi_Stack->Add((TH1*)&h_d_enu_events_jdphi.GetValue());
        //jdphi_Stack->Add((TH1*)&h_bg_enu_events_jdphi.GetValue());
        jdphi_Stack->Add((TH1*)&h_d_munu_events_jdphi.GetValue());
        //jdphi_Stack->Add((TH1*)&h_bg_munu_events_jdphi.GetValue());
        //jdphi_Stack->Add((TH1*)&h_se_enu_events_jdphi.GetValue());
        //jdphi_Stack->Add((TH1*)&h_met_enu_events_jdphi.GetValue());
        //jdphi_Stack->Add((TH1*)&h_smu_munu_events_jdphi.GetValue());
        //jdphi_Stack->Add((TH1*)&h_met_munu_events_jdphi.GetValue());
	jdphi_Stack->Add((TH1*)&h_ww_enu_events_jdphi.GetValue());
	jdphi_Stack->Add((TH1*)&h_wz_enu_events_jdphi.GetValue());
	jdphi_Stack->Add((TH1*)&h_zz_enu_events_jdphi.GetValue());
	jdphi_Stack->Add((TH1*)&h_ttZ_enu_events_jdphi.GetValue());
        jdphi_Stack->Add((TH1*)&h_ww_munu_events_jdphi.GetValue());
        jdphi_Stack->Add((TH1*)&h_wz_munu_events_jdphi.GetValue());
        jdphi_Stack->Add((TH1*)&h_zz_munu_events_jdphi.GetValue());
        jdphi_Stack->Add((TH1*)&h_ttZ_munu_events_jdphi.GetValue());

        auto h_events_jdphi_canvas = new TCanvas("jetdeltaphi", "jetdeltaphi",10,10,900,900);

        h_d_enu_events_jdphi->GetXaxis()->SetTitle("jets delta phi / rad");
        h_d_enu_events_jdphi->GetYaxis()->SetTitle("Events");

        legend_ed.AddEntry(h_d_enu_events_jdphi.GetPtr(),"tZq MC,electron channel jets deltaphi","l");
        //legend_ebg.AddEntry(h_bg_enu_events_jdphi.GetPtr(),"tZq Bakground,electron channel jets deltaphi","l");
        legend_md.AddEntry(h_d_munu_events_jdphi.GetPtr(),"tZq MC,muon channel jets deltaphi","l");
        //legend_mbg.AddEntry(h_bg_munu_events_jdphi.GetPtr(),"tZq Background,muon channel jets deltaphi","l");
        //legend_sed.AddEntry(h_se_enu_events_jdphi.GetPtr(),"Single electron, jets deltaphi","l");
        //legend_emet.AddEntry(h_se_enu_events_jdphi.GetPtr(),"MET data,electron channel jets deltaphi","l");
        //legend_smd.AddEntry(h_se_munu_events_jdphi.GetPtr(),"Signal muon, jets deltaphi","l");
        //legend_mmet.AddEntry(h_met_munu_events_jdphi.GetPtr(),"MET data,muon channel jets deltaphi","l");

        h_d_enu_events_jdphi->Draw();
        //h_bg_enu_events_jdphi->Draw("SAME");
        h_d_munu_events_jdphi->Draw("SAME");
        //h_bg_munu_events_jdphi->Draw("SAME");
        //h_se_enu_events_jdphi->Draw("SAME");
        //h_met_enu_events_jdphi->Draw("SAME");
        //h_smu_munu_events_jdphi->Draw("SAME");
        //h_met_munu_events_jdphi->Draw("SAME");
	h_ww_enu_events_jdphi->Draw("SAME");
	h_wz_enu_events_jdphi->Draw("SAME");
	h_zz_enu_events_jdphi->Draw("SAME");
	h_ttZ_enu_events_jdphi->Draw("SAME");
	h_ww_munu_events_jdphi->Draw("SAME");
	h_wz_munu_events_jdphi->Draw("SAME");
	h_zz_munu_events_jdphi->Draw("SAME");
	h_ttZ_munu_events_jdphi->Draw("SAME");

        h_events_jdphi_canvas->cd(2);
        jdphi_Stack->Draw("HIST");

	h_events_jdphi_canvas->BuildLegend();
        h_events_jdphi_canvas->SaveAs("hist_jdphi.root");
        h_events_jdphi_canvas->SaveAs("hist_jdphi.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
        //legend_mbg.Clear();
//        legend_sed.Clear();
//        legend_emet.Clear();
//        legend_smd.Clear();
//        legend_mmet.Clear();


        auto h_d_enu_events_zmetdphi = d_enu_top_selection.Histo1D({"MC zmetdphi_enu_Channel","MC Z met pt deltaphi in electron-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        //auto h_bg_enu_events_zmetdphi = bg_enu_top_selection.Histo1D({"Background zmetdphi_enu_Channel","Background Z met pt deltaphi in electron-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_d_munu_events_zmetdphi = d_munu_top_selection.Histo1D({"MC zmetdphi_munu_Channel","MC Z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        //auto h_bg_munu_events_zmetdphi = bg_munu_top_selection.Histo1D({"BG zmetdphi_munu_Channel","Background Z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
//        auto h_se_enu_events_zmetdphi = se_enu_top_selection.Histo1D({"Single Electron zmetdphi_enu_Channel","Single Electron z met pt deltaphi in electron-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
//       auto h_met_enu_events_zmetdphi = met_enu_top_selection.Histo1D({"MET zmetdphi_enu_Channel","MET z met pt deltaphi in electron-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
//        auto h_smu_munu_events_zmetdphi = sm_munu_top_selection.Histo1D({"Single Muon zmetdphi_munu_Channel","Single Muon z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
//        auto h_met_munu_events_zmetdphi = met_munu_top_selection.Histo1D({"MET zmetdphi_munu_Channel","MET z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_ww_enu_events_zmetdphi = ww_enu_top_selection.Histo1D({"WW zmetdphi_enu_Channel","WW Z met pt deltaphi in electron-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_wz_enu_events_zmetdphi = wz_enu_top_selection.Histo1D({"WZ zmetdphi_enu_Channel","WZ Z met pt deltaphi in electron-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_ttZ_enu_events_zmetdphi = ttZ_enu_top_selection.Histo1D({"ttZ zmetdphi_enu_Channel","ttZ Z met pt jet deltaphi in electron-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_zz_enu_events_zmetdphi = zz_enu_top_selection.Histo1D({"zz zmetdphi_enu_Channel","zz z met pt deltaphi in electron-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_ww_munu_events_zmetdphi = ww_munu_top_selection.Histo1D({"WW zmetdphi_munu_Channel","WW Z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_wz_munu_events_zmetdphi = wz_munu_top_selection.Histo1D({"WZ zmetdphi_munu_Channel","WZ Z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_ttZ_munu_events_zmetdphi = ttZ_munu_top_selection.Histo1D({"ttZ zmetdphi_munu_Channel","ttZ Z met pt jet deltaphi in Muon-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");
        auto h_zz_munu_events_zmetdphi = zz_munu_top_selection.Histo1D({"zz zmetdphi_munu_Channel","zz z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZMet_deltaphi","nw_ZMet_deltaphi");


        THStack *zmetdphi_Stack = new THStack("MC_Stack","Z & MET delta phi ");
        h_d_enu_events_zmetdphi->SetLineColor(kBlack);
        //h_bg_enu_events_zmetdphi->SetLineColor(kRed);
        h_d_munu_events_zmetdphi->SetLineColor(kGreen);
        //h_bg_munu_events_zmetdphi->SetLineColor(kBlue);
//        h_se_enu_events_zmetdphi->SetLineColor(kPink);
//        h_met_enu_events_zmetdphi->SetLineColor(kCherry);
//        h_sm_munu_events_zmetdphi->SetLineColor(kViolet);
//        h_met_munu_events_zmetdphi->SetLineColor(kRose);
        h_ww_enu_events_zmetdphi->SetLineColor(kGray);
        h_wz_enu_events_zmetdphi->SetLineColor(kRed);
        h_zz_enu_events_zmetdphi->SetLineColor(kBlue);
        h_ttZ_enu_events_zmetdphi->SetLineColor(kCyan);
        h_ww_munu_events_zmetdphi->SetLineColor(kOrange);
        h_wz_munu_events_zmetdphi->SetLineColor(kSpring);
        h_zz_munu_events_zmetdphi->SetLineColor(kAzure);
        h_ttZ_munu_events_zmetdphi->SetLineColor(kTeal);

        zmetdphi_Stack->Add((TH1*)&h_d_enu_events_zmetdphi.GetValue());
        //zmetdphi_Stack->Add((TH1*)&h_bg_enu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_d_munu_events_zmetdphi.GetValue());
        //zmetdphi_Stack->Add((TH1*)&h_bg_munu_events_zmetdphi.GetValue());
        //zmetdphi_Stack->Add((TH1*)&h_se_enu_events_zmetdphi.GetValue());
        //zmetdphi_Stack->Add((TH1*)&h_met_enu_events_zmetdphi.GetValue());
        //zmetdphi_Stack->Add((TH1*)&h_smu_munu_events_zmetdphi.GetValue());
        //zmetdphi_Stack->Add((TH1*)&h_met_munu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_ww_enu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_wz_enu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_zz_enu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_ttZ_enu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_ww_munu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_wz_munu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_zz_munu_events_zmetdphi.GetValue());
        zmetdphi_Stack->Add((TH1*)&h_ttZ_munu_events_zmetdphi.GetValue());

        auto h_events_zmetdphi_canvas = new TCanvas("zmetdeltaphi", "zmetdeltaphi",10,10,900,900);

        h_d_enu_events_zmetdphi->GetXaxis()->SetTitle("delta phi between Z and MET/ rad");
        h_d_enu_events_zmetdphi->GetYaxis()->SetTitle("Events");

        legend_ed.AddEntry(h_d_enu_events_zmetdphi.GetPtr(),"tZq MC,electron channel deltaphi between Z and MET","l");
        //legend_ebg.AddEntry(h_bg_enu_events_zmetdphi.GetPtr(),"tZq Background,electron channel deltaphi between Z and MET","l");
        legend_md.AddEntry(h_d_munu_events_zmetdphi.GetPtr(),"tZq MC,muon channel deltaphi between Z and MET","l");
        //legend_mbg.AddEntry(h_bg_munu_events_zmetdphi.GetPtr(),"tZq Background,muon channel deltaphi between Z and MET","l");
        //legend_sed.AddEntry(h_se_enu_events_zmetdphi.GetPtr(),"Single electron, deltaphi between Z and MET","l");
        //legend_emet.AddEntry(h_se_enu_events_zmetdphi.GetPtr(),"MET data,electron channel deltaphi between Z and MET","l");
        //legend_smd.AddEntry(h_se_munu_events_zmetdphi.GetPtr(),"Signal muon,deltaphi between Z and MET","l");
        //legend_mmet.AddEntry(h_met_munu_events_zmetdphi.GetPtr(),"MET data,muon channel deltaphi between Z and MET","l");


        h_d_enu_events_zmetdphi->Draw();
        //h_bg_enu_events_zmetdphi->Draw("SAME");
        h_d_munu_events_zmetdphi->Draw("SAME");
        //h_bg_munu_events_zmetdphi->Draw("SAME");
        //h_se_enu_events_zmetdphi->Draw("SAME");
        //h_met_enu_events_zmetdphi->Draw("SAME");
        //h_smu_munu_events_zmetdphi->Draw("SAME");
        //h_met_munu_events_zmetdphi->Draw("SAME");
        h_ww_enu_events_zmetdphi->Draw("SAME");
        h_wz_enu_events_zmetdphi->Draw("SAME");
        h_zz_enu_events_zmetdphi->Draw("SAME");
        h_ttZ_enu_events_zmetdphi->Draw("SAME");
        h_ww_munu_events_zmetdphi->Draw("SAME");
        h_wz_munu_events_zmetdphi->Draw("SAME");
        h_zz_munu_events_zmetdphi->Draw("SAME");
        h_ttZ_munu_events_zmetdphi->Draw("SAME");

        h_events_zmetdphi_canvas->cd(2);
        zmetdphi_Stack->Draw("HIST");

	h_events_zmetdphi_canvas->BuildLegend();
        h_events_zmetdphi_canvas->SaveAs("hist_zmetdphi.root");
        h_events_zmetdphi_canvas->SaveAs("hist_zmetdphi.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
        //legend_mbg.Clear();
//        legend_sed.Clear();
//        legend_emet.Clear();
//        legend_smd.Clear();
//        legend_mmet.Clear();



        auto h_d_enu_events_zwdphi = d_enu_top_selection.Histo1D({"MC zwdphi_enu_Channel","MC Z w deltaphi in electron-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        //auto h_bg_enu_events_zwdphi = bg_enu_top_selection.Histo1D({"bg zwdphi_enu_Channel","Background Z w deltaphi in electron-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_d_munu_events_zwdphi = d_munu_top_selection.Histo1D({"MC zwdphi_munu_Channel","MC Z w deltaphi in Muon-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        //auto h_bg_munu_events_zwdphi = bg_munu_top_selection.Histo1D({"BG zwdphi_munu_Channel","Background Z w deltaphi in Muon-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
//        auto h_se_enu_events_zwdphi = se_enu_top_selection.Histo1D({"Single Electron zmetdphi_enu_Channel","Single Electron z met pt deltaphi in electron-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
//        auto h_met_enu_events_zwdphi = met_enu_top_selection.Histo1D({"MET zmetdphi_enu_Channel","MET z met pt deltaphi in electron-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
//        auto h_smu_munu_events_zwdphi = sm_munu_top_selection.Histo1D({"Single Muon zmetdphi_munu_Channel","Single Muon z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
//        auto h_met_munu_events_zwdphi = met_munu_top_selection.Histo1D({"MET zmetdphi_munu_Channel","MET z met pt deltaphi in Muon-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_ww_enu_events_zwdphi = ww_enu_top_selection.Histo1D({"WW zwdphi_enu_Channel","WW Z w deltaphi in electron-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_wz_enu_events_zwdphi = wz_enu_top_selection.Histo1D({"WZ zwdphi_enu_Channel","WZ Z w deltaphi in electron-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_ttZ_enu_events_zwdphi = ttZ_enu_top_selection.Histo1D({"ttZ zwdphi_enu_Channel","ttZ  Z w deltaphi in electron-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_zz_enu_events_zwdphi = zz_enu_top_selection.Histo1D({"zz zwdphi_enu_Channel","zz  Z w deltaphi in electron-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_ww_munu_events_zwdphi = ww_munu_top_selection.Histo1D({"WW zwdphi_munu_Channel","WW Z w deltaphi in Muon-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_wz_munu_events_zwdphi = wz_munu_top_selection.Histo1D({"WZ zwdphi_munu_Channel","WZ Z w deltaphi in Muon-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_ttZ_munu_events_zwdphi = ttZ_munu_top_selection.Histo1D({"ttZ zwdphi_munu_Channel","ttZ  Z w deltaphi in Muon-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");
        auto h_zz_munu_events_zwdphi = zz_munu_top_selection.Histo1D({"zz zwdphi_munu_Channel","zz z Z w deltaphi in Muon-neutrino channel",50,0,2},"ZW_deltaphi","nw_ZW_deltaphi");


        THStack *zwdphi_Stack = new THStack("MC_Stack","Z & W delta phi ");
        h_d_enu_events_zwdphi->SetLineColor(kBlack);
        //h_bg_enu_events_zwdphi->SetLineColor(kRed);
        h_d_munu_events_zwdphi->SetLineColor(kGreen);
        //h_bg_munu_events_zwdphi->SetLineColor(kBlue);
//        h_se_enu_events_zwdphi->SetLineColor(kPink);
//        h_met_enu_events_zwdphi->SetLineColor(kCherry);
//        h_sm_munu_events_zwdphi->SetLineColor(kViolet);
//        h_met_munu_events_zwdphi->SetLineColor(kRose);
        h_ww_enu_events_zwdphi->SetLineColor(kGray);
        h_wz_enu_events_zwdphi->SetLineColor(kRed);
        h_zz_enu_events_zwdphi->SetLineColor(kBlue);
        h_ttZ_enu_events_zwdphi->SetLineColor(kCyan);
        h_ww_munu_events_zwdphi->SetLineColor(kOrange);
        h_wz_munu_events_zwdphi->SetLineColor(kSpring);
        h_zz_munu_events_zwdphi->SetLineColor(kAzure);
        h_ttZ_munu_events_zwdphi->SetLineColor(kTeal);


        zwdphi_Stack->Add((TH1*)&h_d_enu_events_zwdphi.GetValue());
        //zwdphi_Stack->Add((TH1*)&h_bg_enu_events_zwdphi.GetValue());
        zwdphi_Stack->Add((TH1*)&h_d_munu_events_zwdphi.GetValue());
        //zwdphi_Stack->Add((TH1*)&h_bg_munu_events_zwdphi.GetValue());
        //zwdphi_Stack->Add((TH1*)&h_se_enu_events_zwdphi.GetValue());
        //zwdhi_Stack->Add((TH1*)&h_met_enu_events_zwdphi.GetValue());
        //zwdphi_Stack->Add((TH1*)&h_smu_munu_events_zwdphi.GetValue());
        //zwdphi_Stack->Add((TH1*)&h_met_munu_events_zwdphi.GetValue());
        zwdphi_Stack->Add((TH1*)&h_ww_enu_events_zwdphi.GetValue());
	zwdphi_Stack->Add((TH1*)&h_wz_enu_events_zwdphi.GetValue());
	zwdphi_Stack->Add((TH1*)&h_zz_enu_events_zwdphi.GetValue());
        zwdphi_Stack->Add((TH1*)&h_ttZ_enu_events_zwdphi.GetValue());
        zwdphi_Stack->Add((TH1*)&h_ww_munu_events_zwdphi.GetValue());
        zwdphi_Stack->Add((TH1*)&h_wz_munu_events_zwdphi.GetValue());
        zwdphi_Stack->Add((TH1*)&h_zz_munu_events_zwdphi.GetValue());
        zwdphi_Stack->Add((TH1*)&h_ttZ_munu_events_zwdphi.GetValue());


        auto h_events_zwdphi_canvas = new TCanvas("zwdeltaphi", "zwdeltaphi",10,10,900,900);

        h_d_enu_events_zwdphi->GetXaxis()->SetTitle("delta phi between Z and W/ rad");
        h_d_enu_events_zwdphi->GetYaxis()->SetTitle("Events");

        legend_ed.AddEntry(h_d_enu_events_zwdphi.GetPtr(),"tZq MC,electron channel deltaphi between Z and W","l");
        //legend_ebg.AddEntry(h_bg_enu_events_zwdphi.GetPtr(),"tZq Background,electron channel deltaphi between Z and W","l");
        legend_md.AddEntry(h_d_munu_events_zwdphi.GetPtr(),"tZq MC,muon channel deltaphi between Z and W","l");
        //legend_mbg.AddEntry(h_bg_munu_events_zwdphi.GetPtr(),"tZq Background,muon channel deltaphi between Z and W","l");
        //legend_sed.AddEntry(h_se_enu_events_zwdphi.GetPtr(),"Single electron, deltaphi between Z and W","l");
        //legend_emet.AddEntry(h_se_enu_events_zwdphi.GetPtr(),"MET data,electron channel deltaphi between Z and W","l");
        //legend_smd.AddEntry(h_se_munu_events_zwdphi.GetPtr(),"Signal muon,deltaphi between Z and W","l");
        //legend_mmet.AddEntry(h_met_munu_events_zwdphi.GetPtr(),"MET data,muon channel deltaphi between Z and W","l");


        h_d_enu_events_zwdphi->Draw();
        //h_bg_enu_events_zwdphi->Draw("SAME");
        h_d_munu_events_zwdphi->Draw("SAME");
        //h_bg_munu_events_zwdphi->Draw("SAME");
        //h_se_enu_events_zwdphi->Draw("SAME");
        //h_met_enu_events_zwdphi->Draw("SAME");
        //h_smu_munu_events_zwdphi->Draw("SAME");
        //h_met_munu_events_zwdphi->Draw("SAME");
        h_ww_enu_events_zwdphi->Draw("SAME");
        h_wz_enu_events_zwdphi->Draw("SAME");
	h_zz_enu_events_zwdphi->Draw("SAME");
	h_ttZ_enu_events_zwdphi->Draw("SAME");
        h_ww_munu_events_zwdphi->Draw("SAME");
        h_wz_munu_events_zwdphi->Draw("SAME");
        h_zz_munu_events_zwdphi->Draw("SAME");
        h_ttZ_munu_events_zwdphi->Draw("SAME");

        h_events_zwdphi_canvas->cd(2);
        zwdphi_Stack->Draw("HIST");

	h_events_zwdphi_canvas->BuildLegend();
        h_events_zwdphi_canvas->SaveAs("hist_zwdphi.root");
        h_events_zwdphi_canvas->SaveAs("hist_zwdphi.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
        //legend_mbg.Clear();
//        legend_sed.Clear();
//        legend_emet.Clear();
//        legend_smd.Clear();
//        legend_mmet.Clear();


        auto h_d_enu_events_ejdr = d_enu_top_selection.Histo1D({"MC ejdr_enu_Channel","MC e and jet deltaR in electron-neutrino channel",50,0,1},"jet_e_min_dR","nw_jet_e_min_dR");
        //auto h_bg_enu_events_ejdr = bg_enu_top_selection.Histo1D({"bg ejdr_enu_Channel","Background e and jet deltaR in electron-neutrino channel",50,0,1},"jet_e_min_dR","nw_jet_e_min_dR");
        auto h_d_munu_events_mujdr = d_munu_top_selection.Histo1D({"MC ejdr_munu_Channel","MC muon and jet deltaR in muon-neutrino channel",50,0,1},"jet_mu_min_dR","nw_jet_mu_min_dR");
        //auto h_bg_munu_events_mujdr = bg_munu_top_selection.Histo1D({"BG ejdr_munu_Channel","Background muon and jet deltaR in muon-neutrino channel",50,0,1},"jet_mu_min_dR","nw_jet_mu_min_dR");
//        auto h_se_enu_events_ejdr = se_enu_top_selection.Histo1D({"Single Electron ejdr_enu_Channel","Single Electron e and jet pt deltaR in electron-neutrino channel",50,0,1},"jet_e_min_dR","nw_jet_e_min_dR");
//        auto h_met_enu_events_ejdr = met_enu_top_selection.Histo1D({"MET ejdr_enu_Channel","MET e and jet deltaR in electron-neutrino channel",50,0,1},"jet_e_min_dR","nw_jet_e_min_dR");
//        auto h_smu_munu_events_mujdr = sm_munu_top_selection.Histo1D({"Single Muon mujdr_munu_Channel","Single Muon muon and jet pt deltaR in Muon-neutrino channel",50,0,1},"jet_mu_min_dR","nw_jet_mu_min_dR");
//        auto h_met_munu_events_mujdr = met_munu_top_selection.Histo1D({"MET mujdr_munu_Channel","MET muon and jet deltaR in muon-neutrino channel",50,0,1},"jet_mu_min_dR","nw_jet_mu_min_dR");
        auto h_ww_enu_events_ejdr = ww_enu_top_selection.Histo1D({"WW ejdr_enu_Channel","WW e and jet deltaR in electron-neutrino channel",50,0,1},"jet_e_min_dR","nw_jet_e_min_dR");
        auto h_wz_enu_events_ejdr = wz_enu_top_selection.Histo1D({"WZ ejdr_enu_Channel","WZ e and jet deltaR in electron-neutrino channel",50,0,1},"jet_e_min_dR","nw_jet_e_min_dR");
        auto h_ttZ_enu_events_ejdr = ttZ_enu_top_selection.Histo1D({"ttZ ejdr_enu_Channel","ttZ  e and jet deltaR in electron-neutrino channel",50,0,1},"jet_e_min_dR","nw_jet_e_min_dR");
        auto h_zz_enu_events_ejdr = zz_enu_top_selection.Histo1D({"zz ejdr_enu_Channel","zz e and jet deltaR in electron-neutrino channel",50,0,1},"jet_e_min_dR","nw_jet_e_min_dR");
        auto h_ww_munu_events_mujdr = ww_munu_top_selection.Histo1D({"WW ejdr_munu_Channel","WW muon and jet deltaR in muon-neutrino channel",50,0,1},"jet_mu_min_dR","nw_jet_mu_min_dR");
        auto h_wz_munu_events_mujdr = wz_munu_top_selection.Histo1D({"WZ ejdr_munu_Channel","WZ muon and jet deltaR in muon-neutrino channel",50,0,1},"jet_mu_min_dR","nw_jet_mu_min_dR");
        auto h_ttZ_munu_events_mujdr = ttZ_munu_top_selection.Histo1D({"ttZ ejdr_munu_Channel","ttZ muon and jet deltaR in muon-neutrino channel",50,0,1},"jet_mu_min_dR","nw_jet_mu_min_dR");
        auto h_zz_munu_events_mujdr = zz_munu_top_selection.Histo1D({"zz ejdr_munu_Channel","zz muon and jet deltaR in muon-neutrino channel",50,0,1},"jet_mu_min_dR","nw_jet_mu_min_dR");



        THStack *lepjdr_Stack = new THStack("MC_Stack","Lepton and Jet Delta R ");
        h_d_enu_events_ejdr->SetLineColor(kBlack);
        //h_bg_enu_events_ejdr->SetLineColor(kRed);
        h_d_munu_events_mujdr->SetLineColor(kGreen);
        //h_bg_munu_events_mujdr->SetLineColor(kBlue);
 //       h_se_enu_events_ejdr->SetLineColor(kPink);
 //       h_met_enu_events_ejdr->SetLineColor(kCherry);
 //       h_sm_munu_events_mujdr->SetLineColor(kViolet);
 //       h_met_munu_events_mujdr->SetLineColor(kRose);
        h_ww_enu_events_ejdr->SetLineColor(kGray);
        h_wz_enu_events_ejdr->SetLineColor(kRed);
        h_zz_enu_events_ejdr->SetLineColor(kBlue);
        h_ttZ_enu_events_ejdr->SetLineColor(kCyan);
        h_ww_munu_events_mujdr->SetLineColor(kOrange);
        h_wz_munu_events_mujdr->SetLineColor(kSpring);
        h_zz_munu_events_mujdr->SetLineColor(kAzure);
        h_ttZ_munu_events_mujdr->SetLineColor(kTeal);


        lepjdr_Stack->Add((TH1*)&h_d_enu_events_ejdr.GetValue());
        //lepjdr_Stack->Add((TH1*)&h_bg_enu_events_ejdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_d_munu_events_mujdr.GetValue());
        //lepjdr_Stack->Add((TH1*)&h_bg_munu_events_mujdr.GetValue());
        //lepjdr_Stack->Add((TH1*)&h_se_enu_events_ejdr.GetValue());
        //lepjdr_Stack->Add((TH1*)&h_met_enu_events_ejdr.GetValue());
        //lepjdr_Stack->Add((TH1*)&h_smu_munu_events_mujdr.GetValue());
        //lepjdr_Stack->Add((TH1*)&h_met_munu_events_mujdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_ww_enu_events_ejdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_wz_enu_events_ejdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_zz_enu_events_ejdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_ttZ_enu_events_ejdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_ww_munu_events_mujdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_wz_munu_events_mujdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_zz_munu_events_mujdr.GetValue());
        lepjdr_Stack->Add((TH1*)&h_ttZ_munu_events_mujdr.GetValue());


        auto h_events_lepjdr_canvas = new TCanvas("lepjdr", "lepjdr",10,10,900,900);

        h_d_enu_events_ejdr->GetXaxis()->SetTitle("delta R between jets and lepton/ rad");
        h_d_enu_events_ejdr->GetYaxis()->SetTitle("Events");

        legend_ed.AddEntry(h_d_enu_events_ejdr.GetPtr(),"tZq MC,electron channel deltaR between jets and electron","l");
        //legend_ebg.AddEntry(h_bg_enu_events_ejdr.GetPtr(),"tZq Background,electron channel deltaR between jets and electron","l");
        legend_md.AddEntry(h_d_munu_events_mujdr.GetPtr(),"tZq MC,muon channel deltaR between jets and muon","l");
        //legend_mbg.AddEntry(h_bg_munu_events_mujdr.GetPtr(),"tZq Background,muon channel deltaR between jets and muon","l");
        //legend_sed.AddEntry(h_se_enu_events_ejdr.GetPtr(),"Single electron, deltaR between jets and electron","l");
        //legend_emet.AddEntry(h_se_enu_events_ejdr.GetPtr(),"MET data,electron channel deltaR between jets and electron","l");
        //legend_smd.AddEntry(h_se_munu_events_mujdr.GetPtr(),"Signal muon,deltaR between jets and muon","l");
        //legend_mmet.AddEntry(h_met_munu_events_mujdr.GetPtr(),"MET data,muon channel deltaR between jets and muon","l");


        h_d_enu_events_ejdr->Draw();
        //h_bg_enu_events_ejdr->Draw("SAME");
        h_d_munu_events_mujdr->Draw("SAME");
        //h_bg_munu_events_mujdr->Draw("SAME");
        //h_se_enu_events_ejdr->Draw("SAME");
        //h_met_enu_events_ejdr->Draw("SAME");
        //h_smu_munu_events_mujdr->Draw("SAME");
        //h_met_munu_events_mujdr->Draw("SAME");
        h_ww_enu_events_ejdr->Draw("SAME");
        h_wz_enu_events_ejdr->Draw("SAME");
        h_zz_enu_events_ejdr->Draw("SAME");
        h_ttZ_enu_events_ejdr->Draw("SAME");
        h_ww_munu_events_mujdr->Draw("SAME");
        h_wz_munu_events_mujdr->Draw("SAME");
        h_zz_munu_events_mujdr->Draw("SAME");
        h_ttZ_munu_events_mujdr->Draw("SAME");

        h_events_lepjdr_canvas->cd(2);
        lepjdr_Stack->Draw("HIST");

	h_events_lepjdr_canvas->BuildLegend();
        h_events_lepjdr_canvas->SaveAs("hist_lepjdr.root");
        h_events_lepjdr_canvas->SaveAs("hist_lepjdr.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
        //legend_mbg.Clear();
 //       legend_sed.Clear();
 //       legend_emet.Clear();
 //       legend_smd.Clear();
 //       legend_mmet.Clear();


        auto h_d_enu_events_ezdr = d_enu_top_selection.Histo1D({"MC ezdr_enu_Channel","MC e and z deltaR in electron-neutrino channel",50,0,1},"z_e_min_dR","nw_z_e_min_dR");
        //auto h_bg_enu_events_ezdr = bg_enu_top_selection.Histo1D({"BG ezdr_enu_Channel","Background e and z deltaR in electron-neutrino channel",50,0,1},"z_e_min_dR","nw_z_e_min_dR");
        auto h_d_munu_events_muzdr = d_munu_top_selection.Histo1D({"MC muzdr_munu_Channel","MC mu and z deltaR in muon-neutrino channel",50,0,1},"z_mu_min_dR","nw_z_mu_min_dR");
        //auto h_bg_munu_events_muzdr = bg_munu_top_selection.Histo1D({"BG muzdr_munu_Channel","Background mu and z deltaR in muon-neutrino channel",50,0,1},"z_mu_min_dR","nw_z_mu_min_dR");
//        auto h_se_enu_events_ezdr = se_enu_top_selection.Histo1D({"Single Electron ezdr_enu_Channel","Single Electron e and z deltaR in electron-neutrino channel",50,0,1},"z_e_min_dR","nw_z_e_min_dR");
//        auto h_met_enu_events_ezdr = met_enu_top_selection.Histo1D({"MET ezdr_enu_Channel","MET e and z deltaR in electron-neutrino channel",50,0,1},"z_e_min_dR","nw_z_e_min_dR");
//        auto h_smu_munu_events_muzdr = sm_munu_top_selection.Histo1D({"Single Muon muzdr_munu_Channel","Single Muon mu and z deltaR in muon-neutrino channel",50,0,1},"z_mu_min_dR","nw_z_mu_min_dR");
//        auto h_met_munu_events_muzdr = met_munu_top_selection.Histo1D({"MET muzdr_munu_Channel","MET mu and z deltaR in muon-neutrino channel",50,0,1},"z_mu_min_dR","nw_z_mu_min_dR");
        auto h_ww_enu_events_ezdr = ww_enu_top_selection.Histo1D({"WW ezdr_enu_Channel","WW e and z deltaR in electron-neutrino channel",50,0,1},"z_e_min_dR","nw_z_e_min_dR");
        auto h_wz_enu_events_ezdr = wz_enu_top_selection.Histo1D({"WZ ezdr_enu_Channel","WZ e and z deltaR in electron-neutrino channel",50,0,1},"z_e_min_dR","nw_z_e_min_dR");
        auto h_ttZ_enu_events_ezdr = ttZ_enu_top_selection.Histo1D({"ttZ ezdr_enu_Channel","ttZ  e and z deltaR in electron-neutrino channel",50,0,1},"z_e_min_dR","nw_z_e_min_dR");
        auto h_zz_enu_events_ezdr = zz_enu_top_selection.Histo1D({"zz ezdr_enu_Channel","zz e and z deltaR in electron-neutrino channel",50,0,1},"z_e_min_dR","nw_z_e_min_dR");
        auto h_ww_munu_events_muzdr = ww_munu_top_selection.Histo1D({"WW muzdr_munu_Channel","WW mu and z deltaR in muon-neutrino channel",50,0,1},"z_mu_min_dR","nw_z_mu_min_dR");
        auto h_wz_munu_events_muzdr = wz_munu_top_selection.Histo1D({"WZ muzdr_munu_Channel","WZ mu and z deltaR in muon-neutrino channel",50,0,1},"z_mu_min_dR","nw_z_mu_min_dR");
        auto h_ttZ_munu_events_muzdr = ttZ_munu_top_selection.Histo1D({"ttZ muzdr_munu_Channel","ttZ  mu and z deltaR in muon-neutrino channel",50,0,1},"z_mu_min_dR","nw_z_mu_min_dR");
        auto h_zz_munu_events_muzdr = zz_munu_top_selection.Histo1D({"zz muzdr_munu_Channel","zz mu and z deltaR in muon-neutrino channel",50,0,1},"z_mu_min_dR","nw_z_mu_min_dR");

        THStack *lepzdr_Stack = new THStack("MC_Stack","Lepton & Z delta R ");
        h_d_enu_events_ezdr->SetLineColor(kBlack);
        //h_bg_enu_events_ezdr->SetLineColor(kRed);
        h_d_munu_events_muzdr->SetLineColor(kGreen);
        //h_bg_munu_events_muzdr->SetLineColor(kBlue);
//        h_se_enu_events_ezdr->SetLineColor(kPink);
//        h_met_enu_events_ezdr->SetLineColor(kCherry);
//        h_sm_munu_events_muzdr->SetLineColor(kViolet);
//        h_met_munu_events_muzdr->SetLineColor(kRose);
        h_ww_enu_events_ezdr->SetLineColor(kGray);
        h_wz_enu_events_ezdr->SetLineColor(kRed);
        h_zz_enu_events_ezdr->SetLineColor(kBlue);
        h_ttZ_enu_events_ezdr->SetLineColor(kCyan);
        h_ww_munu_events_muzdr->SetLineColor(kOrange);
        h_wz_munu_events_muzdr->SetLineColor(kSpring);
        h_zz_munu_events_muzdr->SetLineColor(kAzure);
        h_ttZ_munu_events_muzdr->SetLineColor(kTeal);


        lepzdr_Stack->Add((TH1*)&h_d_enu_events_ezdr.GetValue());
        //lepzdr_Stack->Add((TH1*)&h_bg_enu_events_ezdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_d_munu_events_muzdr.GetValue());
        //lepzdr_Stack->Add((TH1*)&h_bg_munu_events_muzdr.GetValue());
        //lepzdr_Stack->Add((TH1*)&h_se_enu_events_ezdr.GetValue());
        //lepzdr_Stack->Add((TH1*)&h_met_enu_events_ezdr.GetValue());
        //lepzdr_Stack->Add((TH1*)&h_smu_munu_events_muzdr.GetValue());
        //lepzdr_Stack->Add((TH1*)&h_met_munu_events_muzdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_ww_enu_events_ezdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_wz_enu_events_ezdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_zz_enu_events_ezdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_ttZ_enu_events_ezdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_ww_munu_events_muzdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_wz_munu_events_muzdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_zz_munu_events_muzdr.GetValue());
        lepzdr_Stack->Add((TH1*)&h_ttZ_munu_events_muzdr.GetValue());


        auto h_events_lepzdr_canvas = new TCanvas("lepzdr", "lepjdr",10,10,900,900);

        h_d_enu_events_ezdr->GetXaxis()->SetTitle("delta R between Z and lepton/ rad");
        h_d_enu_events_ezdr->GetYaxis()->SetTitle("Events");

        legend_ed.AddEntry(h_d_enu_events_ezdr.GetPtr(),"tZq MC,electron channel deltaR between Z and electron","l");
        //legend_ebg.AddEntry(h_bg_enu_events_ezdr.GetPtr(),"tZq Background,electron channel deltaR between Z and electron","l");
        legend_md.AddEntry(h_d_munu_events_muzdr.GetPtr(),"tZq MC,muon channel deltaR between Z and muon","l");
        //legend_mbg.AddEntry(h_bg_munu_events_muzdr.GetPtr(),"tZq Background,muon channel deltaR between Z and muon","l");
        //legend_sed.AddEntry(h_se_enu_events_ezdr.GetPtr(),"Single electron, deltaR between Z and electron","l");
        //legend_emet.AddEntry(h_se_enu_events_ezdr.GetPtr(),"MET data,electron channel deltaR between Z and electron","l");
        //legend_smd.AddEntry(h_se_munu_events_muzdr.GetPtr(),"Signal muon,deltaR between Z and muon","l");
        //legend_mmet.AddEntry(h_met_munu_events_muzdr.GetPtr(),"MET data,muon channel deltaR between Z and muon","l");


        h_d_enu_events_ezdr->Draw();
        //h_bg_enu_events_ezdr->Draw("SAME");
        h_d_munu_events_muzdr->Draw("SAME");
        //h_bg_munu_events_muzdr->Draw("SAME");
        //h_se_enu_events_ezdr->Draw("SAME");
        //h_met_enu_events_ezdr->Draw("SAME");
        //h_smu_munu_events_muzdr->Draw("SAME");
        //h_met_munu_events_muzdr->Draw("SAME");
        h_ww_enu_events_ezdr->Draw("SAME");
        h_wz_enu_events_ezdr->Draw("SAME");
        h_zz_enu_events_ezdr->Draw("SAME");
        h_ttZ_enu_events_ezdr->Draw("SAME");
        h_ww_munu_events_muzdr->Draw("SAME");
        h_wz_munu_events_muzdr->Draw("SAME");
        h_zz_munu_events_muzdr->Draw("SAME");
        h_ttZ_munu_events_muzdr->Draw("SAME");

        h_events_lepzdr_canvas->cd(2);
        lepzdr_Stack->Draw("HIST");

	h_events_lepzdr_canvas->BuildLegend();
        h_events_lepzdr_canvas->SaveAs("hist_lepzdr.root");
        h_events_lepzdr_canvas->SaveAs("hist_lepzdr.pdf");
        legend_ed.Clear();
        //legend_ebg.Clear();
        legend_md.Clear();
        //legend_mbg.Clear();
//        legend_sed.Clear();
//        legend_emet.Clear();
//        legend_smd.Clear();
//        legend_mmet.Clear();

////////////////////////////////////////////////////////////////////////////////// b jets ///////////////////////////////////////////////////////////////////////////

/*	auto h_d_enu_events_btag_pt = d_enu_top_selection.Histo1D({"MC Signal btag_pt_enu_Channel","MC b tag pt in electron-neutrino channel",50,0,400},"btag_numer_pt");
        auto h_d_enu_events_btag_eta = d_enu_top_selection.Histo1D({"MC btag_eta_enu_Channel","MC btag eta in electron-neutrino channel",50,-3,3}, "btag_numer_eta");
        auto h_d_enu_events_btag_numer_PtVsEta = d_enu_top_selection.Histo2D({"MC btag_Pt_vs_eta_enu_Channel","MC btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3},"btag_numer_pt", "btag_numer_eta"};
        auto h_d_enu_events_non_btag_numer_PtVsEta = d_enu_top_selection.Histo2D({"MC non btag_Pt_vs_eta_enu_Channel","MC non btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3},"non_btag_numer_pt","non_btag_numer_eta"};
        auto h_d_enu_events_btag_denom_PtVsEta = d_enu_top_selection.Histo2D({"MC btag_Pt_vs_eta_enu_Channel","MC btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3},"btag_denom_pt","btag_denom_eta"};
        auto h_d_enu_events_non_btag_denom_PtVsEta = d_enu_top_selection.Histo2D({"MC non btag_Pt_vs_eta_enu_Channel","MC non btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3},"non_btag_denom_pt","non_btag_denom_eta"};






	h_d_enu_events_btag_pt->SetLineColor(kBlack);

	auto h_events_btag_pt_canvas = new TCanvas("b tag pt", "b tag pt",10,10,900,900);

        h_d_enu_events_btag_pt->GetXaxis()->SetTitle("b tag of pt");
        h_d_enu_events_btag_pt->GetYaxis()->SetTitle("Events");
	h_d_enu_events_btag_pt->Draw();

	h_events_btag_pt_canvas->BuildLegend();
        h_events_btag_pt_canvas->SaveAs("hist_btag_pt.root");
        h_events_btag_pt_canvas->SaveAs("hist_btag_pt.pdf");
	legend_ed.Clear();


        h_d_enu_events_btag_eta->SetLineColor(kBlack);

        auto h_events_btag_eta_canvas = new TCanvas("b tag eta", "b tag eta",10,10,900,900);

        h_d_enu_events_btag_eta->GetXaxis()->SetTitle("b tag of eta");
        h_d_enu_events_btag_eta->GetYaxis()->SetTitle("Events");


        h_d_enu_events_btag_eta->Draw();

        h_events_btag_eta_canvas->BuildLegend();
        h_events_btag_eta_canvas->SaveAs("hist_btag_eta.root");
        h_events_btag_eta_canvas->SaveAs("hist_btag_eta.pdf");
        legend_ed.Clear();


        auto h_events_btag_numer_PtVsEta_canvas = new TCanvas("b tag pt Vs eta", "b tag pt Vs eta",10,10,900,900);

        h_d_enu_events_btag_numer_PtVsEta->GetXaxis()->SetTitle("b tag numer pt");
        h_d_enu_events_btag_numer_PtVsEta->GetYaxis()->SetTitle("b tag numer eta");

        h_events_btag_numer_PtVsEta_canvas->BuildLegend();
        h_d_enu_events_btag_numer_PtVsEta->Draw("COLZ");

        h_events_btag_numer_PtVsEta_canvas->SaveAs("hist_btag_numer_PtVsEta.root");
        h_events_btag_numer_PtVsEta_canvas->SaveAs("hist_btag_numer_PtVsEta.pdf");


        auto h_events_non_btag_numer_PtVsEta_canvas = new TCanvas("non b tag pt Vs eta", "non b tag pt Vs eta",10,10,900,900);

        h_d_enu_events_non_btag_numer_PtVsEta->GetXaxis()->SetTitle("non b tag numer pt");
        h_d_enu_events_non_btag_numer_PtVsEta->GetYaxis()->SetTitle("non b tag numer eta");

        h_events_non_btag_numer_PtVsEta_canvas->BuildLegend();
        h_d_enu_events_non_btag_numer_PtVsEta->Draw("COLZ");

        h_events_non_btag_numer_PtVsEta_canvas->SaveAs("hist_non_btag_numer_PtVsEta.root");
        h_events_non_btag_numer_PtVsEta_canvas->SaveAs("hist_non_btag_numer_PtVsEta.pdf");


	auto h_d_enu_events_btag_denom_PtVsEta = d_enu_top_selection.Histo2D({"MC btag_Pt_vs_eta_enu_Channel","MC btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3}, "btag_denom_pt" , "btag_denom_eta");

        auto h_events_btag_denom_PtVsEta_canvas = new TCanvas("b tag pt Vs eta", "b tag pt Vs eta",10,10,900,900);

        h_d_enu_events_btag_denom_PtVsEta->GetXaxis()->SetTitle("b tag denom pt");
        h_d_enu_events_btag_denom_PtVsEta->GetYaxis()->SetTitle("b tag denom eta");

        h_events_btag_denom_PtVsEta_canvas->BuildLegend();
        h_d_enu_events_btag_denom_PtVsEta->Draw("COLZ");

        h_events_btag_denom_PtVsEta_canvas->SaveAs("hist_btag_denom_PtVsEta.root");
        h_events_btag_denom_PtVsEta_canvas->SaveAs("hist_btag_denom_PtVsEta.pdf");


	auto h_d_enu_events_non_btag_denom_PtVsEta = d_enu_top_selection.Histo2D({"MC non btag_Pt_vs_eta_enu_Channel","MC non btag pt Vs eta in electron-neutrino channel",50,0,400,50,-3,3}, "non_btag_denom_pt" , "non_btag_denom_eta");

        auto h_events_non_btag_denom_PtVsEta_canvas = new TCanvas("non b tag pt Vs eta", "non b tag pt Vs eta",10,10,900,900);

        h_d_enu_events_non_btag_denom_PtVsEta->GetXaxis()->SetTitle("non b tag denom pt");
        h_d_enu_events_non_btag_denom_PtVsEta->GetYaxis()->SetTitle("non b tag denom eta");

        h_events_non_btag_denom_PtVsEta_canvas->BuildLegend();
        h_d_enu_events_non_btag_denom_PtVsEta->Draw("COLZ");

        h_events_non_btag_denom_PtVsEta_canvas->SaveAs("hist_non_btag_denom_PtVsEta.root");
        h_events_non_btag_denom_PtVsEta_canvas->SaveAs("hist_non_btag_denom_PtVsEta.pdf");


	auto h_events_btag_PtVsEta_canvas = new TCanvas("b tag pt Vs eta", "b tag pt Vs eta",10,10,900,900);
	TH2D *btag_ratio = new TH2D("ei","b tag ei",50,0,400,50,-3,3);
	btag_ratio = (TH2D*)h_d_enu_events_btag_numer_PtVsEta->Clone();
	btag_ratio->GetXaxis()->SetTitle(" b tag Pt");
	btag_ratio->GetYaxis()->SetTitle("b tag eta");
	btag_ratio->Divide(h_d_enu_events_btag_denom_PtVsEta.GetPtr());
	h_events_btag_PtVsEta_canvas->BuildLegend();
	btag_ratio->Draw("COLZ");
	h_events_btag_PtVsEta_canvas->SaveAs("h_events_btag_PtVsEta_canvas.root");
        h_events_btag_PtVsEta_canvas->SaveAs("h_events_btag_PtVsEta_canvas.pdf");

        auto h_events_non_btag_PtVsEta_canvas = new TCanvas("non b tag pt Vs eta", "non b tag pt Vs eta",10,10,900,900);
        TH2D *non_btag_ratio = new TH2D("ej","non b tag ei",50,0,400,50,-3,3);
        non_btag_ratio = (TH2D*)h_d_enu_events_non_btag_numer_PtVsEta->Clone();
        non_btag_ratio->GetXaxis()->SetTitle("non b tag Pt");
        non_btag_ratio->GetYaxis()->SetTitle("non b tag eta");
        non_btag_ratio->Divide(h_d_enu_events_non_btag_denom_PtVsEta.GetPtr());
        h_events_non_btag_PtVsEta_canvas->BuildLegend();
        non_btag_ratio->Draw("COLZ");
        h_events_non_btag_PtVsEta_canvas->SaveAs("h_events_non_btag_PtVsEta_canvas.root");
        h_events_non_btag_PtVsEta_canvas->SaveAs("h_events_non_btag_PtVsEta_canvas.pdf");


*/

/////////////////////////////////////////////////////////////////////// E-Nu Channel histograms AND Canvases /////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////PT//////////////////////////////////////////////////////////////////////////////////
/*        auto h_d_enu_events_ept = d_enu_top_selection.Histo1D({"MC electron_pt_enu_Channel","MC electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
/*        auto h_ww_enu_events_ept = ww_enu_top_selection.Histo1D({"WW electron_pt_enu_Channel","WW electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
        auto h_wz_enu_events_ept = wz_enu_top_selection.Histo1D({"WZ electron_pt_enu_Channel","WZ electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
        auto h_zz_enu_events_ept = zz_enu_top_selection.Histo1D({"ZZ electron_pt_enu_Channel","ZZ electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
        auto h_ttZ_enu_events_ept = ttZ_enu_top_selection.Histo1D({"ttZ electron_pt_enu_Channel","ttZ electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
	auto h_bg_enu_events_ept = bg_enu_top_selection.Histo1D({"bg electron_pt_enu_Channel","Back ground electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
//	auto h_se_enu_events_ept = se_enu_top_selection.Histo1D({"Single Electron electron_pt_enu_Channel","Single Electron electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
//	auto h_met_enu_events_ept = met_enu_top_selection.Histo1D({"MET electron_pt_enu_Channel","MET electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");

        auto h_events_ept_canvas = new TCanvas("electron pt", "electron pt",10,10,900,900);

	h_d_enu_events_ept->GetXaxis()->SetTitle("Pt/GeV");
        h_d_enu_events_ept->GetYaxis()->SetTitle("Events");

        h_d_enu_events_ept->SetLineColor(kBlack);
	h_bg_enu_events_ept->SetLineColor(kPink-6);
/*        h_ww_enu_events_ept->SetLineColor(kRed);
        h_wz_enu_events_ept->SetLineColor(kOrange);
        h_ttZ_enu_events_ept->SetLineColor(kYellow);
        h_zz_enu_events_ept->SetLineColor(kTeal);
//	h_se_enu_events_ept->SetLineColor(kViolet);
//	h_met_enu_events_ept->SetLineColor(kPink);

	h_d_enu_events_ept->Scale(NWS_F);
/*	h_ww_enu_events_ept->Scale(NWS_F);
	h_wz_enu_events_ept->Scale(NWS_F);
	h_ttZ_enu_events_ept->Scale(NWS_F);
	h_zz_enu_events_ept->Scale(NWS_F);
	h_bg_enu_events_ept->Scale(NWS_F);


	h_d_enu_events_ept->Draw();
/*        h_ww_enu_events_ept->Draw("SAME");
        h_wz_enu_events_ept->Draw("SAME");
        h_ttZ_enu_events_ept->Draw("SAME");
        h_zz_enu_events_ept->Draw("SAME");
	h_bg_enu_events_ept->Draw("SAME");
//	h_se_enu_events_ept->Draw("SAME");
//	h_met_enu_events_ept->Draw("SAME");

        h_events_ept_canvas->BuildLegend();
	h_events_ept_canvas->SaveAs("enu_pt.root");
	h_events_ept_canvas->SaveAs("enu_pt.pdf");


	auto h_d_enu_events_jpt = d_enu_top_selection.Histo1D({"MC jet_pt_enu_Channel","MC jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
/*        auto h_ww_enu_events_jpt = ww_enu_top_selection.Histo1D({"WW jet_pt_enu_Channel","WW jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_wz_enu_events_jpt = wz_enu_top_selection.Histo1D({"WZ jet_pt_enu_Channel","WZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_zz_enu_events_jpt = zz_enu_top_selection.Histo1D({"ZZ jet_pt_enu_Channel","ZZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_ttZ_enu_events_jpt = ttZ_enu_top_selection.Histo1D({"ttZ jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
	auto h_bg_enu_events_jpt = bg_enu_top_selection.Histo1D({"bg jet_pt_enu_Channel","Back ground jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
//	auto h_se_enu_events_jpt = se_enu_top_selection.Histo1D({"Single Electron jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
//	auto h_met_enu_events_jpt = met_enu_top_selection.Histo1D({"MET jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");

	auto h_events_jetpt_canvas = new TCanvas("e nu jet pt ", "enu jet pt",10,10,900,900);

        h_d_enu_events_jpt->GetXaxis()->SetTitle("pt / GeV");
        h_d_enu_events_jpt->GetYaxis()->SetTitle("Events");

        h_d_enu_events_jpt->SetLineColor(kBlack);
	h_bg_enu_events_jpt->SetLineColor(kPink-6);
/*        h_ww_enu_events_jpt->SetLineColor(kRed);
        h_wz_enu_events_jpt->SetLineColor(kOrange);
        h_ttZ_enu_events_jpt->SetLineColor(kYellow);
        h_zz_enu_events_jpt->SetLineColor(kTeal);
//        h_se_enu_events_jpt->SetLineColor(kViolet);
//        h_met_enu_events_jpt->SetLineColor(kPink);

        h_d_enu_events_jpt->Scale(NWS_F);
/*        h_ww_enu_events_jpt->Scale(NWS_F);
        h_wz_enu_events_jpt->Scale(NWS_F);
        h_ttZ_enu_events_jpt->Scale(NWS_F);
        h_zz_enu_events_jpt->Scale(NWS_F);
	h_bg_enu_events_jpt->Scale(NWS_F);

        h_d_enu_events_jpt->Draw();
/*        h_ww_enu_events_jpt->Draw("SAME");
        h_wz_enu_events_jpt->Draw("SAME");
        h_ttZ_enu_events_jpt->Draw("SAME");
        h_zz_enu_events_jpt->Draw("SAME");
	h_bg_enu_events_jpt->Draw("SAME");
//        h_se_enu_events_jpt->Draw("SAME");
//        h_met_enu_events_jpt->Draw("SAME");

        h_events_jetpt_canvas->BuildLegend();
        h_events_jetpt_canvas->SaveAs("enu_jetpt.root");
        h_events_jetpt_canvas->SaveAs("enu_jetpt.pdf");




        auto h_d_enu_events_eeta = d_enu_top_selection.Histo1D({"MC electron_eta_enu_Channel","MC electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
	auto h_bg_enu_events_eeta = bg_enu_top_selection.Histo1D({"bg electron_eta_enu_Channel","Back ground electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
/*        auto h_ww_enu_events_eeta = ww_enu_top_selection.Histo1D({"WW electron_eta_enu_Channel","WW electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
        auto h_wz_enu_events_eeta = wz_enu_top_selection.Histo1D({"WZ electron_eta_enu_Channel","WZ electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
        auto h_zz_enu_events_eeta = zz_enu_top_selection.Histo1D({"ZZ electron_eta_enu_Channel","ZZ electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
        auto h_ttZ_enu_events_eeta = ttZ_enu_top_selection.Histo1D({"ttZ electron_eta_enu_Channel","ttZ electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
//        auto h_se_enu_events_eeta = se_enu_top_selection.Histo1D({"Single Electron electron_eta_enu_Channel","Single Electron electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
//        auto h_met_enu_events_eeta = met_enu_top_selection.Histo1D({"MET electron_eta_enu_Channel","MET electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");


        auto h_events_eeta_canvas = new TCanvas("electron eta", "electron eta",10,10,900,900);

        h_d_enu_events_eeta->GetXaxis()->SetTitle("eta");
        h_d_enu_events_eeta->GetYaxis()->SetTitle("Events");

        h_d_enu_events_eeta->SetLineColor(kBlack);
	h_bg_enu_events_eeta->SetLineColor(kPink-6);
/*        h_ww_enu_events_eeta->SetLineColor(kRed);
        h_wz_enu_events_eeta->SetLineColor(kOrange);
        h_ttZ_enu_events_eeta->SetLineColor(kYellow);
        h_zz_enu_events_eeta->SetLineColor(kTeal);
//        h_se_enu_events_eeta->SetLineColor(kViolet);
//        h_met_enu_events_eeta->SetLineColor(kPink);

        h_d_enu_events_eeta->Scale(NWS_F);
/*        h_ww_enu_events_eeta->Scale(NWS_F);
        h_wz_enu_events_eeta->Scale(NWS_F);
        h_ttZ_enu_events_eeta->Scale(NWS_F);
        h_zz_enu_events_eeta->Scale(NWS_F);
	h_bg_enu_events_eeta->Scale(NWS_F);


        h_d_enu_events_eeta->Draw();
/*        h_ww_enu_events_eeta->Draw("SAME");
        h_wz_enu_events_eeta->Draw("SAME");
        h_ttZ_enu_events_eeta->Draw("SAME");
        h_zz_enu_events_eeta->Draw("SAME");
	h_bg_enu_events_eeta->Draw("SAME");
//        h_se_enu_events_eeta->Draw("SAME");
//        h_met_enu_events_eeta->Draw("SAME");

        h_events_eeta_canvas->BuildLegend();
        h_events_eeta_canvas->SaveAs("enu_eta.root");
        h_events_eeta_canvas->SaveAs("enu_eta.pdf");



        auto h_d_enu_events_jeta = d_enu_top_selection.Histo1D({"MC jet_eta_enu_Channel","MC jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
	auto h_bg_enu_events_jeta = bg_enu_top_selection.Histo1D({"bg jet_eta_enu_Channel","Back ground jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
/*        auto h_ww_enu_events_jeta = ww_enu_top_selection.Histo1D({"WW jet_eta_enu_Channel","WW jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_wz_enu_events_jeta = wz_enu_top_selection.Histo1D({"WZ jet_eta_enu_Channel","WZ jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_zz_enu_events_jeta = zz_enu_top_selection.Histo1D({"ZZ jet_eta_enu_Channel","ZZ jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_ttZ_enu_events_jeta = ttZ_enu_top_selection.Histo1D({"ttZ jet_eta_enu_Channel","ttZ jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
//        auto h_se_enu_events_jeta = se_enu_top_selection.Histo1D({"Single Electron jet_eta_enu_Channel","Single Electron jet eta in electron-neutrino channels",50,-4,4},"tight_jets_eta");
//        auto h_met_enu_events_jeta = met_enu_top_selection.Histo1D({"MET electron_jet_enu_Channel","MET jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");


        auto h_events_jeta_canvas = new TCanvas("jet eta", "jet eta",10,10,900,900);

        h_d_enu_events_jeta->GetXaxis()->SetTitle("eta");
        h_d_enu_events_jeta->GetYaxis()->SetTitle("Events");

        h_d_enu_events_jeta->SetLineColor(kBlack);
	h_bg_enu_events_jeta->SetLineColor(kPink-6);
/*        h_ww_enu_events_jeta->SetLineColor(kRed);
        h_wz_enu_events_jeta->SetLineColor(kOrange);
        h_ttZ_enu_events_jeta->SetLineColor(kYellow);
        h_zz_enu_events_jeta->SetLineColor(kTeal);
//        h_se_enu_events_jeta->SetLineColor(kViolet);
//        h_met_enu_events_jeta->SetLineColor(kPink);


        h_d_enu_events_jeta->Scale(NWS_F);
/*        h_ww_enu_events_jeta->Scale(NWS_F);
        h_wz_enu_events_jeta->Scale(NWS_F);
        h_ttZ_enu_events_jeta->Scale(NWS_F);
        h_zz_enu_events_jeta->Scale(NWS_F);
	h_bg_enu_events_jeta->Scale(NWS_F);

        h_d_enu_events_jeta->Draw();
/*        h_ww_enu_events_jeta->Draw("SAME");
        h_wz_enu_events_jeta->Draw("SAME");
        h_ttZ_enu_events_jeta->Draw("SAME");
        h_zz_enu_events_jeta->Draw("SAME");
	h_bg_enu_events_jeta->Draw("SAME");
//        h_se_enu_events_jeta->Draw("SAME");
//        h_met_enu_events_jeta->Draw("SAME");

        h_events_jeta_canvas->BuildLegend();
        h_events_jeta_canvas->SaveAs("enu_jeta.root");
        h_events_jeta_canvas->SaveAs("enu_jeta.pdf");



        auto h_d_enu_events_wmass = d_enu_top_selection.Histo1D({"MC enu_w_mass","MC electron-neutrino transverse w mass",50,0,500},"w_e_mass");
/*        auto h_ww_enu_events_wmass = ww_enu_top_selection.Histo1D({"WW enu_w_mass","WW electron-neutrino transverse w mass",50,0,500},"w_e_mass");
        auto h_wz_enu_events_wmass = wz_enu_top_selection.Histo1D({"WZ enu_w_mass","WZ electron-neutrino transverse w mass",50,0,500},"w_e_mass");
        auto h_ttZ_enu_events_wmass = ttZ_enu_top_selection.Histo1D({"ttZ enu_w_mass","ttZ electron-neutrino transverse w mass",50,0,500},"w_e_mass");
        auto h_zz_enu_events_wmass = zz_enu_top_selection.Histo1D({"ZZ enu_w_mass","ZZ electron-neutrino transverse w mass",50,0,500},"w_e_mass");
	auto h_bg_enu_events_wmass = bg_enu_top_selection.Histo1D({"BG enu_w_mass","Back ground electron-neutrino transverse w mass",50,0,500},"w_e_mass");
//        auto h_se_enu_events_wmass = se_enu_top_selection.Histo1D({"Single Electron enu_w_mass","Single Electron electron-neutrino transverse w mass",50,0,500},"w_e_mass");
//        auto h_met_enu_events_wmass = met_enu_top_selection.Histo1D({"MET enu_w_mass","MET electron-neutrino transverse w mass",50,0,500},"w_e_mass");


        auto h_events_wmass_canvas = new TCanvas("enu w mass ", "enu w mass",10,10,900,900);

        h_d_enu_events_wmass->GetXaxis()->SetTitle("mass/GeV/C^2");
        h_d_enu_events_wmass->GetYaxis()->SetTitle("Events");

        h_d_enu_events_wmass->SetLineColor(kBlack);
	h_bg_enu_events_wmass->SetLineColor(kPink-6);
/*        h_ww_enu_events_wmass->SetLineColor(kRed);
        h_wz_enu_events_wmass->SetLineColor(kOrange);
        h_ttZ_enu_events_wmass->SetLineColor(kYellow);
        h_zz_enu_events_wmass->SetLineColor(kTeal);
//        h_se_enu_events_wmass->SetLineColor(kViolet);
//        h_met_enu_events_wmass->SetLineColor(kPink);

        h_d_enu_events_wmass->Scale(NWS_F);
/*        h_ww_enu_events_wmass->Scale(NWS_F);
        h_wz_enu_events_wmass->Scale(NWS_F);
        h_ttZ_enu_events_wmass->Scale(NWS_F);
        h_zz_enu_events_wmass->Scale(NWS_F);
	h_bg_enu_events_wmass->Scale(NWS_F);



        h_d_enu_events_wmass->Draw();
/*        h_ww_enu_events_wmass->Draw("SAME");
        h_wz_enu_events_wmass->Draw("SAME");
        h_ttZ_enu_events_wmass->Draw("SAME");
        h_zz_enu_events_wmass->Draw("SAME");
	h_bg_enu_events_wmass->Draw("SAME");
//        h_se_enu_events_wmass->Draw("SAME");
//        h_met_enu_events_wmass->Draw("SAME");

        h_events_wmass_canvas->BuildLegend();
        h_events_wmass_canvas->SaveAs("enu_transverse_Wmass.root");
        h_events_wmass_canvas->SaveAs("enu_transverse_Wmass.pdf");




        auto h_d_enu_events_zmass = d_enu_top_selection.Histo1D({"MC Z_mass_enu_Channel","MC Z mass in electron-neutrino channel",50,0,500},"z_mass");
	auto h_bg_enu_events_zmass = bg_enu_top_selection.Histo1D({"BG Z_mass_enu_Channel","Background Z mass in electron-neutrino channel",50,0,500},"z_mass");
/*        auto h_ww_enu_events_zmass = ww_enu_top_selection.Histo1D({"WW Z_mass_enu_Channel","WW Z mass in electron-neutrino channel",50,0,500},"z_mass");
        auto h_wz_enu_events_zmass = wz_enu_top_selection.Histo1D({"WZ Z_mass_enu_Channel","WZ Z mass in electron-neutrino channel",50,0,500},"z_mass");
        auto h_ttZ_enu_events_zmass = ttZ_enu_top_selection.Histo1D({"ttZ Z_mass_enu_Channel","ttZ Z mass in electron-neutrino channel",50,0,500},"z_mass");
        auto h_zz_enu_events_zmass = zz_enu_top_selection.Histo1D({"zz Z_mass_enu_Channel","zz Z mass in electron-neutrino channel",50,0,500},"z_mass");
//	auto h_se_enu_events_zmass = se_enu_top_selection.Histo1D({"Single Electron Z_mass_enu_Channel","Single Electron Z mass in electron-neutrino channel",50,0,500},"z_mass");
//	auto h_met_enu_events_zmass = met_enu_top_selection.Histo1D({"MET Z_mass_enu_Channel","MET Z mass in electron-neutrino channel",50,0,500},"z_mass");

        auto h_events_zmass_canvas = new TCanvas("Z mass", "Z mass",10,10,900,900);

        h_d_enu_events_zmass->GetXaxis()->SetTitle("mass/GeVC^2");
        h_d_enu_events_zmass->GetYaxis()->SetTitle("Events");

        h_d_enu_events_zmass->SetLineColor(kBlack);
	h_bg_enu_events_zmass->SetLineColor(kPink-6);
/*        h_ww_enu_events_zmass->SetLineColor(kRed);
        h_wz_enu_events_zmass->SetLineColor(kOrange);
        h_ttZ_enu_events_zmass->SetLineColor(kYellow);
        h_zz_enu_events_zmass->SetLineColor(kTeal);
//	h_se_enu_events_zmass->SetLineColor(kViolet);
//	h_met_enu_events_zmass->SetLineColor(kPink);

	h_d_enu_events_zmass->Scale(NWS_F);
/*        h_ww_enu_events_zmass->Scale(NWS_F);
        h_wz_enu_events_zmass->Scale(NWS_F);
        h_ttZ_enu_events_zmass->Scale(NWS_F);
        h_zz_enu_events_zmass->Scale(NWS_F);
	h_bg_enu_events_zmass->Scale(NWS_F);


	h_d_enu_events_zmass->Draw();
/*        h_ww_enu_events_zmass->Draw("SAME");
        h_wz_enu_events_zmass->Draw("SAME");
        h_ttZ_enu_events_zmass->Draw("SAME");
        h_zz_enu_events_zmass->Draw("SAME");
	h_bg_enu_events_zmass->Draw("SAME");
//	h_se_enu_events_zmass->Draw("SAME");
//	h_met_enu_events_zmass->Draw("SAME");

	h_events_zmass_canvas->BuildLegend();
        h_events_zmass_canvas->SaveAs("en_Z_mass.root");
	h_events_zmass_canvas->SaveAs("en_Z_mass.pdf");



        auto h_d_enu_events_jdphi = d_enu_top_selection.Histo1D({"MC jdphi_enu_Channel","MC jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
	auto h_bg_enu_events_jdphi = bg_enu_top_selection.Histo1D({"Background jdphi_enu_Channel","Back ground jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
/*        auto h_ww_enu_events_jdphi = ww_enu_top_selection.Histo1D({"WW jdphi_enu_Channel","WW jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_wz_enu_events_jdphi = wz_enu_top_selection.Histo1D({"WZ jdphi_enu_Channel","WZ jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_ttZ_enu_events_jdphi = ttZ_enu_top_selection.Histo1D({"ttZ jdphi_enu_Channel","ttZ jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_zz_enu_events_jdphi = zz_enu_top_selection.Histo1D({"zz jdphi_enu_Channel","zz jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
//        auto h_se_enu_events_jdphi = se_enu_top_selection.Histo1D({"Single Electron jdphi_enu_Channel","Single Electron jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
//        auto h_met_enu_events_jdphi = met_enu_top_selection.Histo1D({"MET jdphi_enu_Channel","MET jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");


        auto h_events_jdphi_canvas = new TCanvas("jetdeltaphi", "jetdeltaphi",10,10,900,900);

        h_d_enu_events_jdphi->GetXaxis()->SetTitle("jets delta phi / rad");
        h_d_enu_events_jdphi->GetYaxis()->SetTitle("Events");

        h_d_enu_events_jdphi->SetLineColor(kBlack);
	h_bg_enu_events_jdphi->SetLineColor(kPink-6);
/*        h_ww_enu_events_jdphi->SetLineColor(kRed);
        h_wz_enu_events_jdphi->SetLineColor(kOrange);
        h_ttZ_enu_events_jdphi->SetLineColor(kYellow);
        h_zz_enu_events_jdphi->SetLineColor(kTeal);
//        h_se_enu_events_jdphi->SetLineColor(kViolet);
//        h_met_enu_events_jdphi->SetLineColor(kPink);

        h_d_enu_events_jdphi->Scale(NWS_F);
/*        h_ww_enu_events_jdphi->Scale(NWS_F);
        h_wz_enu_events_jdphi->Scale(NWS_F);
        h_ttZ_enu_events_jdphi->Scale(NWS_F);
        h_zz_enu_events_jdphi->Scale(NWS_F);
	h_bg_enu_events_jdphi->Scale(NWS_F);

        h_d_enu_events_jdphi->Draw();
/*        h_ww_enu_events_jdphi->Draw("SAME");
        h_wz_enu_events_jdphi->Draw("SAME");
        h_ttZ_enu_events_jdphi->Draw("SAME");
        h_zz_enu_events_jdphi->Draw("SAME");
	h_bg_enu_events_jdphi->Draw("SAME");
//        h_se_enu_events_jdphi->Draw("SAME");
//        h_met_enu_events_jdphi->Draw("SAME");

        h_events_jdphi_canvas->BuildLegend();
        h_events_jdphi_canvas->SaveAs("en_jet_dphi.root");
        h_events_jdphi_canvas->SaveAs("en_jet_dphi.pdf");



	auto h_d_enu_events_zmetdphi = d_enu_top_selection.Histo1D({"MC zmetdphi_enu_Channel","MC Z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_bg_enu_events_zmetdphi = bg_enu_top_selection.Histo1D({"Background zmetdphi_enu_Channel","Background Z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
/*	auto h_ww_enu_events_zmetdphi = ww_enu_top_selection.Histo1D({"WW zmetdphi_enu_Channel","WW Z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_wz_enu_events_zmetdphi = wz_enu_top_selection.Histo1D({"WZ zmetdphi_enu_Channel","WZ Z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_ttZ_enu_events_zmetdphi = ttZ_enu_top_selection.Histo1D({"ttZ zmetdphi_enu_Channel","ttZ Z met pt jet deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_zz_enu_events_zmetdphi = zz_enu_top_selection.Histo1D({"zz zmetdphi_enu_Channel","zz z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
//        auto h_se_enu_events_zmetdphi = se_enu_top_selection.Histo1D({"Single Electron zmetdphi_enu_Channel","Single Electron z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
//       auto h_met_enu_events_zmetdphi = met_enu_top_selection.Histo1D({"MET zmetdphi_enu_Channel","MET z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");



        auto h_events_zmetdphi_canvas = new TCanvas("zmetdeltaphi", "zmetdeltaphi",10,10,900,900);

        h_d_enu_events_zmetdphi->GetXaxis()->SetTitle("z and met delta phi / rad");
        h_d_enu_events_zmetdphi->GetYaxis()->SetTitle("Events");

        h_d_enu_events_zmetdphi->SetLineColor(kBlack);
	h_bg_enu_events_zmetdphi->SetLineColor(kPink-6);
/*        h_ww_enu_events_zmetdphi->SetLineColor(kRed);
        h_wz_enu_events_zmetdphi->SetLineColor(kOrange);
        h_ttZ_enu_events_zmetdphi->SetLineColor(kYellow);
        h_zz_enu_events_zmetdphi->SetLineColor(kTeal);
//        h_se_enu_events_zmetdphi->SetLineColor(kViolet);
//        h_met_enu_events_zmetdphi->SetLineColor(kPink);

        h_d_enu_events_zmetdphi->Scale(NWS_F);
/*        h_ww_enu_events_zmetdphi->Scale(NWS_F);
        h_wz_enu_events_zmetdphi->Scale(NWS_F);
        h_ttZ_enu_events_zmetdphi->Scale(NWS_F);
        h_zz_enu_events_zmetdphi->Scale(NWS_F);
	h_bg_enu_events_zmetdphi->Scale(NWS_F);

        h_d_enu_events_zmetdphi->Draw();
/*        h_ww_enu_events_zmetdphi->Draw("SAME");
        h_wz_enu_events_zmetdphi->Draw("SAME");
        h_ttZ_enu_events_zmetdphi->Draw("SAME");
        h_zz_enu_events_zmetdphi->Draw("SAME");
	h_bg_enu_events_zmetdphi->Draw("SAME");
//        h_se_enu_events_zmetdphi->Draw("SAME");
//        h_met_enu_events_zmetdphi->Draw("SAME");

        h_events_zmetdphi_canvas->BuildLegend();
        h_events_zmetdphi_canvas->SaveAs("en_zmetpt_dphi.root");
        h_events_zmetdphi_canvas->SaveAs("en_zmetpt_dphi.pdf");


        auto h_d_enu_events_zwdphi = d_enu_top_selection.Histo1D({"MC zwdphi_enu_Channel","MC Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
	auto h_bg_enu_events_zwdphi = bg_enu_top_selection.Histo1D({"bg zwdphi_enu_Channel","Background Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
/*        auto h_ww_enu_events_zwdphi = ww_enu_top_selection.Histo1D({"WW zwdphi_enu_Channel","WW Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_wz_enu_events_zwdphi = wz_enu_top_selection.Histo1D({"WZ zwdphi_enu_Channel","WZ Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_ttZ_enu_events_zwdphi = ttZ_enu_top_selection.Histo1D({"ttZ zwdphi_enu_Channel","ttZ  Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_zz_enu_events_zwdphi = zz_enu_top_selection.Histo1D({"zz zwdphi_enu_Channel","zz z Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
//        auto h_se_enu_events_zwdphi = se_enu_top_selection.Histo1D({"Single Electron zmetdphi_enu_Channel","Single Electron z met pt deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
//        auto h_met_enu_events_zwdphi = met_enu_top_selection.Histo1D({"MET zmetdphi_enu_Channel","MET z met pt deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");


        auto h_events_zwdphi_canvas = new TCanvas("zwdeltaphi", "zwdeltaphi",10,10,900,900);

        h_d_enu_events_zwdphi->GetXaxis()->SetTitle("z and w delta phi / rad");
        h_d_enu_events_zwdphi->GetYaxis()->SetTitle("Events");

        h_d_enu_events_zwdphi->SetLineColor(kBlack);
	h_bg_enu_events_zwdphi->SetLineColor(kPink-6);
/*        h_ww_enu_events_zwdphi->SetLineColor(kRed);
        h_wz_enu_events_zwdphi->SetLineColor(kOrange);
        h_ttZ_enu_events_zwdphi->SetLineColor(kYellow);
        h_zz_enu_events_zwdphi->SetLineColor(kTeal);
//        h_se_enu_events_zwdphi->SetLineColor(kViolet);
//        h_met_enu_events_zwdphi->SetLineColor(kPink);

        h_d_enu_events_zwdphi->Scale(NWS_F);
/*        h_ww_enu_events_zwdphi->Scale(NWS_F);
        h_wz_enu_events_zwdphi->Scale(NWS_F);
        h_ttZ_enu_events_zwdphi->Scale(NWS_F);
        h_zz_enu_events_zwdphi->Scale(NWS_F);
	h_bg_enu_events_zwdphi->Scale(NWS_F);

        h_d_enu_events_zwdphi->Draw();
/*        h_ww_enu_events_zwdphi->Draw("SAME");
        h_wz_enu_events_zwdphi->Draw("SAME");
        h_ttZ_enu_events_zwdphi->Draw("SAME");
        h_zz_enu_events_zmetdphi->Draw("SAME");
	h_bg_enu_events_zmetdphi->Draw("SAME");
//        h_se_enu_events_zwdphi->Draw("SAME");
//        h_met_enu_events_zwdphi->Draw("SAME");

        h_events_zwdphi_canvas->BuildLegend();
        h_events_zwdphi_canvas->SaveAs("en_zwd_dphi.root");
        h_events_zwdphi_canvas->SaveAs("en_zw_dphi.pdf");




        auto h_d_enu_events_ejdr = d_enu_top_selection.Histo1D({"MC ejdr_enu_Channel","MC e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
        auto h_bg_enu_events_ejdr = bg_enu_top_selection.Histo1D({"bg ejdr_enu_Channel","Background e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
/*	auto h_ww_enu_events_ejdr = ww_enu_top_selection.Histo1D({"WW ejdr_enu_Channel","WW e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
        auto h_wz_enu_events_ejdr = wz_enu_top_selection.Histo1D({"WZ ejdr_enu_Channel","WZ e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
        auto h_ttZ_enu_events_ejdr = ttZ_enu_top_selection.Histo1D({"ttZ ejdr_enu_Channel","ttZ  e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
        auto h_zz_enu_events_ejdr = zz_enu_top_selection.Histo1D({"zz ejdr_enu_Channel","zz e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
//        auto h_se_enu_events_ejdr = se_enu_top_selection.Histo1D({"Single Electron ejdr_enu_Channel","Single Electron e and jet pt deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
//        auto h_met_enu_events_ejdr = met_enu_top_selection.Histo1D({"MET ejdr_enu_Channel","MET e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");

        auto h_events_ejdr_canvas = new TCanvas("ejdeltar", "ejdeltar",10,10,900,900);

        h_d_enu_events_ejdr->GetXaxis()->SetTitle("e and jet deltaR");
        h_d_enu_events_ejdr->GetYaxis()->SetTitle("Events");



        h_d_enu_events_ejdr->SetLineColor(kBlack);
	h_bg_enu_events_ejdr->SetLineColor(kPink-6);
/*        h_ww_enu_events_ejdr->SetLineColor(kRed);
        h_wz_enu_events_ejdr->SetLineColor(kOrange);
        h_ttZ_enu_events_ejdr->SetLineColor(kYellow);
        h_zz_enu_events_ejdr->SetLineColor(kTeal);
//        h_se_enu_events_ejdr->SetLineColor(kViolet);
//        h_met_enu_events_ejdr->SetLineColor(kPink);

        h_d_enu_events_ejdr->Scale(NWS_F);
/*        h_ww_enu_events_ejdr->Scale(NWS_F);
        h_wz_enu_events_ejdr->Scale(NWS_F);
        h_ttZ_enu_events_ejdr->Scale(NWS_F);
        h_zz_enu_events_ejdr->Scale(NWS_F);
	h_bg_enu_events_ejdr->Scale(NWS_F);

        h_d_enu_events_ejdr->Draw();
/*        h_ww_enu_events_ejdr->Draw("SAME");
        h_wz_enu_events_ejdr->Draw("SAME");
        h_ttZ_enu_events_ejdr->Draw("SAME");
        h_zz_enu_events_ejdr->Draw("SAME");
	h_bg_enu_events_ejdr->Draw("SAME");
//        h_se_enu_events_ejdr->Draw("SAME");
//        h_met_enu_events_ejdr->Draw("SAME");

        h_events_ejdr_canvas->BuildLegend();
        h_events_ejdr_canvas->SaveAs("en_ej_dr.root");
        h_events_ejdr_canvas->SaveAs("en_ej_dr.pdf");



        auto h_d_enu_events_ezdr = d_enu_top_selection.Histo1D({"MC ezdr_enu_Channel","MC e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
        auto h_bg_enu_events_ezdr = bg_enu_top_selection.Histo1D({"BG ezdr_enu_Channel","Background e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
/*	auto h_ww_enu_events_ezdr = ww_enu_top_selection.Histo1D({"WW ezdr_enu_Channel","WW e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
        auto h_wz_enu_events_ezdr = wz_enu_top_selection.Histo1D({"WZ ezdr_enu_Channel","WZ e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
        auto h_ttZ_enu_events_ezdr = ttZ_enu_top_selection.Histo1D({"ttZ ezdr_enu_Channel","ttZ  e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
        auto h_zz_enu_events_ezdr = zz_enu_top_selection.Histo1D({"zz ezdr_enu_Channel","zz e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
//        auto h_se_enu_events_ezdr = se_enu_top_selection.Histo1D({"Single Electron ezdr_enu_Channel","Single Electron e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
//        auto h_met_enu_events_ezdr = met_enu_top_selection.Histo1D({"MET ezdr_enu_Channel","MET e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");

        auto h_events_ezdr_canvas = new TCanvas("ezdeltar", "ezdeltar",10,10,900,900);

        h_d_enu_events_ezdr->GetXaxis()->SetTitle("e and z deltaR");
        h_d_enu_events_ezdr->GetYaxis()->SetTitle("Events");


        h_d_enu_events_ezdr->SetLineColor(kBlack);
	h_bg_enu_events_ezdr->SetLineColor(kPink-6);
/*        h_ww_enu_events_ezdr->SetLineColor(kRed);
        h_wz_enu_events_ezdr->SetLineColor(kOrange);
        h_ttZ_enu_events_ezdr->SetLineColor(kYellow);
        h_zz_enu_events_ezdr->SetLineColor(kTeal);
//        h_se_enu_events_ezdr->SetLineColor(kViolet);
//        h_met_enu_events_ezdr->SetLineColor(kPink);


        h_d_enu_events_ezdr->Scale(NWS_F);
/*        h_ww_enu_events_ezdr->Scale(NWS_F);
        h_wz_enu_events_ezdr->Scale(NWS_F);
        h_ttZ_enu_events_ezdr->Scale(NWS_F);
        h_zz_enu_events_ezdr->Scale(NWS_F);
	h_bg_enu_events_ezdr->Scale(NWS_F);

        h_d_enu_events_ezdr->Draw();
/*        h_ww_enu_events_ezdr->Draw("SAME");
        h_wz_enu_events_ezdr->Draw("SAME");
        h_ttZ_enu_events_ezdr->Draw("SAME");
        h_zz_enu_events_ezdr->Draw("SAME");
	h_bg_enu_events_ezdr->Draw("SAME");
//        h_se_enu_events_ezdr->Draw("SAME");
//        h_met_enu_events_ezdr->Draw("SAME");

        h_events_ezdr_canvas->BuildLegend();
        h_events_ezdr_canvas->SaveAs("en_ez_dr.root");
        h_events_ezdr_canvas->SaveAs("en_ez_dr.pdf");





//////////////////////////////////////////////////////////////////////// Mu-Nu MCs //////////////////////////////////////////////////////////////////////////////////
        auto h_d_munu_events_mupt = d_munu_top_selection.Histo1D({"MC muon_pt_munu_Channel","MC muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
	auto h_bg_munu_events_mupt = bg_munu_top_selection.Histo1D({"bg muon_pt_munu_Channel","Background muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
/*        auto h_ww_munu_events_mupt = ww_munu_top_selection.Histo1D({"WW muon_pt_munu_Channel","WW muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
        auto h_wz_munu_events_mupt = wz_munu_top_selection.Histo1D({"WZ muon_pt_munu_Channel","WZ muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
        auto h_zz_munu_events_mupt = zz_munu_top_selection.Histo1D({"ZZ muon_pt_munu_Channel","ZZ muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
        auto h_ttZ_munu_events_mupt = ttZ_munu_top_selection.Histo1D({"ttZ muon_pt_munu_Channel","ttZ muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
//	auto h_smu_munu_events_mupt = sm_munu_top_selection.Histo1D({"Single Muon muon_pt_munu_Channel","Single muon muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
//	auto h_met_munu_events_mupt = met_munu_top_selection.Histo1D({"MET muon_pt_munu_Channel","MET muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");

        auto h_events_mupt_canvas = new TCanvas("muon pt", "muon pt",10,10,900,900);

	h_d_munu_events_mupt->GetXaxis()->SetTitle("Pt/GeV");
        h_d_munu_events_mupt->GetYaxis()->SetTitle("Events");

        h_d_munu_events_mupt->SetLineColor(kBlack);
/*        h_ww_munu_events_mupt->SetLineColor(kRed);
        h_wz_munu_events_mupt->SetLineColor(kOrange);
        h_ttZ_munu_events_mupt->SetLineColor(kYellow);
        h_zz_munu_events_mupt->SetLineColor(kTeal);
	h_bg_munu_events_mupt->SetLineColor(kPink-6);
//	h_smu_munu_events_mupt->SetLineColor(kViolet);
//	h_met_munu_events_mupt->SetLineColor(kPink);

	h_d_munu_events_mupt->Scale(NWS_F);
/*	h_ww_munu_events_mupt->Scale(NWS_F);
	h_wz_munu_events_mupt->Scale(NWS_F);
	h_ttZ_munu_events_mupt->Scale(NWS_F);
	h_zz_munu_events_mupt->Scale(NWS_F);
	h_bg_munu_events_mupt->Scale(NWS_F);

	h_d_munu_events_mupt->Draw();
/*        h_ww_munu_events_mupt->Draw("SAME");
        h_wz_munu_events_mupt->Draw("SAME");
        h_ttZ_munu_events_mupt->Draw("SAME");
        h_zz_munu_events_mupt->Draw("SAME");
	h_bg_munu_events_mupt->Draw("SAME");
//	h_smu_munu_events_mupt->Draw("SAME");
//	h_met_munu_events_mupt->Draw("SAME");

        h_events_mupt_canvas->BuildLegend();
	h_events_mupt_canvas->SaveAs("munu_pt.root");
	h_events_mupt_canvas->SaveAs("munu_pt.pdf");


	auto h_d_munu_events_jpt = d_munu_top_selection.Histo1D({"MC jet_pt_munu_Channel","MC jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_bg_munu_events_jpt = bg_munu_top_selection.Histo1D({"Bg jet_pt_munu_Channel","Background jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
/*	auto h_ww_munu_events_jpt = ww_munu_top_selection.Histo1D({"WW jet_pt_munu_Channel","WW jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_wz_munu_events_jpt = wz_munu_top_selection.Histo1D({"WZ jet_pt_munu_Channel","WZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_zz_munu_events_jpt = zz_munu_top_selection.Histo1D({"ZZ jet_pt_munu_Channel","ZZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_ttZ_munu_events_jpt = ttZ_munu_top_selection.Histo1D({"ttZ jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
//	auto h_smu_munu_events_jpt = sm_munu_top_selection.Histo1D({"Single muon jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
//	auto h_met_munu_events_jpt = met_munu_top_selection.Histo1D({"MET jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");

	auto h_munu_events_jetpt_canvas = new TCanvas("mu nu jet pt ", "munu jet pt",10,10,900,900);

        h_d_munu_events_jpt->GetXaxis()->SetTitle("pt / GeV");
        h_d_munu_events_jpt->GetYaxis()->SetTitle("Events");

        h_d_munu_events_jpt->SetLineColor(kBlack);
/*        h_ww_munu_events_jpt->SetLineColor(kRed);
        h_wz_munu_events_jpt->SetLineColor(kOrange);
        h_ttZ_munu_events_jpt->SetLineColor(kYellow);
        h_zz_munu_events_jpt->SetLineColor(kTeal);
	h_bg_munu_events_jpt->SetLineColor(kPink-6);
//        h_smu_munu_events_jpt->SetLineColor(kViolet);
//        h_met_munu_events_jpt->SetLineColor(kPink);

        h_d_munu_events_jpt->Scale(NWS_F);
/*        h_ww_munu_events_jpt->Scale(NWS_F);
        h_wz_munu_events_jpt->Scale(NWS_F);
        h_ttZ_munu_events_jpt->Scale(NWS_F);
        h_zz_munu_events_jpt->Scale(NWS_F);
	h_bg_munu_events_jpt->Scale(NWS_F);


        h_d_munu_events_jpt->Draw();
/*        h_ww_munu_events_jpt->Draw("SAME");
        h_wz_munu_events_jpt->Draw("SAME");
        h_ttZ_munu_events_jpt->Draw("SAME");
        h_zz_munu_events_jpt->Draw("SAME");
	h_bg_munu_events_jpt->Draw("SAME");
//        h_smu_munu_events_jpt->Draw("SAME");
//        h_met_munu_events_jpt->Draw("SAME");

        h_munu_events_jetpt_canvas->BuildLegend();
        h_munu_events_jetpt_canvas->SaveAs("munu_jetpt.root");
        h_munu_events_jetpt_canvas->SaveAs("munu_jetpt.pdf");




        auto h_d_munu_events_mueta = d_munu_top_selection.Histo1D({"MC muon_eta_enu_Channel","MC muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
        auto h_bg_munu_events_mueta = bg_munu_top_selection.Histo1D({"Bg muon_eta_enu_Channel","Background muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
/*	auto h_ww_munu_events_mueta = ww_munu_top_selection.Histo1D({"WW muon_eta_enu_Channel","WW muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
        auto h_wz_munu_events_mueta = wz_munu_top_selection.Histo1D({"WZ muon_eta_enu_Channel","WZ muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
        auto h_zz_munu_events_mueta = zz_munu_top_selection.Histo1D({"ZZ muon_eta_enu_Channel","ZZ muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
        auto h_ttZ_munu_events_mueta = ttZ_munu_top_selection.Histo1D({"ttZ muon_eta_enu_Channel","ttZ muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
//        auto h_smu_munu_events_mueta = sm_munu_top_selection.Histo1D({"Single Muon muon_eta_Channel","Single Muon Muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
//        auto h_met_munu_events_mueta = met_munu_top_selection.Histo1D({"MET Muon_eta_Channel","MET muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");


        auto h_events_mueta_canvas = new TCanvas("Muon eta", "Muon eta",10,10,900,900);

        h_d_munu_events_mueta->GetXaxis()->SetTitle("eta");
        h_d_munu_events_mueta->GetYaxis()->SetTitle("Events");

        h_d_munu_events_mueta->SetLineColor(kBlack);
	h_bg_munu_events_mueta->SetLineColor(kPink-6);
/*        h_ww_munu_events_mueta->SetLineColor(kRed);
        h_wz_munu_events_mueta->SetLineColor(kOrange);
        h_ttZ_munu_events_mueta->SetLineColor(kYellow);
        h_zz_munu_events_mueta->SetLineColor(kTeal);
//        h_smu_munu_events_mueta->SetLineColor(kViolet);
//        h_met_munu_events_mueta->SetLineColor(kPink);

        h_d_munu_events_mueta->Scale(NWS_F);
	h_bg_munu_events_mueta->Scale(NWS_F);
/*        h_ww_munu_events_mueta->Scale(NWS_F);
        h_wz_munu_events_mueta->Scale(NWS_F);
        h_ttZ_munu_events_mueta->Scale(NWS_F);
        h_zz_munu_events_mueta->Scale(NWS_F);


        h_d_munu_events_mueta->Draw();
/*        h_ww_munu_events_mueta->Draw("SAME");
        h_wz_munu_events_mueta->Draw("SAME");
        h_ttZ_munu_events_mueta->Draw("SAME");
        h_zz_munu_events_mueta->Draw("SAME");
	h_bg_munu_events_mueta->Draw("SAME");
//        h_smu_munu_events_mueta->Draw("SAME");
//        h_met_munu_events_mueta->Draw("SAME");

        h_events_mueta_canvas->BuildLegend();
        h_events_mueta_canvas->SaveAs("munu_eta.root");
        h_events_mueta_canvas->SaveAs("munu_eta.pdf");



        auto h_d_munu_events_jeta = d_munu_top_selection.Histo1D({"MC jet_eta_munu_Channel","MC jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
	auto h_bg_munu_events_jeta = bg_munu_top_selection.Histo1D({"Bg jet_eta_munu_Channel","Background jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
/*        auto h_ww_munu_events_jeta = ww_munu_top_selection.Histo1D({"WW jet_eta_munu_Channel","WW jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_wz_munu_events_jeta = wz_munu_top_selection.Histo1D({"WZ jet_eta_munu_Channel","WZ jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_zz_munu_events_jeta = zz_munu_top_selection.Histo1D({"ZZ jet_eta_munu_Channel","ZZ jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_ttZ_munu_events_jeta = ttZ_munu_top_selection.Histo1D({"ttZ jet_eta_munu_Channel","ttZ jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
//        auto h_smu_munu_events_jeta = sm_munu_top_selection.Histo1D({"Single Muon jet_eta_munu_Channel","Single Muon jet eta in Muon-neutrino channels",50,-4,4},"tight_jets_eta");
//        auto h_met_munu_events_jeta = met_munu_top_selection.Histo1D({"MET muon_jet_munu_Channel","MET jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");


        auto h_munu_events_jeta_canvas = new TCanvas("jet eta", "jet eta",10,10,900,900);

        h_d_munu_events_jeta->GetXaxis()->SetTitle("eta");
        h_d_munu_events_jeta->GetYaxis()->SetTitle("Events");

        h_d_munu_events_jeta->SetLineColor(kBlack);
/*        h_ww_munu_events_jeta->SetLineColor(kRed);
        h_wz_munu_events_jeta->SetLineColor(kOrange);
        h_ttZ_munu_events_jeta->SetLineColor(kYellow);
        h_zz_munu_events_jeta->SetLineColor(kTeal);
	h_bg_munu_events_jeta->SetLineColor(kPink-6);
//        h_smu_munu_events_jeta->SetLineColor(kViolet);
//        h_met_munu_events_jeta->SetLineColor(kPink);


        h_d_munu_events_jeta->Scale(NWS_F);
/*        h_ww_munu_events_jeta->Scale(NWS_F);
        h_wz_munu_events_jeta->Scale(NWS_F);
        h_ttZ_munu_events_jeta->Scale(NWS_F);
        h_zz_munu_events_jeta->Scale(NWS_F);
	h_bg_munu_events_jeta->Scale(NWS_F);

        h_d_munu_events_jeta->Draw();
/*        h_ww_munu_events_jeta->Draw("SAME");
        h_wz_munu_events_jeta->Draw("SAME");
        h_ttZ_munu_events_jeta->Draw("SAME");
        h_zz_munu_events_jeta->Draw("SAME");
	h_bg_munu_events_jeta->Draw("SAME");
//        h_smu_munu_events_jeta->Draw("SAME");
//        h_met_munu_events_jeta->Draw("SAME");

        h_munu_events_jeta_canvas->BuildLegend();
        h_munu_events_jeta_canvas->SaveAs("munu_jeta.root");
        h_munu_events_jeta_canvas->SaveAs("munu_jeta.pdf");



        auto h_d_munu_events_wmass = d_munu_top_selection.Histo1D({"MC munu_w_mass","MC Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
	auto h_bg_munu_events_wmass = bg_munu_top_selection.Histo1D({"BG munu_w_mass","Background Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
/*        auto h_ww_munu_events_wmass = ww_munu_top_selection.Histo1D({"WW munu_w_mass","WW Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
        auto h_wz_munu_events_wmass = wz_munu_top_selection.Histo1D({"WZ munu_w_mass","WZ Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
        auto h_ttZ_munu_events_wmass = ttZ_munu_top_selection.Histo1D({"ttZ munu_w_mass","ttZ Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
        auto h_zz_munu_events_wmass = zz_munu_top_selection.Histo1D({"ZZ munu_w_mass","ZZ Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
//        auto h_smu_munu_events_wmass = sm_munu_top_selection.Histo1D({"Single Muon munu_w_mass","Single Muon Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
//        auto h_met_munu_events_wmass = met_munu_top_selection.Histo1D({"MET munu_w_mass","MET Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");


        auto h_munu_events_wmass_canvas = new TCanvas("munu w mass ", "munu w mass",10,10,900,900);

        h_d_munu_events_wmass->GetXaxis()->SetTitle("mass/GeV/C^2");
        h_d_munu_events_wmass->GetYaxis()->SetTitle("Events");

        h_d_munu_events_wmass->SetLineColor(kBlack);
/*        h_ww_munu_events_wmass->SetLineColor(kRed);
        h_wz_munu_events_wmass->SetLineColor(kOrange);
        h_ttZ_munu_events_wmass->SetLineColor(kYellow);
        h_zz_munu_events_wmass->SetLineColor(kTeal);
	h_bg_munu_events_wmass->SetLineColor(kPink-6);
//        h_smu_munu_events_wmass->SetLineColor(kViolet);
//        h_met_munu_events_wmass->SetLineColor(kPink);

        h_d_munu_events_wmass->Scale(NWS_F);
/*        h_ww_munu_events_wmass->Scale(NWS_F);
        h_wz_munu_events_wmass->Scale(NWS_F);
        h_ttZ_munu_events_wmass->Scale(NWS_F);
        h_zz_munu_events_wmass->Scale(NWS_F);
	h_bg_munu_events_wmass->Scale(NWS_F);

        h_d_munu_events_wmass->Draw();
/*        h_ww_munu_events_wmass->Draw("SAME");
        h_wz_munu_events_wmass->Draw("SAME");
        h_ttZ_munu_events_wmass->Draw("SAME");
        h_zz_munu_events_wmass->Draw("SAME");
	h_bg_munu_events_wmass->Draw("SAME");
//        h_smu_munu_events_wmass->Draw("SAME");
//        h_met_munu_events_wmass->Draw("SAME");

        h_munu_events_wmass_canvas->BuildLegend();
        h_munu_events_wmass_canvas->SaveAs("munu_transverse_Wmass.root");
        h_munu_events_wmass_canvas->SaveAs("munu_transverse_Wmass.pdf");




        auto h_d_munu_events_zmass = d_munu_top_selection.Histo1D({"MC Z_mass_munu_Channel","MC Z mass in Muon-neutrino channel",50,0,500},"z_mass");
/*        auto h_ww_munu_events_zmass = ww_munu_top_selection.Histo1D({"WW Z_mass_munu_Channel","WW Z mass in Muon-neutrino channel",50,0,500},"z_mass");
        auto h_wz_munu_events_zmass = wz_munu_top_selection.Histo1D({"WZ Z_mass_munu_Channel","WZ Z mass in Muon-neutrino channel",50,0,500},"z_mass");
        auto h_ttZ_munu_events_zmass = ttZ_munu_top_selection.Histo1D({"ttZ Z_mass_munu_Channel","ttZ Z mass in Muon-neutrino channel",50,0,500},"z_mass");
        auto h_zz_munu_events_zmass = zz_munu_top_selection.Histo1D({"zz Z_mass_munu_Channel","zz Z mass in Muon-neutrino channel",50,0,500},"z_mass");
	auto h_bg_munu_events_zmass = bg_munu_top_selection.Histo1D({"bg Z_mass_munu_Channel","Background Z mass in Muon-neutrino channel",50,0,500},"z_mass");
//	auto h_smu_munu_events_zmass = sm_munu_top_selection.Histo1D({"Single Muon Z_mass_enu_Channel","Single Muon Z mass in Muon-neutrino channel",50,0,500},"z_mass");
//	auto h_met_munu_events_zmass = met_munu_top_selection.Histo1D({"MET Z_mass_munu_Channel","MET Z mass in Muon-neutrino channel",50,0,500},"z_mass");

        auto h_munu_events_zmass_canvas = new TCanvas("Z mass", "Z mass",10,10,900,900);

        h_d_munu_events_zmass->GetXaxis()->SetTitle("mass/GeVC^2");
        h_d_munu_events_zmass->GetYaxis()->SetTitle("Events");

        h_d_munu_events_zmass->SetLineColor(kBlack);
/*        h_ww_munu_events_zmass->SetLineColor(kRed);
        h_wz_munu_events_zmass->SetLineColor(kOrange);
        h_ttZ_munu_events_zmass->SetLineColor(kYellow);
        h_bg_munu_events_zmass->SetLineColor(kPink-6);
//	h_zz_munu_events_zmass->SetLineColor(kTeal);
//	h_smu_munu_events_zmass->SetLineColor(kViolet);
//	h_met_munu_events_zmass->SetLineColor(kPink);

	h_d_munu_events_zmass->Scale(NWS_F);
/*        h_ww_munu_events_zmass->Scale(NWS_F);
        h_wz_munu_events_zmass->Scale(NWS_F);
        h_ttZ_munu_events_zmass->Scale(NWS_F);
        h_zz_munu_events_zmass->Scale(NWS_F);
	h_bg_munu_events_zmass->Scale(NWS_F);


	h_d_munu_events_zmass->Draw();
/*        h_ww_munu_events_zmass->Draw("SAME");
        h_wz_munu_events_zmass->Draw("SAME");
        h_ttZ_munu_events_zmass->Draw("SAME");
        h_zz_munu_events_zmass->Draw("SAME");
	h_bg_munu_events_zmass->Draw("SAME");
//	h_smu_munu_events_zmass->Draw("SAME");
//	h_met_munu_events_zmass->Draw("SAME");

	h_munu_events_zmass_canvas->BuildLegend();
        h_munu_events_zmass_canvas->SaveAs("munu_Z_mass.root");
	h_munu_events_zmass_canvas->SaveAs("munu_Z_mass.pdf");



        auto h_d_munu_events_jdphi = d_munu_top_selection.Histo1D({"MC jdphi_munu_Channel","MC jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_bg_munu_events_jdphi = bg_munu_top_selection.Histo1D({"bg jdphi_munu_Channel","Backgrouund jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
/*	auto h_ww_munu_events_jdphi = ww_munu_top_selection.Histo1D({"WW jdphi_munu_Channel","WW jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_wz_munu_events_jdphi = wz_munu_top_selection.Histo1D({"WZ jdphi_munu_Channel","WZ jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_ttZ_munu_events_jdphi = ttZ_munu_top_selection.Histo1D({"ttZ jdphi_munu_Channel","ttZ jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_zz_munu_events_jdphi = zz_munu_top_selection.Histo1D({"zz jdphi_munu_Channel","zz jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
//        auto h_smu_munu_events_jdphi = sm_munu_top_selection.Histo1D({"Single Muon jdphi_munu_Channel","Single Muon jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
//        auto h_met_munu_events_jdphi = met_munu_top_selection.Histo1D({"MET jdphi_munu_Channel","MET jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");


        auto h_munu_events_jdphi_canvas = new TCanvas("jetdeltaphi", "jetdeltaphi",10,10,900,900);

        h_d_munu_events_jdphi->GetXaxis()->SetTitle("jets delta phi / rad");
        h_d_munu_events_jdphi->GetYaxis()->SetTitle("Events");

        h_d_munu_events_jdphi->SetLineColor(kBlack);
/*        h_ww_munu_events_jdphi->SetLineColor(kRed);
        h_wz_munu_events_jdphi->SetLineColor(kOrange);
        h_ttZ_munu_events_jdphi->SetLineColor(kYellow);
        h_zz_munu_events_jdphi->SetLineColor(kTeal);
	h_bg_munu_events_jdphi->SetLineColor(kPink-6);
//        h_smu_munu_events_jdphi->SetLineColor(kViolet);
//        h_met_munu_events_jdphi->SetLineColor(kPink);

        h_d_munu_events_jdphi->Scale(NWS_F);
/*        h_ww_munu_events_jdphi->Scale(NWS_F);
        h_wz_munu_events_jdphi->Scale(NWS_F);
        h_ttZ_munu_events_jdphi->Scale(NWS_F);
        h_zz_munu_events_jdphi->Scale(NWS_F);
	h_bg_munu_events_jdphi->Scale(NWS_F);

        h_d_munu_events_jdphi->Draw();
/*        h_ww_munu_events_jdphi->Draw("SAME");
        h_wz_munu_events_jdphi->Draw("SAME");
        h_ttZ_munu_events_jdphi->Draw("SAME");
        h_zz_munu_events_jdphi->Draw("SAME");
	h_bg_munu_events_jdphi->Draw("SAME");
//        h_smu_munu_events_jdphi->Draw("SAME");
//        h_met_munu_events_jdphi->Draw("SAME");

        h_munu_events_jdphi_canvas->BuildLegend();
        h_munu_events_jdphi_canvas->SaveAs("munu_jet_dphi.root");
        h_munu_events_jdphi_canvas->SaveAs("munu_jet_dphi.pdf");



        auto h_d_munu_events_zmetdphi = d_munu_top_selection.Histo1D({"MC zmetdphi_munu_Channel","MC Z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
/*        auto h_ww_munu_events_zmetdphi = ww_munu_top_selection.Histo1D({"WW zmetdphi_munu_Channel","WW Z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_wz_munu_events_zmetdphi = wz_munu_top_selection.Histo1D({"WZ zmetdphi_munu_Channel","WZ Z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_ttZ_munu_events_zmetdphi = ttZ_munu_top_selection.Histo1D({"ttZ zmetdphi_munu_Channel","ttZ Z met pt jet deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_zz_munu_events_zmetdphi = zz_munu_top_selection.Histo1D({"zz zmetdphi_munu_Channel","zz z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
	auto h_bg_munu_events_zmetdphi = bg_munu_top_selection.Histo1D({"BG zmetdphi_munu_Channel","Background Z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
//        auto h_smu_munu_events_zmetdphi = sm_munu_top_selection.Histo1D({"Single Muon zmetdphi_munu_Channel","Single Muon z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
//        auto h_met_munu_events_zmetdphi = met_munu_top_selection.Histo1D({"MET zmetdphi_munu_Channel","MET z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");



        auto h_munu_events_zmetdphi_canvas = new TCanvas("zmetdeltaphi", "zmetdeltaphi",10,10,900,900);

        h_d_munu_events_zmetdphi->GetXaxis()->SetTitle("z and met delta phi / rad");
        h_d_munu_events_zmetdphi->GetYaxis()->SetTitle("Events");

        h_d_munu_events_zmetdphi->SetLineColor(kBlack);
/*        h_ww_munu_events_zmetdphi->SetLineColor(kRed);
        h_wz_munu_events_zmetdphi->SetLineColor(kOrange);
        h_ttZ_munu_events_zmetdphi->SetLineColor(kYellow);
        h_zz_munu_events_zmetdphi->SetLineColor(kTeal);
	h_bg_munu_events_zmetdphi->SetLineColor(kPink-6);
//        h_smu_munu_events_zmetdphi->SetLineColor(kViolet);
//        h_met_munu_events_zmetdphi->SetLineColor(kPink);

        h_d_munu_events_zmetdphi->Scale(NWS_F);
        h_ww_munu_events_zmetdphi->Scale(NWS_F);
        h_wz_munu_events_zmetdphi->Scale(NWS_F);
        h_ttZ_munu_events_zmetdphi->Scale(NWS_F);
        h_zz_munu_events_zmetdphi->Scale(NWS_F);
	h_bg_munu_events_zmetdphi->Scale(NWS_F);

        h_d_munu_events_zmetdphi->Draw();
/*        h_ww_munu_events_zmetdphi->Draw("SAME");
        h_wz_munu_events_zmetdphi->Draw("SAME");
        h_ttZ_munu_events_zmetdphi->Draw("SAME");
        h_zz_munu_events_zmetdphi->Draw("SAME");
	h_bg_munu_events_zmetdphi->Draw("SAME");
//        h_smu_munu_events_zmetdphi->Draw("SAME");
//        h_met_munu_events_zmetdphi->Draw("SAME");

        h_munu_events_zmetdphi_canvas->BuildLegend();
        h_munu_events_zmetdphi_canvas->SaveAs("munu_zmetpt_dphi.root");
        h_munu_events_zmetdphi_canvas->SaveAs("munu_zmetpt_dphi.pdf");


        auto h_d_munu_events_zwdphi = d_munu_top_selection.Histo1D({"MC zwdphi_munu_Channel","MC Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_bg_munu_events_zwdphi = bg_munu_top_selection.Histo1D({"BG zwdphi_munu_Channel","Background Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
/*	auto h_ww_munu_events_zwdphi = ww_munu_top_selection.Histo1D({"WW zwdphi_munu_Channel","WW Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_wz_munu_events_zwdphi = wz_munu_top_selection.Histo1D({"WZ zwdphi_munu_Channel","WZ Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_ttZ_munu_events_zwdphi = ttZ_munu_top_selection.Histo1D({"ttZ zwdphi_munu_Channel","ttZ  Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_zz_munu_events_zwdphi = zz_munu_top_selection.Histo1D({"zz zwdphi_munu_Channel","zz z Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
//        auto h_smu_munu_events_zwdphi = sm_munu_top_selection.Histo1D({"Single Muon zmetdphi_munu_Channel","Single Muon z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
//        auto h_met_munu_events_zwdphi = met_munu_top_selection.Histo1D({"MET zmetdphi_munu_Channel","MET z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");


        auto h_munu_events_zwdphi_canvas = new TCanvas("zwdeltaphi", "zwdeltaphi",10,10,900,900);

        h_d_munu_events_zwdphi->GetXaxis()->SetTitle("z and w delta phi / rad");
        h_d_munu_events_zwdphi->GetYaxis()->SetTitle("Events");

        h_d_munu_events_zwdphi->SetLineColor(kBlack);
/*        h_ww_munu_events_zwdphi->SetLineColor(kRed);
        h_wz_munu_events_zwdphi->SetLineColor(kOrange);
        h_ttZ_munu_events_zwdphi->SetLineColor(kYellow);
        h_zz_munu_events_zwdphi->SetLineColor(kTeal);
	h_bg_munu_events_zwdphi->SetLineColor(kPink-6);
//        h_smu_munu_events_zwdphi->SetLineColor(kViolet);
//        h_met_munu_events_zwdphi->SetLineColor(kPink);

        h_d_munu_events_zwdphi->Scale(NWS_F);
/*        h_ww_munu_events_zwdphi->Scale(NWS_F);
        h_wz_munu_events_zwdphi->Scale(NWS_F);
        h_ttZ_munu_events_zwdphi->Scale(NWS_F);
        h_zz_munu_events_zwdphi->Scale(NWS_F);
	h_bg_munu_events_zwdphi->Scale(NWS_F);

        h_d_munu_events_zwdphi->Draw();
/*        h_ww_munu_events_zwdphi->Draw("SAME");
        h_wz_munu_events_zwdphi->Draw("SAME");
        h_ttZ_munu_events_zwdphi->Draw("SAME");
        h_zz_munu_events_zmetdphi->Draw("SAME");
	h_bg_munu_events_zmetdphi->Draw("SAME");
//        h_smu_munu_events_zwdphi->Draw("SAME");
//        h_met_munu_events_zwdphi->Draw("SAME");

        h_munu_events_zwdphi_canvas->BuildLegend();
        h_munu_events_zwdphi_canvas->SaveAs("munu_zwd_dphi.root");
        h_munu_events_zwdphi_canvas->SaveAs("munu_zw_dphi.pdf");




        auto h_d_munu_events_mujdr = d_munu_top_selection.Histo1D({"MC ejdr_munu_Channel","MC muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
        auto h_bg_munu_events_mujdr = bg_munu_top_selection.Histo1D({"BG ejdr_munu_Channel","Background muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
/*	auto h_ww_munu_events_mujdr = ww_munu_top_selection.Histo1D({"WW ejdr_munu_Channel","WW muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
        auto h_wz_munu_events_mujdr = wz_munu_top_selection.Histo1D({"WZ ejdr_munu_Channel","WZ muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
        auto h_ttZ_munu_events_mujdr = ttZ_munu_top_selection.Histo1D({"ttZ ejdr_munu_Channel","ttZ muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
        auto h_zz_munu_events_mujdr = zz_munu_top_selection.Histo1D({"zz ejdr_munu_Channel","zz muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
//        auto h_smu_munu_events_mujdr = sm_munu_top_selection.Histo1D({"Single Muon mujdr_munu_Channel","Single Muon muon and jet pt deltaR in Muon-neutrino channel",50,0,10},"jet_mu_min_dR");
//        auto h_met_munu_events_mujdr = met_munu_top_selection.Histo1D({"MET mujdr_munu_Channel","MET muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");

        auto h_events_mujdr_canvas = new TCanvas("mujdeltar", "mujdeltar",10,10,900,900);

        h_d_munu_events_mujdr->GetXaxis()->SetTitle("mu and jet deltaR");
        h_d_munu_events_mujdr->GetYaxis()->SetTitle("Events");



        h_d_munu_events_mujdr->SetLineColor(kBlack);
/*        h_ww_munu_events_mujdr->SetLineColor(kRed);
        h_wz_munu_events_mujdr->SetLineColor(kOrange);
        h_ttZ_munu_events_mujdr->SetLineColor(kYellow);
        h_zz_munu_events_mujdr->SetLineColor(kTeal);
	h_bg_munu_events_mujdr->SetLineColor(kPink-6);
//        h_smu_munu_events_mujdr->SetLineColor(kViolet);
//        h_met_munu_events_mujdr->SetLineColor(kPink);

        h_d_munu_events_mujdr->Scale(NWS_F);
/*        h_ww_munu_events_mujdr->Scale(NWS_F);
        h_wz_munu_events_mujdr->Scale(NWS_F);
        h_ttZ_munu_events_mujdr->Scale(NWS_F);
        h_zz_munu_events_mujdr->Scale(NWS_F);
	h_bg_munu_events_mujdr->Scale(NWS_F);


        h_d_munu_events_mujdr->Draw();
/*        h_ww_munu_events_mujdr->Draw("SAME");
        h_wz_munu_events_mujdr->Draw("SAME");
        h_ttZ_munu_events_mujdr->Draw("SAME");
        h_zz_munu_events_mujdr->Draw("SAME");
	h_bg_munu_events_mujdr->Draw("SAME");
//        h_smu_munu_events_mujdr->Draw("SAME");
//        h_met_munu_events_mujdr->Draw("SAME");

        h_events_mujdr_canvas->BuildLegend();
        h_events_mujdr_canvas->SaveAs("munu_muj_dr.root");
        h_events_mujdr_canvas->SaveAs("munu_muj_dr.pdf");



        auto h_d_munu_events_muzdr = d_munu_top_selection.Histo1D({"MC muzdr_munu_Channel","MC mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
	auto h_bg_munu_events_muzdr = bg_munu_top_selection.Histo1D({"BG muzdr_munu_Channel","Background mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
/*        auto h_ww_munu_events_muzdr = ww_munu_top_selection.Histo1D({"WW muzdr_munu_Channel","WW mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
        auto h_wz_munu_events_muzdr = wz_munu_top_selection.Histo1D({"WZ muzdr_munu_Channel","WZ mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
        auto h_ttZ_munu_events_muzdr = ttZ_munu_top_selection.Histo1D({"ttZ muzdr_munu_Channel","ttZ  mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
        auto h_zz_munu_events_muzdr = zz_munu_top_selection.Histo1D({"zz muzdr_munu_Channel","zz mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
//        auto h_smu_munu_events_muzdr = sm_munu_top_selection.Histo1D({"Single Muon muzdr_munu_Channel","Single Muon mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
//        auto h_met_munu_events_muzdr = met_munu_top_selection.Histo1D({"MET muzdr_munu_Channel","MET mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");

        auto h_events_muzdr_canvas = new TCanvas("muzdeltar", "muzdeltar",10,10,900,900);

        h_d_munu_events_muzdr->GetXaxis()->SetTitle("muon and z deltaR");
        h_d_munu_events_muzdr->GetYaxis()->SetTitle("Events");


        h_d_munu_events_muzdr->SetLineColor(kBlack);
/*        h_ww_munu_events_muzdr->SetLineColor(kRed);
        h_wz_munu_events_muzdr->SetLineColor(kOrange);
        h_ttZ_munu_events_muzdr->SetLineColor(kYellow);
        h_zz_munu_events_muzdr->SetLineColor(kTeal);
	h_bg_munu_events_muzdr->SetLineColor(kPink-6);
//        h_smu_munu_events_muzdr->SetLineColor(kViolet);
//        h_met_munu_events_muzdr->SetLineColor(kPink);


        h_d_munu_events_muzdr->Scale(NWS_F);
/*        h_ww_munu_events_muzdr->Scale(NWS_F);
        h_wz_munu_events_muzdr->Scale(NWS_F);
        h_ttZ_munu_events_muzdr->Scale(NWS_F);
        h_zz_munu_events_muzdr->Scale(NWS_F);
	h_bg_munu_events_muzdr->Scale(NWS_F);


        h_d_munu_events_muzdr->Draw();
/*        h_ww_munu_events_muzdr->Draw("SAME");
        h_wz_munu_events_muzdr->Draw("SAME");
        h_ttZ_munu_events_muzdr->Draw("SAME");
        h_zz_munu_events_muzdr->Draw("SAME");
	h_bg_munu_events_muzdr->Draw("SAME");
//        h_smu_munu_events_muzdr->Draw("SAME");
//        h_met_munu_events_muzdr->Draw("SAME");

        h_events_muzdr_canvas->BuildLegend();
        h_events_muzdr_canvas->SaveAs("munu_ez_dr.root");
        h_events_muzdr_canvas->SaveAs("munu_ez_dr.pdf");




////////////////////////////////////////////////////////////////////////////////// Pt And ETa Test/////////////////////////////////////////////////////////////////
/*
	auto electron_pT = d.Histo1D({"Electron pT","Electron_pT",20,0,200},"Electron_pt");
	auto electron_pTCanvas = new TCanvas("Electron pT", "Electron_pT",10,10,700,700);
	electron_pT->GetXaxis()->SetTitle("pT/GeV");
        electron_pT->GetYaxis()->SetTitle("Events");
	electron_pT->SetLineColor(kRed);
        electron_pTCanvas->BuildLegend();
        electron_pT->Draw();
	electron_pTCanvas->SaveAs("Pure_electron_pT.root");


	auto electron_eta = d.Histo1D({"Electron eta","Electron_eta",20,-5,5},"Electron_eta");
        auto electron_etaCanvas = new TCanvas("Electron eta", "Electron_eta",10,10,700,700);
        electron_eta->GetXaxis()->SetTitle("eta");
        electron_eta->GetYaxis()->SetTitle("Events");
        electron_eta->SetLineColor(kRed);
        electron_etaCanvas->BuildLegend();
        electron_eta->Draw();
        electron_etaCanvas->SaveAs("Pure_electron_eta.root");
*/
/////////////////////////////////////////////////////////////////  Z Recon from U Histo2D /////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TO DO
// you need correlation of delta phi vs z mass reconstruction.
//later on W reconstruction : and lepton and MET, delta phi vs w recons, real data : reconstruct trans. mass

///////////////////////////////////////////////////////// THSTACK IN PROGRESS //////////////////////////////////////////////////////////////////////////////////////



	//auto NumDist = new THStack("NumDist","Quark Dist. Per event");//defining a stack that is made up many histograms
   	//auto temp_bnumHist = (TH1D*)&bnumHist.GetValue();//need to copy bnumHist on another hist to make it work with thstack
	//temp_bnumHist->SetLineColor(kRed); // saying to stack use colour red for the histogram
	//bnumHist->SetLineColor(kRed);
	//NumDist->Add(TH1 &bnumHist);
   	//NumDist->Add(temp_bnumHist); //adding the histogram to my stack
   	//auto temp_cnumHist = (TH1D*)&cnumHist.GetValue();
	//temp_cnumHist->SetLineColor(kBlue); //adding another one and so on....
   	//NumDist->Add(temp_cnumHist);
   	//h3->SetLineColor(kGreen);
   	//hs->Add((TH1*)&h3);
   	//auto QDist = new TCanvas("Quark Dist.","Quark Dist",20,20,400,500); //making a canvas for the stack
   	//auto T = new TText(.5,.5,"TStack");
	//T->SetTextFont(42); T->SetTextAlign(21);
   	//QDist->Divide(2,2);
   	//QDist->cd(2); NumDist->Draw("Quark Dist."); T->DrawTextNDC(.5,.95,"Quark Dist");//this is the stack canvas style used in phd
	//QDist->SaveAs("QDist.root");// save it


	//auto hist1_ptr = (TH1D*)&hist1.GetValue();
	//hs->Add((TH1*)&hist1_ptr)



  	// std::cout << "Number of Z Bosons "<<*Znum<<std::endl;
	// std::cout << "Number of W Bosons "<<*Wnum<<std::endl;
	// std::cout << "Number of t "<<*tnum<<std::endl;
	// std::cout << "Number of tbar "<<*tbarnum<<std::endl;
	// std::cout << "Number of c "<<*cnum<<std::endl;
	// std::cout << "Number of cbar "<<*cbarnum<<std::endl;
	// std::cout << "Number of u "<<*unum<<std::endl;
	// std::cout << "Number of ubar "<<*ubarnum<<std::endl;
	// std::cout << "Number of b "<<*bnum<<std::endl;
	// std::cout << "Number of bbar "<<*bbarnum<<std::endl;
	// std::cout << "Number of d "<<*dnum<<std::endl;
	// std::cout << "Number of dbar "<<*dbarnum<<std::endl;
	// std::cout << "Number of s "<<*snum<<std::endl;
	// std::cout << "Number of sbar "<<*sbarnum<<std::endl;
	// std::cout << "Number of ele "<<*elenum<<std::endl;
	// std::cout << "Number of pos "<<*posnum<<std::endl;
	// std::cout << "Number of muon "<<*muonnum<<std::endl;
	// std::cout << "Number of muonbar "<<*muonbarnum<<std::endl;
	// std::cout << "Number of elenu "<<*elenunum<<std::endl;
	// std::cout << "Number of posnu "<<*posnunum<<std::endl;
	// std::cout << "Number of muonnu "<<*muonnunum<<std::endl;
	// std::cout << "Number of muonbarnu "<<*muonbarnunum<<std::endl;
	//	std::cout << "Number of mean b "<<*bnumMean<<std::endl;

    // Print cut report
    //auto allCutsReport{d.Report()};
    //for (auto&& cutInfo: allCutsReport)
    //{
    //    std::cout << cutInfo.GetName() << '\t' << cutInfo.GetAll() << '\t' << cutInfo.GetPass() << '\t' << cutInfo.GetEff() << " %" << std::endl;
    //}
}

