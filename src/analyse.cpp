#include "analyse.hpp"
#include <algorithm>
#include "TLorentzVector.h"
#include "sf.hpp"
#include <TCanvas.h>
#include <TText.h>
#include <THStack.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <TStyle.h>

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

namespace
{

constexpr double MAX_ELE_NUM{1};
constexpr double MIN_ELE_PT{45}; //{15}//min 12, AP 45, 
constexpr float MIN_ELE_LEADING_PT{35.f};
constexpr double MAX_ELE_ETA{2.5};
constexpr double ENDCAP_MIN_ETA{1.566};
constexpr double BARREL_MAX_ETA{1.4442};

constexpr double MAX_MU_NUM{1};
constexpr double MIN_MU_PT{40};//min 33 ,40 AP
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
constexpr double NWS_F = TZQ_W*WWLNQQ_W*WZLNQQ_W*TTZQQ_W*ZZLLQQ_W; //normalization scale factor is the product of all scale factors.

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

   	ROOT::RDataFrame d{"Events", "/data/disk3/nanoAOD_2017/tZqlvqq/*.root"};
	ROOT::RDataFrame ww{"Events", "/data/disk0/nanoAOD_2017/WWToLNuQQ/*.root"};
	ROOT::RDataFrame wz{"Events", "/data/disk0/nanoAOD_2017/WZTo1L1Nu2Q/*.root"};
	ROOT::RDataFrame ttZ{"Events", "/data/disk0/nanoAOD_2017/ttZToQQ/*.root"};
	ROOT::RDataFrame zz{"Events", "/data/disk0/nanoAOD_2017/ZZTo2L2Q/*.root"};
	//ROOT::RDataFrame Se{"Events","/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_Run*/*.root"};
	//ROOT::RDataFrame Sm{"Events","/data/disk3/nanoAOD_2017/SingleMuon_*/*.root"};
	//ROOT::RDataFrame met{"Events","/data/disk0/nanoAOD_2017/MET*/*.root"};
	/*auto d = dc.Range(0, 100);
	/auto ww = wwc.Range(0, 100);
	auto wz = wzc.Range(0, 100);
	auto ttZ = ttZc.Range(0, 100);
	auto zz = zzc.Range(0, 100);
	auto se = Sec.Range(0, 100);
	auto sm = Smc.Range(0, 100);
	auto met = metc.Range(0, 100);
	*/
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
         bptHist->SetFillColor(kBlue);
 	bptHist->DrawClone();
         bptcanvas->SaveAs("bptHist.root");


// /////////////////////////////////////////////////////////////////// Eta Histograms //////////////////////////////////////////////////////////////////////////////////
 	std::cout << "eta histogram"<<std::endl;
         auto betafunc = [](const ints id, const floats eta) {return eta[id ==5];};
         auto betaHist = d.Define("betas",betafunc,{"GenPart_pdgId","GenPart_eta"})
                        .Filter([](floats betas){return std::all_of(betas.cbegin(), betas.cend(), [](float eta){return abs(eta) <= 2.5;});},{"betas"})
                         .Histo1D({"beta","Eta of b quarks per event",20,-5,5},"betas"); // This is how a histogtam is done, sayin$
         auto betacanvas = new TCanvas("Eta of b quarks per event", "Eta of b quarks per event",10,10,700,700);
         betaHist->SetFillColor(kBlue);
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
	dphihist->SetFillColor(kRed);
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


	auto bjet_id{[](const ints& tight_jets,	const floats& btags, const floats& etas){
		return  tight_jets && (btags > 0.8838f) && (etas < 2.4f);
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
		} //maybe use select<floats> from the column?
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


	auto ScaleFact_func{[](const double& i){// this function make the variable which is equal to one and will be used for all scale factors
		return 1;
	}};

	auto NormScaleFact_func{[](const int& i){// this function calculates the weight scale factor
		double weight{TZQ_W*WWLNQQ_W*WZLNQQ_W*TTZQQ_W*ZZLLQQ_W};
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
							.Define("tight_jets_deltaphi", jet_deltaphi_func, {"tight_jets_phi"})
                   					.Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

	auto d_enu_jets_bjets_selection = d_enu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                    					        .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto d_enu_z_rec_selection = d_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});

///////////////////////////////////////////////////////////////////////// Muon Channel/////////////////////////////////////////////////////////////////////////////
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

        auto d_munu_jets_bjets_selection = d_munu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto d_munu_z_rec_selection = d_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});



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

	auto ww_enu_jets_bjets_selection = ww_enu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto ww_enu_z_rec_selection = ww_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});

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

        auto ww_munu_jets_bjets_selection = ww_munu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto ww_munu_z_rec_selection = ww_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});

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

	auto wz_enu_jets_bjets_selection = wz_enu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto wz_enu_z_rec_selection = wz_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});

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

        auto wz_munu_jets_bjets_selection = wz_munu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto wz_munu_z_rec_selection = wz_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});


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

	auto ttZ_enu_jets_bjets_selection = ttZ_enu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto ttZ_enu_z_rec_selection = ttZ_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});

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

        auto ttZ_munu_jets_bjets_selection = ttZ_munu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto ttZ_munu_z_rec_selection = ttZ_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});



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

	auto zz_enu_jets_bjets_selection = zz_enu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto zz_enu_z_rec_selection = zz_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
							.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});

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

        auto zz_munu_jets_bjets_selection = zz_munu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto zz_munu_z_rec_selection = zz_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});


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

	auto se_enu_jets_bjets_selection = se_enu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto se_enu_z_rec_selection = se_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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

        auto sm_munu_jets_bjets_selection = sm_munu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto sm_munu_z_rec_selection = sm_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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

	auto met_enu_jets_bjets_selection = met_enu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                    					      .Filter(bjet_cut, {"bjets"}, "b jet cut");


	auto met_enu_z_rec_selection = met_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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

        auto met_munu_jets_bjets_selection = met_munu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                                                              .Filter(bjet_cut, {"bjets"}, "b jet cut");

  	auto met_munu_z_rec_selection = met_munu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                                                                .Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
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
/////////////////////////////////////////////////////////////////////// E-Nu Channel histograms AND Canvases /////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////PT//////////////////////////////////////////////////////////////////////////////////
        auto h_d_enu_events_ept = d_enu_top_selection.Histo1D({"MC electron_pt_enu_Channel","MC electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
        auto h_ww_enu_events_ept = ww_enu_top_selection.Histo1D({"WW electron_pt_enu_Channel","WW electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
        auto h_wz_enu_events_ept = wz_enu_top_selection.Histo1D({"WZ electron_pt_enu_Channel","WZ electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
        auto h_zz_enu_events_ept = zz_enu_top_selection.Histo1D({"ZZ electron_pt_enu_Channel","ZZ electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
        auto h_ttZ_enu_events_ept = ttZ_enu_top_selection.Histo1D({"ttZ electron_pt_enu_Channel","ttZ electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
//	auto h_se_enu_events_ept = se_enu_top_selection.Histo1D({"Single Electron electron_pt_enu_Channel","Single Electron electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");
//	auto h_met_enu_events_ept = met_enu_top_selection.Histo1D({"MET electron_pt_enu_Channel","MET electron pt in electron-neutrino channel",50,0,300},"tight_ele_pt");

        auto h_events_ept_canvas = new TCanvas("electron pt", "electron pt",10,10,900,900);

	h_d_enu_events_ept->GetXaxis()->SetTitle("Pt/GeV");
        h_d_enu_events_ept->GetYaxis()->SetTitle("Events");

        h_d_enu_events_ept->SetFillColor(kBlack);
        h_ww_enu_events_ept->SetFillColor(kRed);
        h_wz_enu_events_ept->SetFillColor(kOrange);
        h_ttZ_enu_events_ept->SetFillColor(kYellow);
        h_zz_enu_events_ept->SetFillColor(kTeal);
//	h_se_enu_events_ept->SetFillColor(kViolet);
//	h_met_enu_events_ept->SetFillColor(kPink);

	h_d_enu_events_ept->Scale(NWS_F);
	h_ww_enu_events_ept->Scale(NWS_F);
	h_wz_enu_events_ept->Scale(NWS_F);
	h_ttZ_enu_events_ept->Scale(NWS_F);
	h_zz_enu_events_ept->Scale(NWS_F);


	h_d_enu_events_ept->Draw();
        h_ww_enu_events_ept->Draw("SAME");
        h_wz_enu_events_ept->Draw("SAME");
        h_ttZ_enu_events_ept->Draw("SAME");
        h_zz_enu_events_ept->Draw("SAME");
//	h_se_enu_events_ept->Draw("SAME");
//	h_met_enu_events_ept->Draw("SAME");

        h_events_ept_canvas->BuildLegend();
	h_events_ept_canvas->SaveAs("enu_pt.root");
	h_events_ept_canvas->SaveAs("enu_pt.pdf");


	auto h_d_enu_events_jpt = d_enu_top_selection.Histo1D({"MC jet_pt_enu_Channel","MC jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_ww_enu_events_jpt = ww_enu_top_selection.Histo1D({"WW jet_pt_enu_Channel","WW jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_wz_enu_events_jpt = wz_enu_top_selection.Histo1D({"WZ jet_pt_enu_Channel","WZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_zz_enu_events_jpt = zz_enu_top_selection.Histo1D({"ZZ jet_pt_enu_Channel","ZZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_ttZ_enu_events_jpt = ttZ_enu_top_selection.Histo1D({"ttZ jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
//	auto h_se_enu_events_jpt = se_enu_top_selection.Histo1D({"Single Electron jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");
//	auto h_met_enu_events_jpt = met_enu_top_selection.Histo1D({"MET jet_pt_enu_Channel","ttZ jet pt in electron-neutrino channel",50,0,300},"tight_jets_pt");

	auto h_events_jetpt_canvas = new TCanvas("e nu jet pt ", "enu jet pt",10,10,900,900);

        h_d_enu_events_jpt->GetXaxis()->SetTitle("pt / GeV");
        h_d_enu_events_jpt->GetYaxis()->SetTitle("Events");
        h_d_enu_events_jpt->SetFillColor(kBlack);
        h_ww_enu_events_jpt->SetFillColor(kRed);
        h_wz_enu_events_jpt->SetFillColor(kOrange);
        h_ttZ_enu_events_jpt->SetFillColor(kYellow);
        h_zz_enu_events_jpt->SetFillColor(kTeal);
//        h_se_enu_events_jpt->SetFillColor(kViolet);
//        h_met_enu_events_jpt->SetFillColor(kPink);

        h_d_enu_events_jpt->Scale(NWS_F);
        h_ww_enu_events_jpt->Scale(NWS_F);
        h_wz_enu_events_jpt->Scale(NWS_F);
        h_ttZ_enu_events_jpt->Scale(NWS_F);
        h_zz_enu_events_jpt->Scale(NWS_F);


        h_d_enu_events_jpt->Draw();
        h_ww_enu_events_jpt->Draw("SAME");
        h_wz_enu_events_jpt->Draw("SAME");
        h_ttZ_enu_events_jpt->Draw("SAME");
        h_zz_enu_events_jpt->Draw("SAME");
//        h_se_enu_events_jpt->Draw("SAME");
//        h_met_enu_events_jpt->Draw("SAME");

        h_events_jetpt_canvas->BuildLegend();
        h_events_jetpt_canvas->SaveAs("enu_jetpt.root");
        h_events_jetpt_canvas->SaveAs("enu_jetpt.pdf");




        auto h_d_enu_events_eeta = d_enu_top_selection.Histo1D({"MC electron_eta_enu_Channel","MC electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
        auto h_ww_enu_events_eeta = ww_enu_top_selection.Histo1D({"WW electron_eta_enu_Channel","WW electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
        auto h_wz_enu_events_eeta = wz_enu_top_selection.Histo1D({"WZ electron_eta_enu_Channel","WZ electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
        auto h_zz_enu_events_eeta = zz_enu_top_selection.Histo1D({"ZZ electron_eta_enu_Channel","ZZ electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
        auto h_ttZ_enu_events_eeta = ttZ_enu_top_selection.Histo1D({"ttZ electron_eta_enu_Channel","ttZ electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
//        auto h_se_enu_events_eeta = se_enu_top_selection.Histo1D({"Single Electron electron_eta_enu_Channel","Single Electron electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");
//        auto h_met_enu_events_eeta = met_enu_top_selection.Histo1D({"MET electron_eta_enu_Channel","MET electron eta in electron-neutrino channel",50,-4,4},"tight_ele_eta");


        auto h_events_eeta_canvas = new TCanvas("electron eta", "electron eta",10,10,900,900);

        h_d_enu_events_eeta->GetXaxis()->SetTitle("eta");
        h_d_enu_events_eeta->GetYaxis()->SetTitle("Events");

        h_d_enu_events_eeta->SetFillColor(kBlack);
        h_ww_enu_events_eeta->SetFillColor(kRed);
        h_wz_enu_events_eeta->SetFillColor(kOrange);
        h_ttZ_enu_events_eeta->SetFillColor(kYellow);
        h_zz_enu_events_eeta->SetFillColor(kTeal);
//        h_se_enu_events_eeta->SetFillColor(kViolet);
//        h_met_enu_events_eeta->SetFillColor(kPink);

        h_d_enu_events_eeta->Scale(NWS_F);
        h_ww_enu_events_eeta->Scale(NWS_F);
        h_wz_enu_events_eeta->Scale(NWS_F);
        h_ttZ_enu_events_eeta->Scale(NWS_F);
        h_zz_enu_events_eeta->Scale(NWS_F);


        h_d_enu_events_eeta->Draw();
        h_ww_enu_events_eeta->Draw("SAME");
        h_wz_enu_events_eeta->Draw("SAME");
        h_ttZ_enu_events_eeta->Draw("SAME");
        h_zz_enu_events_eeta->Draw("SAME");
//        h_se_enu_events_eeta->Draw("SAME");
//        h_met_enu_events_eeta->Draw("SAME");

        h_events_eeta_canvas->BuildLegend();
        h_events_eeta_canvas->SaveAs("enu_eta.root");
        h_events_eeta_canvas->SaveAs("enu_eta.pdf");



        auto h_d_enu_events_jeta = d_enu_top_selection.Histo1D({"MC jet_eta_enu_Channel","MC jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_ww_enu_events_jeta = ww_enu_top_selection.Histo1D({"WW jet_eta_enu_Channel","WW jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_wz_enu_events_jeta = wz_enu_top_selection.Histo1D({"WZ jet_eta_enu_Channel","WZ jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_zz_enu_events_jeta = zz_enu_top_selection.Histo1D({"ZZ jet_eta_enu_Channel","ZZ jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_ttZ_enu_events_jeta = ttZ_enu_top_selection.Histo1D({"ttZ jet_eta_enu_Channel","ttZ jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");
//        auto h_se_enu_events_jeta = se_enu_top_selection.Histo1D({"Single Electron jet_eta_enu_Channel","Single Electron jet eta in electron-neutrino channels",50,-4,4},"tight_jets_eta");
//        auto h_met_enu_events_jeta = met_enu_top_selection.Histo1D({"MET electron_jet_enu_Channel","MET jet eta in electron-neutrino channel",50,-4,4},"tight_jets_eta");


        auto h_events_jeta_canvas = new TCanvas("jet eta", "jet eta",10,10,900,900);

        h_d_enu_events_jeta->GetXaxis()->SetTitle("eta");
        h_d_enu_events_jeta->GetYaxis()->SetTitle("Events");

        h_d_enu_events_jeta->SetFillColor(kBlack);
        h_ww_enu_events_jeta->SetFillColor(kRed);
        h_wz_enu_events_jeta->SetFillColor(kOrange);
        h_ttZ_enu_events_jeta->SetFillColor(kYellow);
        h_zz_enu_events_jeta->SetFillColor(kTeal);
//        h_se_enu_events_jeta->SetFillColor(kViolet);
//        h_met_enu_events_jeta->SetFillColor(kPink);


        h_d_enu_events_jeta->Scale(NWS_F);
        h_ww_enu_events_jeta->Scale(NWS_F);
        h_wz_enu_events_jeta->Scale(NWS_F);
        h_ttZ_enu_events_jeta->Scale(NWS_F);
        h_zz_enu_events_jeta->Scale(NWS_F);


        h_d_enu_events_jeta->Draw();
        h_ww_enu_events_jeta->Draw("SAME");
        h_wz_enu_events_jeta->Draw("SAME");
        h_ttZ_enu_events_jeta->Draw("SAME");
        h_zz_enu_events_jeta->Draw("SAME");
//        h_se_enu_events_jeta->Draw("SAME");
//        h_met_enu_events_jeta->Draw("SAME");

        h_events_jeta_canvas->BuildLegend();
        h_events_jeta_canvas->SaveAs("enu_jeta.root");
        h_events_jeta_canvas->SaveAs("enu_jeta.pdf");



        auto h_d_enu_events_wmass = d_enu_top_selection.Histo1D({"MC enu_w_mass","MC electron-neutrino transverse w mass",50,0,500},"w_e_mass");
        auto h_ww_enu_events_wmass = ww_enu_top_selection.Histo1D({"WW enu_w_mass","WW electron-neutrino transverse w mass",50,0,500},"w_e_mass");
        auto h_wz_enu_events_wmass = wz_enu_top_selection.Histo1D({"WZ enu_w_mass","WZ electron-neutrino transverse w mass",50,0,500},"w_e_mass");
        auto h_ttZ_enu_events_wmass = ttZ_enu_top_selection.Histo1D({"ttZ enu_w_mass","ttZ electron-neutrino transverse w mass",50,0,500},"w_e_mass");
        auto h_zz_enu_events_wmass = zz_enu_top_selection.Histo1D({"ZZ enu_w_mass","ZZ electron-neutrino transverse w mass",50,0,500},"w_e_mass");
//        auto h_se_enu_events_wmass = se_enu_top_selection.Histo1D({"Single Electron enu_w_mass","Single Electron electron-neutrino transverse w mass",50,0,500},"w_e_mass");
//        auto h_met_enu_events_wmass = met_enu_top_selection.Histo1D({"MET enu_w_mass","MET electron-neutrino transverse w mass",50,0,500},"w_e_mass");


        auto h_events_wmass_canvas = new TCanvas("enu w mass ", "enu w mass",10,10,900,900);

        h_d_enu_events_wmass->GetXaxis()->SetTitle("mass/GeV/C^2");
        h_d_enu_events_wmass->GetYaxis()->SetTitle("Events");
        h_d_enu_events_wmass->SetFillColor(kBlack);
        h_ww_enu_events_wmass->SetFillColor(kRed);
        h_wz_enu_events_wmass->SetFillColor(kOrange);
        h_ttZ_enu_events_wmass->SetFillColor(kYellow);
        h_zz_enu_events_wmass->SetFillColor(kTeal);
//        h_se_enu_events_wmass->SetFillColor(kViolet);
//        h_met_enu_events_wmass->SetFillColor(kPink);

        h_d_enu_events_wmass->Scale(NWS_F);
        h_ww_enu_events_wmass->Scale(NWS_F);
        h_wz_enu_events_wmass->Scale(NWS_F);
        h_ttZ_enu_events_wmass->Scale(NWS_F);
        h_zz_enu_events_wmass->Scale(NWS_F);



        h_d_enu_events_wmass->Draw();
        h_ww_enu_events_wmass->Draw("SAME");
        h_wz_enu_events_wmass->Draw("SAME");
        h_ttZ_enu_events_wmass->Draw("SAME");
        h_zz_enu_events_wmass->Draw("SAME");
//        h_se_enu_events_wmass->Draw("SAME");
//        h_met_enu_events_wmass->Draw("SAME");

        h_events_wmass_canvas->BuildLegend();
        h_events_wmass_canvas->SaveAs("enu_transverse_Wmass.root");
        h_events_wmass_canvas->SaveAs("enu_transverse_Wmass.pdf");




        auto h_d_enu_events_zmass = d_enu_z_rec_selection.Histo1D({"MC Z_mass_enu_Channel","MC Z mass in electron-neutrino channel",50,0,500},"z_mass");
        auto h_ww_enu_events_zmass = ww_enu_z_rec_selection.Histo1D({"WW Z_mass_enu_Channel","WW Z mass in electron-neutrino channel",50,0,500},"z_mass");
        auto h_wz_enu_events_zmass = wz_enu_z_rec_selection.Histo1D({"WZ Z_mass_enu_Channel","WZ Z mass in electron-neutrino channel",50,0,500},"z_mass");
        auto h_ttZ_enu_events_zmass = ttZ_enu_z_rec_selection.Histo1D({"ttZ Z_mass_enu_Channel","ttZ Z mass in electron-neutrino channel",50,0,500},"z_mass");
        auto h_zz_enu_events_zmass = zz_enu_z_rec_selection.Histo1D({"zz Z_mass_enu_Channel","zz Z mass in electron-neutrino channel",50,0,500},"z_mass");
//	auto h_se_enu_events_zmass = se_enu_z_rec_selection.Histo1D({"Single Electron Z_mass_enu_Channel","Single Electron Z mass in electron-neutrino channel",50,0,500},"z_mass");
//	auto h_met_enu_events_zmass = met_enu_z_rec_selection.Histo1D({"MET Z_mass_enu_Channel","MET Z mass in electron-neutrino channel",50,0,500},"z_mass");

        auto h_events_zmass_canvas = new TCanvas("Z mass", "Z mass",10,10,900,900);

        h_d_enu_events_zmass->GetXaxis()->SetTitle("mass/GeVC^2");
        h_d_enu_events_zmass->GetYaxis()->SetTitle("Events");

        h_d_enu_events_zmass->SetFillColor(kBlack);
        h_ww_enu_events_zmass->SetFillColor(kRed);
        h_wz_enu_events_zmass->SetFillColor(kOrange);
        h_ttZ_enu_events_zmass->SetFillColor(kYellow);
        h_zz_enu_events_zmass->SetFillColor(kTeal);
//	h_se_enu_events_zmass->SetFillColor(kViolet);
//	h_met_enu_events_zmass->SetFillColor(kPink);

	h_d_enu_events_zmass->Scale(NWS_F);
        h_ww_enu_events_zmass->Scale(NWS_F);
        h_wz_enu_events_zmass->Scale(NWS_F);
        h_ttZ_enu_events_zmass->Scale(NWS_F);
        h_zz_enu_events_zmass->Scale(NWS_F);


	h_d_enu_events_zmass->Draw();
        h_ww_enu_events_zmass->Draw("SAME");
        h_wz_enu_events_zmass->Draw("SAME");
        h_ttZ_enu_events_zmass->Draw("SAME");
        h_zz_enu_events_zmass->Draw("SAME");
//	h_se_enu_events_zmass->Draw("SAME");
//	h_met_enu_events_zmass->Draw("SAME");

	h_events_zmass_canvas->BuildLegend();
        h_events_zmass_canvas->SaveAs("en_Z_mass.root");
	h_events_zmass_canvas->SaveAs("en_Z_mass.pdf");



        auto h_d_enu_events_jdphi = d_enu_z_rec_selection.Histo1D({"MC jdphi_enu_Channel","MC jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_ww_enu_events_jdphi = ww_enu_z_rec_selection.Histo1D({"WW jdphi_enu_Channel","WW jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_wz_enu_events_jdphi = wz_enu_z_rec_selection.Histo1D({"WZ jdphi_enu_Channel","WZ jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_ttZ_enu_events_jdphi = ttZ_enu_z_rec_selection.Histo1D({"ttZ jdphi_enu_Channel","ttZ jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_zz_enu_events_jdphi = zz_enu_z_rec_selection.Histo1D({"zz jdphi_enu_Channel","zz jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
//        auto h_se_enu_events_jdphi = se_enu_z_rec_selection.Histo1D({"Single Electron jdphi_enu_Channel","Single Electron jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");
//        auto h_met_enu_events_jdphi = met_enu_z_rec_selection.Histo1D({"MET jdphi_enu_Channel","MET jet deltaphi in electron-neutrino channel",50,0,5},"tight_jets_deltaphi");


        auto h_events_jdphi_canvas = new TCanvas("jetdeltaphi", "jetdeltaphi",10,10,900,900);

        h_d_enu_events_jdphi->GetXaxis()->SetTitle("jets delta phi / rad");
        h_d_enu_events_jdphi->GetYaxis()->SetTitle("Events");

        h_d_enu_events_jdphi->SetFillColor(kBlack);
        h_ww_enu_events_jdphi->SetFillColor(kRed);
        h_wz_enu_events_jdphi->SetFillColor(kOrange);
        h_ttZ_enu_events_jdphi->SetFillColor(kYellow);
        h_zz_enu_events_jdphi->SetFillColor(kTeal);
//        h_se_enu_events_jdphi->SetFillColor(kViolet);
//        h_met_enu_events_jdphi->SetFillColor(kPink);

        h_d_enu_events_jdphi->Scale(NWS_F);
        h_ww_enu_events_jdphi->Scale(NWS_F);
        h_wz_enu_events_jdphi->Scale(NWS_F);
        h_ttZ_enu_events_jdphi->Scale(NWS_F);
        h_zz_enu_events_jdphi->Scale(NWS_F);


        h_d_enu_events_jdphi->Draw();
        h_ww_enu_events_jdphi->Draw("SAME");
        h_wz_enu_events_jdphi->Draw("SAME");
        h_ttZ_enu_events_jdphi->Draw("SAME");
        h_zz_enu_events_jdphi->Draw("SAME");
//        h_se_enu_events_jdphi->Draw("SAME");
//        h_met_enu_events_jdphi->Draw("SAME");

        h_events_jdphi_canvas->BuildLegend();
        h_events_jdphi_canvas->SaveAs("en_jet_dphi.root");
        h_events_jdphi_canvas->SaveAs("en_jet_dphi.pdf");



        auto h_d_enu_events_zmetdphi = d_enu_z_rec_selection.Histo1D({"MC zmetdphi_enu_Channel","MC Z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_ww_enu_events_zmetdphi = ww_enu_z_rec_selection.Histo1D({"WW zmetdphi_enu_Channel","WW Z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_wz_enu_events_zmetdphi = wz_enu_z_rec_selection.Histo1D({"WZ zmetdphi_enu_Channel","WZ Z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_ttZ_enu_events_zmetdphi = ttZ_enu_z_rec_selection.Histo1D({"ttZ zmetdphi_enu_Channel","ttZ Z met pt jet deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_zz_enu_events_zmetdphi = zz_enu_z_rec_selection.Histo1D({"zz zmetdphi_enu_Channel","zz z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
//        auto h_se_enu_events_zmetdphi = se_enu_z_rec_selection.Histo1D({"Single Electron zmetdphi_enu_Channel","Single Electron z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");
//       auto h_met_enu_events_zmetdphi = met_enu_z_rec_selection.Histo1D({"MET zmetdphi_enu_Channel","MET z met pt deltaphi in electron-neutrino channel",50,0,5},"ZMet_deltaphi");



        auto h_events_zmetdphi_canvas = new TCanvas("zmetdeltaphi", "zmetdeltaphi",10,10,900,900);

        h_d_enu_events_zmetdphi->GetXaxis()->SetTitle("z and met delta phi / rad");
        h_d_enu_events_zmetdphi->GetYaxis()->SetTitle("Events");

        h_d_enu_events_zmetdphi->SetFillColor(kBlack);
        h_ww_enu_events_zmetdphi->SetFillColor(kRed);
        h_wz_enu_events_zmetdphi->SetFillColor(kOrange);
        h_ttZ_enu_events_zmetdphi->SetFillColor(kYellow);
        h_zz_enu_events_zmetdphi->SetFillColor(kTeal);
//        h_se_enu_events_zmetdphi->SetFillColor(kViolet);
//        h_met_enu_events_zmetdphi->SetFillColor(kPink);

        h_d_enu_events_zmetdphi->Scale(NWS_F);
        h_ww_enu_events_zmetdphi->Scale(NWS_F);
        h_wz_enu_events_zmetdphi->Scale(NWS_F);
        h_ttZ_enu_events_zmetdphi->Scale(NWS_F);
        h_zz_enu_events_zmetdphi->Scale(NWS_F);


        h_d_enu_events_zmetdphi->Draw();
        h_ww_enu_events_zmetdphi->Draw("SAME");
        h_wz_enu_events_zmetdphi->Draw("SAME");
        h_ttZ_enu_events_zmetdphi->Draw("SAME");
        h_zz_enu_events_zmetdphi->Draw("SAME");
//        h_se_enu_events_zmetdphi->Draw("SAME");
//        h_met_enu_events_zmetdphi->Draw("SAME");

        h_events_zmetdphi_canvas->BuildLegend();
        h_events_zmetdphi_canvas->SaveAs("en_zmetpt_dphi.root");
        h_events_zmetdphi_canvas->SaveAs("en_zmetpt_dphi.pdf");


        auto h_d_enu_events_zwdphi = d_enu_z_rec_selection.Histo1D({"MC zwdphi_enu_Channel","MC Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_ww_enu_events_zwdphi = ww_enu_z_rec_selection.Histo1D({"WW zwdphi_enu_Channel","WW Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_wz_enu_events_zwdphi = wz_enu_z_rec_selection.Histo1D({"WZ zwdphi_enu_Channel","WZ Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_ttZ_enu_events_zwdphi = ttZ_enu_z_rec_selection.Histo1D({"ttZ zwdphi_enu_Channel","ttZ  Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_zz_enu_events_zwdphi = zz_enu_z_rec_selection.Histo1D({"zz zwdphi_enu_Channel","zz z Z w deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
//        auto h_se_enu_events_zwdphi = se_enu_z_rec_selection.Histo1D({"Single Electron zmetdphi_enu_Channel","Single Electron z met pt deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");
//        auto h_met_enu_events_zwdphi = met_enu_z_rec_selection.Histo1D({"MET zmetdphi_enu_Channel","MET z met pt deltaphi in electron-neutrino channel",50,0,5},"ZW_deltaphi");


        auto h_events_zwdphi_canvas = new TCanvas("zwdeltaphi", "zwdeltaphi",10,10,900,900);

        h_d_enu_events_zwdphi->GetXaxis()->SetTitle("z and w delta phi / rad");
        h_d_enu_events_zwdphi->GetYaxis()->SetTitle("Events");

        h_d_enu_events_zwdphi->SetFillColor(kBlack);
        h_ww_enu_events_zwdphi->SetFillColor(kRed);
        h_wz_enu_events_zwdphi->SetFillColor(kOrange);
        h_ttZ_enu_events_zwdphi->SetFillColor(kYellow);
        h_zz_enu_events_zwdphi->SetFillColor(kTeal);
//        h_se_enu_events_zwdphi->SetFillColor(kViolet);
//        h_met_enu_events_zwdphi->SetFillColor(kPink);

        h_d_enu_events_zwdphi->Scale(NWS_F);
        h_ww_enu_events_zwdphi->Scale(NWS_F);
        h_wz_enu_events_zwdphi->Scale(NWS_F);
        h_ttZ_enu_events_zwdphi->Scale(NWS_F);
        h_zz_enu_events_zwdphi->Scale(NWS_F);


        h_d_enu_events_zwdphi->Draw();
        h_ww_enu_events_zwdphi->Draw("SAME");
        h_wz_enu_events_zwdphi->Draw("SAME");
        h_ttZ_enu_events_zwdphi->Draw("SAME");
        h_zz_enu_events_zmetdphi->Draw("SAME");
//        h_se_enu_events_zwdphi->Draw("SAME");
//        h_met_enu_events_zwdphi->Draw("SAME");

        h_events_zwdphi_canvas->BuildLegend();
        h_events_zwdphi_canvas->SaveAs("en_zwd_dphi.root");
        h_events_zwdphi_canvas->SaveAs("en_zw_dphi.pdf");




        auto h_d_enu_events_ejdr = d_enu_z_rec_selection.Histo1D({"MC ejdr_enu_Channel","MC e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
        auto h_ww_enu_events_ejdr = ww_enu_z_rec_selection.Histo1D({"WW ejdr_enu_Channel","WW e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
        auto h_wz_enu_events_ejdr = wz_enu_z_rec_selection.Histo1D({"WZ ejdr_enu_Channel","WZ e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
        auto h_ttZ_enu_events_ejdr = ttZ_enu_z_rec_selection.Histo1D({"ttZ ejdr_enu_Channel","ttZ  e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
        auto h_zz_enu_events_ejdr = zz_enu_z_rec_selection.Histo1D({"zz ejdr_enu_Channel","zz e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
//        auto h_se_enu_events_ejdr = se_enu_z_rec_selection.Histo1D({"Single Electron ejdr_enu_Channel","Single Electron e and jet pt deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");
//        auto h_met_enu_events_ejdr = met_enu_z_rec_selection.Histo1D({"MET ejdr_enu_Channel","MET e and jet deltaR in electron-neutrino channel",50,0,10},"jet_e_min_dR");

        auto h_events_ejdr_canvas = new TCanvas("ejdeltar", "ejdeltar",10,10,900,900);

        h_d_enu_events_ejdr->GetXaxis()->SetTitle("e and jet deltaR");
        h_d_enu_events_ejdr->GetYaxis()->SetTitle("Events");



        h_d_enu_events_ejdr->SetFillColor(kBlack);
        h_ww_enu_events_ejdr->SetFillColor(kRed);
        h_wz_enu_events_ejdr->SetFillColor(kOrange);
        h_ttZ_enu_events_ejdr->SetFillColor(kYellow);
        h_zz_enu_events_ejdr->SetFillColor(kTeal);
//        h_se_enu_events_ejdr->SetFillColor(kViolet);
//        h_met_enu_events_ejdr->SetFillColor(kPink);

        h_d_enu_events_ejdr->Scale(NWS_F);
        h_ww_enu_events_ejdr->Scale(NWS_F);
        h_wz_enu_events_ejdr->Scale(NWS_F);
        h_ttZ_enu_events_ejdr->Scale(NWS_F);
        h_zz_enu_events_ejdr->Scale(NWS_F);


        h_d_enu_events_ejdr->Draw();
        h_ww_enu_events_ejdr->Draw("SAME");
        h_wz_enu_events_ejdr->Draw("SAME");
        h_ttZ_enu_events_ejdr->Draw("SAME");
        h_zz_enu_events_ejdr->Draw("SAME");
//        h_se_enu_events_ejdr->Draw("SAME");
//        h_met_enu_events_ejdr->Draw("SAME");

        h_events_ejdr_canvas->BuildLegend();
        h_events_ejdr_canvas->SaveAs("en_ej_dr.root");
        h_events_ejdr_canvas->SaveAs("en_ej_dr.pdf");



        auto h_d_enu_events_ezdr = d_enu_z_rec_selection.Histo1D({"MC ezdr_enu_Channel","MC e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
        auto h_ww_enu_events_ezdr = ww_enu_z_rec_selection.Histo1D({"WW ezdr_enu_Channel","WW e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
        auto h_wz_enu_events_ezdr = wz_enu_z_rec_selection.Histo1D({"WZ ezdr_enu_Channel","WZ e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
        auto h_ttZ_enu_events_ezdr = ttZ_enu_z_rec_selection.Histo1D({"ttZ ezdr_enu_Channel","ttZ  e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
        auto h_zz_enu_events_ezdr = zz_enu_z_rec_selection.Histo1D({"zz ezdr_enu_Channel","zz e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
//        auto h_se_enu_events_ezdr = se_enu_z_rec_selection.Histo1D({"Single Electron ezdr_enu_Channel","Single Electron e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");
//        auto h_met_enu_events_ezdr = met_enu_z_rec_selection.Histo1D({"MET ezdr_enu_Channel","MET e and z deltaR in electron-neutrino channel",50,0,10},"z_e_min_dR");

        auto h_events_ezdr_canvas = new TCanvas("ezdeltar", "ezdeltar",10,10,900,900);

        h_d_enu_events_ezdr->GetXaxis()->SetTitle("e and z deltaR");
        h_d_enu_events_ezdr->GetYaxis()->SetTitle("Events");


        h_d_enu_events_ezdr->SetFillColor(kBlack);
        h_ww_enu_events_ezdr->SetFillColor(kRed);
        h_wz_enu_events_ezdr->SetFillColor(kOrange);
        h_ttZ_enu_events_ezdr->SetFillColor(kYellow);
        h_zz_enu_events_ezdr->SetFillColor(kTeal);
//        h_se_enu_events_ezdr->SetFillColor(kViolet);
//        h_met_enu_events_ezdr->SetFillColor(kPink);


        h_d_enu_events_ezdr->Scale(NWS_F);
        h_ww_enu_events_ezdr->Scale(NWS_F);
        h_wz_enu_events_ezdr->Scale(NWS_F);
        h_ttZ_enu_events_ezdr->Scale(NWS_F);
        h_zz_enu_events_ezdr->Scale(NWS_F);


        h_d_enu_events_ezdr->Draw();
        h_ww_enu_events_ezdr->Draw("SAME");
        h_wz_enu_events_ezdr->Draw("SAME");
        h_ttZ_enu_events_ezdr->Draw("SAME");
        h_zz_enu_events_ezdr->Draw("SAME");
//        h_se_enu_events_ezdr->Draw("SAME");
//        h_met_enu_events_ezdr->Draw("SAME");

        h_events_ezdr_canvas->BuildLegend();
        h_events_ezdr_canvas->SaveAs("en_ez_dr.root");
        h_events_ezdr_canvas->SaveAs("en_ez_dr.pdf");





//////////////////////////////////////////////////////////////////////// Mu-Nu MCs //////////////////////////////////////////////////////////////////////////////////
        auto h_d_munu_events_mupt = d_munu_top_selection.Histo1D({"MC muon_pt_munu_Channel","MC muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
        auto h_ww_munu_events_mupt = ww_munu_top_selection.Histo1D({"WW muon_pt_munu_Channel","WW muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
        auto h_wz_munu_events_mupt = wz_munu_top_selection.Histo1D({"WZ muon_pt_munu_Channel","WZ muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
        auto h_zz_munu_events_mupt = zz_munu_top_selection.Histo1D({"ZZ muon_pt_munu_Channel","ZZ muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
        auto h_ttZ_munu_events_mupt = ttZ_munu_top_selection.Histo1D({"ttZ muon_pt_munu_Channel","ttZ muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
//	auto h_smu_munu_events_mupt = sm_munu_top_selection.Histo1D({"Single Muon muon_pt_munu_Channel","Single muon muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");
//	auto h_met_munu_events_mupt = met_munu_top_selection.Histo1D({"MET muon_pt_munu_Channel","MET muon pt in muon-neutrino channel",50,0,300},"tight_mu_pt");

        auto h_events_mupt_canvas = new TCanvas("muon pt", "muon pt",10,10,900,900);

	h_d_munu_events_mupt->GetXaxis()->SetTitle("Pt/GeV");
        h_d_munu_events_mupt->GetYaxis()->SetTitle("Events");

        h_d_munu_events_mupt->SetFillColor(kBlack);
        h_ww_munu_events_mupt->SetFillColor(kRed);
        h_wz_munu_events_mupt->SetFillColor(kOrange);
        h_ttZ_munu_events_mupt->SetFillColor(kYellow);
        h_zz_munu_events_mupt->SetFillColor(kTeal);
//	h_smu_munu_events_mupt->SetFillColor(kViolet);
//	h_met_munu_events_mupt->SetFillColor(kPink);

	h_d_munu_events_mupt->Scale(NWS_F);
	h_ww_munu_events_mupt->Scale(NWS_F);
	h_wz_munu_events_mupt->Scale(NWS_F);
	h_ttZ_munu_events_mupt->Scale(NWS_F);
	h_zz_munu_events_mupt->Scale(NWS_F);


	h_d_munu_events_mupt->Draw();
        h_ww_munu_events_mupt->Draw("SAME");
        h_wz_munu_events_mupt->Draw("SAME");
        h_ttZ_munu_events_mupt->Draw("SAME");
        h_zz_munu_events_mupt->Draw("SAME");
//	h_smu_munu_events_mupt->Draw("SAME");
//	h_met_munu_events_mupt->Draw("SAME");

        h_events_mupt_canvas->BuildLegend();
	h_events_mupt_canvas->SaveAs("munu_pt.root");
	h_events_mupt_canvas->SaveAs("munu_pt.pdf");


	auto h_d_munu_events_jpt = d_munu_top_selection.Histo1D({"MC jet_pt_munu_Channel","MC jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_ww_munu_events_jpt = ww_munu_top_selection.Histo1D({"WW jet_pt_munu_Channel","WW jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_wz_munu_events_jpt = wz_munu_top_selection.Histo1D({"WZ jet_pt_munu_Channel","WZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_zz_munu_events_jpt = zz_munu_top_selection.Histo1D({"ZZ jet_pt_munu_Channel","ZZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
        auto h_ttZ_munu_events_jpt = ttZ_munu_top_selection.Histo1D({"ttZ jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
//	auto h_smu_munu_events_jpt = sm_munu_top_selection.Histo1D({"Single muon jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");
//	auto h_met_munu_events_jpt = met_munu_top_selection.Histo1D({"MET jet_pt_munu_Channel","ttZ jet pt in muon-neutrino channel",50,0,300},"tight_jets_pt");

	auto h_munu_events_jetpt_canvas = new TCanvas("mu nu jet pt ", "munu jet pt",10,10,900,900);

        h_d_munu_events_jpt->GetXaxis()->SetTitle("pt / GeV");
        h_d_munu_events_jpt->GetYaxis()->SetTitle("Events");
        h_d_munu_events_jpt->SetFillColor(kBlack);
        h_ww_munu_events_jpt->SetFillColor(kRed);
        h_wz_munu_events_jpt->SetFillColor(kOrange);
        h_ttZ_munu_events_jpt->SetFillColor(kYellow);
        h_zz_munu_events_jpt->SetFillColor(kTeal);
//        h_smu_munu_events_jpt->SetFillColor(kViolet);
//        h_met_munu_events_jpt->SetFillColor(kPink);

        h_d_munu_events_jpt->Scale(NWS_F);
        h_ww_munu_events_jpt->Scale(NWS_F);
        h_wz_munu_events_jpt->Scale(NWS_F);
        h_ttZ_munu_events_jpt->Scale(NWS_F);
        h_zz_munu_events_jpt->Scale(NWS_F);



        h_d_munu_events_jpt->Draw();
        h_ww_munu_events_jpt->Draw("SAME");
        h_wz_munu_events_jpt->Draw("SAME");
        h_ttZ_munu_events_jpt->Draw("SAME");
        h_zz_munu_events_jpt->Draw("SAME");
//        h_smu_munu_events_jpt->Draw("SAME");
//        h_met_munu_events_jpt->Draw("SAME");

        h_munu_events_jetpt_canvas->BuildLegend();
        h_munu_events_jetpt_canvas->SaveAs("munu_jetpt.root");
        h_munu_events_jetpt_canvas->SaveAs("munu_jetpt.pdf");




        auto h_d_munu_events_mueta = d_munu_top_selection.Histo1D({"MC muon_eta_enu_Channel","MC muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
        auto h_ww_munu_events_mueta = ww_munu_top_selection.Histo1D({"WW muon_eta_enu_Channel","WW muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
        auto h_wz_munu_events_mueta = wz_munu_top_selection.Histo1D({"WZ muon_eta_enu_Channel","WZ muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
        auto h_zz_munu_events_mueta = zz_munu_top_selection.Histo1D({"ZZ muon_eta_enu_Channel","ZZ muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
        auto h_ttZ_munu_events_mueta = ttZ_munu_top_selection.Histo1D({"ttZ muon_eta_enu_Channel","ttZ muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
//        auto h_smu_munu_events_mueta = sm_munu_top_selection.Histo1D({"Single Muon muon_eta_Channel","Single Muon Muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");
//        auto h_met_munu_events_mueta = met_munu_top_selection.Histo1D({"MET Muon_eta_Channel","MET muon eta in muon-neutrino channel",50,-4,4},"tight_mu_eta");


        auto h_events_mueta_canvas = new TCanvas("Muon eta", "Muon eta",10,10,900,900);

        h_d_munu_events_mueta->GetXaxis()->SetTitle("eta");
        h_d_munu_events_mueta->GetYaxis()->SetTitle("Events");

        h_d_munu_events_mueta->SetFillColor(kBlack);
        h_ww_munu_events_mueta->SetFillColor(kRed);
        h_wz_munu_events_mueta->SetFillColor(kOrange);
        h_ttZ_munu_events_mueta->SetFillColor(kYellow);
        h_zz_munu_events_mueta->SetFillColor(kTeal);
//        h_smu_munu_events_mueta->SetFillColor(kViolet);
//        h_met_munu_events_mueta->SetFillColor(kPink);

        h_d_munu_events_mueta->Scale(NWS_F);
        h_ww_munu_events_mueta->Scale(NWS_F);
        h_wz_munu_events_mueta->Scale(NWS_F);
        h_ttZ_munu_events_mueta->Scale(NWS_F);
        h_zz_munu_events_mueta->Scale(NWS_F);


        h_d_munu_events_mueta->Draw();
        h_ww_munu_events_mueta->Draw("SAME");
        h_wz_munu_events_mueta->Draw("SAME");
        h_ttZ_munu_events_mueta->Draw("SAME");
        h_zz_munu_events_mueta->Draw("SAME");
//        h_smu_munu_events_mueta->Draw("SAME");
//        h_met_munu_events_mueta->Draw("SAME");

        h_events_mueta_canvas->BuildLegend();
        h_events_mueta_canvas->SaveAs("munu_eta.root");
        h_events_mueta_canvas->SaveAs("munu_eta.pdf");



        auto h_d_munu_events_jeta = d_munu_top_selection.Histo1D({"MC jet_eta_munu_Channel","MC jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_ww_munu_events_jeta = ww_munu_top_selection.Histo1D({"WW jet_eta_munu_Channel","WW jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_wz_munu_events_jeta = wz_munu_top_selection.Histo1D({"WZ jet_eta_munu_Channel","WZ jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_zz_munu_events_jeta = zz_munu_top_selection.Histo1D({"ZZ jet_eta_munu_Channel","ZZ jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
        auto h_ttZ_munu_events_jeta = ttZ_munu_top_selection.Histo1D({"ttZ jet_eta_munu_Channel","ttZ jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");
//        auto h_smu_munu_events_jeta = sm_munu_top_selection.Histo1D({"Single Muon jet_eta_munu_Channel","Single Muon jet eta in Muon-neutrino channels",50,-4,4},"tight_jets_eta");
//        auto h_met_munu_events_jeta = met_munu_top_selection.Histo1D({"MET muon_jet_munu_Channel","MET jet eta in Muon-neutrino channel",50,-4,4},"tight_jets_eta");


        auto h_munu_events_jeta_canvas = new TCanvas("jet eta", "jet eta",10,10,900,900);

        h_d_munu_events_jeta->GetXaxis()->SetTitle("eta");
        h_d_munu_events_jeta->GetYaxis()->SetTitle("Events");

        h_d_munu_events_jeta->SetFillColor(kBlack);
        h_ww_munu_events_jeta->SetFillColor(kRed);
        h_wz_munu_events_jeta->SetFillColor(kOrange);
        h_ttZ_munu_events_jeta->SetFillColor(kYellow);
        h_zz_munu_events_jeta->SetFillColor(kTeal);
//        h_smu_munu_events_jeta->SetFillColor(kViolet);
//        h_met_munu_events_jeta->SetFillColor(kPink);


        h_d_munu_events_jeta->Scale(NWS_F);
        h_ww_munu_events_jeta->Scale(NWS_F);
        h_wz_munu_events_jeta->Scale(NWS_F);
        h_ttZ_munu_events_jeta->Scale(NWS_F);
        h_zz_munu_events_jeta->Scale(NWS_F);


        h_d_munu_events_jeta->Draw();
        h_ww_munu_events_jeta->Draw("SAME");
        h_wz_munu_events_jeta->Draw("SAME");
        h_ttZ_munu_events_jeta->Draw("SAME");
        h_zz_munu_events_jeta->Draw("SAME");
//        h_smu_munu_events_jeta->Draw("SAME");
//        h_met_munu_events_jeta->Draw("SAME");

        h_munu_events_jeta_canvas->BuildLegend();
        h_munu_events_jeta_canvas->SaveAs("munu_jeta.root");
        h_munu_events_jeta_canvas->SaveAs("munu_jeta.pdf");



        auto h_d_munu_events_wmass = d_munu_top_selection.Histo1D({"MC munu_w_mass","MC Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
        auto h_ww_munu_events_wmass = ww_munu_top_selection.Histo1D({"WW munu_w_mass","WW Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
        auto h_wz_munu_events_wmass = wz_munu_top_selection.Histo1D({"WZ munu_w_mass","WZ Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
        auto h_ttZ_munu_events_wmass = ttZ_munu_top_selection.Histo1D({"ttZ munu_w_mass","ttZ Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
        auto h_zz_munu_events_wmass = zz_munu_top_selection.Histo1D({"ZZ munu_w_mass","ZZ Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
//        auto h_smu_munu_events_wmass = sm_munu_top_selection.Histo1D({"Single Muon munu_w_mass","Single Muon Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");
//        auto h_met_munu_events_wmass = met_munu_top_selection.Histo1D({"MET munu_w_mass","MET Muon-neutrino transverse w mass",50,0,500},"w_mu_mass");


        auto h_munu_events_wmass_canvas = new TCanvas("munu w mass ", "munu w mass",10,10,900,900);

        h_d_munu_events_wmass->GetXaxis()->SetTitle("mass/GeV/C^2");
        h_d_munu_events_wmass->GetYaxis()->SetTitle("Events");
        h_d_munu_events_wmass->SetFillColor(kBlack);
        h_ww_munu_events_wmass->SetFillColor(kRed);
        h_wz_munu_events_wmass->SetFillColor(kOrange);
        h_ttZ_munu_events_wmass->SetFillColor(kYellow);
        h_zz_munu_events_wmass->SetFillColor(kTeal);
//        h_smu_munu_events_wmass->SetFillColor(kViolet);
//        h_met_munu_events_wmass->SetFillColor(kPink);

        h_d_munu_events_wmass->Scale(NWS_F);
        h_ww_munu_events_wmass->Scale(NWS_F);
        h_wz_munu_events_wmass->Scale(NWS_F);
        h_ttZ_munu_events_wmass->Scale(NWS_F);
        h_zz_munu_events_wmass->Scale(NWS_F);


        h_d_munu_events_wmass->Draw();
        h_ww_munu_events_wmass->Draw("SAME");
        h_wz_munu_events_wmass->Draw("SAME");
        h_ttZ_munu_events_wmass->Draw("SAME");
        h_zz_munu_events_wmass->Draw("SAME");
//        h_smu_munu_events_wmass->Draw("SAME");
//        h_met_munu_events_wmass->Draw("SAME");

        h_munu_events_wmass_canvas->BuildLegend();
        h_munu_events_wmass_canvas->SaveAs("munu_transverse_Wmass.root");
        h_munu_events_wmass_canvas->SaveAs("munu_transverse_Wmass.pdf");




        auto h_d_munu_events_zmass = d_munu_z_rec_selection.Histo1D({"MC Z_mass_munu_Channel","MC Z mass in Muon-neutrino channel",50,0,500},"z_mass");
        auto h_ww_munu_events_zmass = ww_munu_z_rec_selection.Histo1D({"WW Z_mass_munu_Channel","WW Z mass in Muon-neutrino channel",50,0,500},"z_mass");
        auto h_wz_munu_events_zmass = wz_munu_z_rec_selection.Histo1D({"WZ Z_mass_munu_Channel","WZ Z mass in Muon-neutrino channel",50,0,500},"z_mass");
        auto h_ttZ_munu_events_zmass = ttZ_munu_z_rec_selection.Histo1D({"ttZ Z_mass_munu_Channel","ttZ Z mass in Muon-neutrino channel",50,0,500},"z_mass");
        auto h_zz_munu_events_zmass = zz_munu_z_rec_selection.Histo1D({"zz Z_mass_munu_Channel","zz Z mass in Muon-neutrino channel",50,0,500},"z_mass");
//	auto h_smu_munu_events_zmass = sm_munu_z_rec_selection.Histo1D({"Single Muon Z_mass_enu_Channel","Single Muon Z mass in Muon-neutrino channel",50,0,500},"z_mass");
//	auto h_met_munu_events_zmass = met_munu_z_rec_selection.Histo1D({"MET Z_mass_munu_Channel","MET Z mass in Muon-neutrino channel",50,0,500},"z_mass");

        auto h_munu_events_zmass_canvas = new TCanvas("Z mass", "Z mass",10,10,900,900);

        h_d_munu_events_zmass->GetXaxis()->SetTitle("mass/GeVC^2");
        h_d_munu_events_zmass->GetYaxis()->SetTitle("Events");

        h_d_munu_events_zmass->SetFillColor(kBlack);
        h_ww_munu_events_zmass->SetFillColor(kRed);
        h_wz_munu_events_zmass->SetFillColor(kOrange);
        h_ttZ_munu_events_zmass->SetFillColor(kYellow);
        h_zz_munu_events_zmass->SetFillColor(kTeal);
//	h_smu_munu_events_zmass->SetFillColor(kViolet);
//	h_met_munu_events_zmass->SetFillColor(kPink);

	h_d_munu_events_zmass->Scale(NWS_F);
        h_ww_munu_events_zmass->Scale(NWS_F);
        h_wz_munu_events_zmass->Scale(NWS_F);
        h_ttZ_munu_events_zmass->Scale(NWS_F);
        h_zz_munu_events_zmass->Scale(NWS_F);


	h_d_munu_events_zmass->Draw();
        h_ww_munu_events_zmass->Draw("SAME");
        h_wz_munu_events_zmass->Draw("SAME");
        h_ttZ_munu_events_zmass->Draw("SAME");
        h_zz_munu_events_zmass->Draw("SAME");
//	h_smu_munu_events_zmass->Draw("SAME");
//	h_met_munu_events_zmass->Draw("SAME");

	h_munu_events_zmass_canvas->BuildLegend();
        h_munu_events_zmass_canvas->SaveAs("munu_Z_mass.root");
	h_munu_events_zmass_canvas->SaveAs("munu_Z_mass.pdf");



        auto h_d_munu_events_jdphi = d_munu_z_rec_selection.Histo1D({"MC jdphi_munu_Channel","MC jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_ww_munu_events_jdphi = ww_munu_z_rec_selection.Histo1D({"WW jdphi_munu_Channel","WW jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_wz_munu_events_jdphi = wz_munu_z_rec_selection.Histo1D({"WZ jdphi_munu_Channel","WZ jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_ttZ_munu_events_jdphi = ttZ_munu_z_rec_selection.Histo1D({"ttZ jdphi_munu_Channel","ttZ jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
        auto h_zz_munu_events_jdphi = zz_munu_z_rec_selection.Histo1D({"zz jdphi_munu_Channel","zz jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
//        auto h_smu_munu_events_jdphi = sm_munu_z_rec_selection.Histo1D({"Single Muon jdphi_munu_Channel","Single Muon jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");
//        auto h_met_munu_events_jdphi = met_munu_z_rec_selection.Histo1D({"MET jdphi_munu_Channel","MET jet deltaphi in Muon-neutrino channel",50,0,5},"tight_jets_deltaphi");


        auto h_munu_events_jdphi_canvas = new TCanvas("jetdeltaphi", "jetdeltaphi",10,10,900,900);

        h_d_munu_events_jdphi->GetXaxis()->SetTitle("jets delta phi / rad");
        h_d_munu_events_jdphi->GetYaxis()->SetTitle("Events");

        h_d_munu_events_jdphi->SetFillColor(kBlack);
        h_ww_munu_events_jdphi->SetFillColor(kRed);
        h_wz_munu_events_jdphi->SetFillColor(kOrange);
        h_ttZ_munu_events_jdphi->SetFillColor(kYellow);
        h_zz_munu_events_jdphi->SetFillColor(kTeal);
//        h_smu_munu_events_jdphi->SetFillColor(kViolet);
//        h_met_munu_events_jdphi->SetFillColor(kPink);

        h_d_munu_events_jdphi->Scale(NWS_F);
        h_ww_munu_events_jdphi->Scale(NWS_F);
        h_wz_munu_events_jdphi->Scale(NWS_F);
        h_ttZ_munu_events_jdphi->Scale(NWS_F);
        h_zz_munu_events_jdphi->Scale(NWS_F);


        h_d_munu_events_jdphi->Draw();
        h_ww_munu_events_jdphi->Draw("SAME");
        h_wz_munu_events_jdphi->Draw("SAME");
        h_ttZ_munu_events_jdphi->Draw("SAME");
        h_zz_munu_events_jdphi->Draw("SAME");
//        h_smu_munu_events_jdphi->Draw("SAME");
//        h_met_munu_events_jdphi->Draw("SAME");

        h_munu_events_jdphi_canvas->BuildLegend();
        h_munu_events_jdphi_canvas->SaveAs("munu_jet_dphi.root");
        h_munu_events_jdphi_canvas->SaveAs("munu_jet_dphi.pdf");



        auto h_d_munu_events_zmetdphi = d_munu_z_rec_selection.Histo1D({"MC zmetdphi_munu_Channel","MC Z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_ww_munu_events_zmetdphi = ww_munu_z_rec_selection.Histo1D({"WW zmetdphi_munu_Channel","WW Z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_wz_munu_events_zmetdphi = wz_munu_z_rec_selection.Histo1D({"WZ zmetdphi_munu_Channel","WZ Z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_ttZ_munu_events_zmetdphi = ttZ_munu_z_rec_selection.Histo1D({"ttZ zmetdphi_munu_Channel","ttZ Z met pt jet deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
        auto h_zz_munu_events_zmetdphi = zz_munu_z_rec_selection.Histo1D({"zz zmetdphi_munu_Channel","zz z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
//        auto h_smu_munu_events_zmetdphi = sm_munu_z_rec_selection.Histo1D({"Single Muon zmetdphi_munu_Channel","Single Muon z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");
//        auto h_met_munu_events_zmetdphi = met_munu_z_rec_selection.Histo1D({"MET zmetdphi_munu_Channel","MET z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZMet_deltaphi");



        auto h_munu_events_zmetdphi_canvas = new TCanvas("zmetdeltaphi", "zmetdeltaphi",10,10,900,900);

        h_d_munu_events_zmetdphi->GetXaxis()->SetTitle("z and met delta phi / rad");
        h_d_munu_events_zmetdphi->GetYaxis()->SetTitle("Events");

        h_d_munu_events_zmetdphi->SetFillColor(kBlack);
        h_ww_munu_events_zmetdphi->SetFillColor(kRed);
        h_wz_munu_events_zmetdphi->SetFillColor(kOrange);
        h_ttZ_munu_events_zmetdphi->SetFillColor(kYellow);
        h_zz_munu_events_zmetdphi->SetFillColor(kTeal);
//        h_smu_munu_events_zmetdphi->SetFillColor(kViolet);
//        h_met_munu_events_zmetdphi->SetFillColor(kPink);

        h_d_munu_events_zmetdphi->Scale(NWS_F);
        h_ww_munu_events_zmetdphi->Scale(NWS_F);
        h_wz_munu_events_zmetdphi->Scale(NWS_F);
        h_ttZ_munu_events_zmetdphi->Scale(NWS_F);
        h_zz_munu_events_zmetdphi->Scale(NWS_F);


        h_d_munu_events_zmetdphi->Draw();
        h_ww_munu_events_zmetdphi->Draw("SAME");
        h_wz_munu_events_zmetdphi->Draw("SAME");
        h_ttZ_munu_events_zmetdphi->Draw("SAME");
        h_zz_munu_events_zmetdphi->Draw("SAME");
//        h_smu_munu_events_zmetdphi->Draw("SAME");
//        h_met_munu_events_zmetdphi->Draw("SAME");

        h_munu_events_zmetdphi_canvas->BuildLegend();
        h_munu_events_zmetdphi_canvas->SaveAs("munu_zmetpt_dphi.root");
        h_munu_events_zmetdphi_canvas->SaveAs("munu_zmetpt_dphi.pdf");


        auto h_d_munu_events_zwdphi = d_munu_z_rec_selection.Histo1D({"MC zwdphi_munu_Channel","MC Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_ww_munu_events_zwdphi = ww_munu_z_rec_selection.Histo1D({"WW zwdphi_munu_Channel","WW Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_wz_munu_events_zwdphi = wz_munu_z_rec_selection.Histo1D({"WZ zwdphi_munu_Channel","WZ Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_ttZ_munu_events_zwdphi = ttZ_munu_z_rec_selection.Histo1D({"ttZ zwdphi_munu_Channel","ttZ  Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
        auto h_zz_munu_events_zwdphi = zz_munu_z_rec_selection.Histo1D({"zz zwdphi_munu_Channel","zz z Z w deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
//        auto h_smu_munu_events_zwdphi = sm_munu_z_rec_selection.Histo1D({"Single Muon zmetdphi_munu_Channel","Single Muon z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");
//        auto h_met_munu_events_zwdphi = met_munu_z_rec_selection.Histo1D({"MET zmetdphi_munu_Channel","MET z met pt deltaphi in Muon-neutrino channel",50,0,5},"ZW_deltaphi");


        auto h_munu_events_zwdphi_canvas = new TCanvas("zwdeltaphi", "zwdeltaphi",10,10,900,900);

        h_d_munu_events_zwdphi->GetXaxis()->SetTitle("z and w delta phi / rad");
        h_d_munu_events_zwdphi->GetYaxis()->SetTitle("Events");

        h_d_munu_events_zwdphi->SetFillColor(kBlack);
        h_ww_munu_events_zwdphi->SetFillColor(kRed);
        h_wz_munu_events_zwdphi->SetFillColor(kOrange);
        h_ttZ_munu_events_zwdphi->SetFillColor(kYellow);
        h_zz_munu_events_zwdphi->SetFillColor(kTeal);
//        h_smu_munu_events_zwdphi->SetFillColor(kViolet);
//        h_met_munu_events_zwdphi->SetFillColor(kPink);

        h_d_munu_events_zwdphi->Scale(NWS_F);
        h_ww_munu_events_zwdphi->Scale(NWS_F);
        h_wz_munu_events_zwdphi->Scale(NWS_F);
        h_ttZ_munu_events_zwdphi->Scale(NWS_F);
        h_zz_munu_events_zwdphi->Scale(NWS_F);


        h_d_munu_events_zwdphi->Draw();
        h_ww_munu_events_zwdphi->Draw("SAME");
        h_wz_munu_events_zwdphi->Draw("SAME");
        h_ttZ_munu_events_zwdphi->Draw("SAME");
        h_zz_munu_events_zmetdphi->Draw("SAME");
//        h_smu_munu_events_zwdphi->Draw("SAME");
//        h_met_munu_events_zwdphi->Draw("SAME");

        h_munu_events_zwdphi_canvas->BuildLegend();
        h_munu_events_zwdphi_canvas->SaveAs("munu_zwd_dphi.root");
        h_munu_events_zwdphi_canvas->SaveAs("munu_zw_dphi.pdf");




        auto h_d_munu_events_mujdr = d_munu_z_rec_selection.Histo1D({"MC ejdr_munu_Channel","MC muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
        auto h_ww_munu_events_mujdr = ww_munu_z_rec_selection.Histo1D({"WW ejdr_munu_Channel","WW muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
        auto h_wz_munu_events_mujdr = wz_munu_z_rec_selection.Histo1D({"WZ ejdr_munu_Channel","WZ muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
        auto h_ttZ_munu_events_mujdr = ttZ_munu_z_rec_selection.Histo1D({"ttZ ejdr_munu_Channel","ttZ muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
        auto h_zz_munu_events_mujdr = zz_munu_z_rec_selection.Histo1D({"zz ejdr_munu_Channel","zz muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");
//        auto h_smu_munu_events_mujdr = sm_munu_z_rec_selection.Histo1D({"Single Muon mujdr_munu_Channel","Single Muon muon and jet pt deltaR in Muon-neutrino channel",50,0,10},"jet_mu_min_dR");
//        auto h_met_munu_events_mujdr = met_munu_z_rec_selection.Histo1D({"MET mujdr_munu_Channel","MET muon and jet deltaR in muon-neutrino channel",50,0,10},"jet_mu_min_dR");

        auto h_events_mujdr_canvas = new TCanvas("mujdeltar", "mujdeltar",10,10,900,900);

        h_d_munu_events_mujdr->GetXaxis()->SetTitle("mu and jet deltaR");
        h_d_munu_events_mujdr->GetYaxis()->SetTitle("Events");



        h_d_munu_events_mujdr->SetFillColor(kBlack);
        h_ww_munu_events_mujdr->SetFillColor(kRed);
        h_wz_munu_events_mujdr->SetFillColor(kOrange);
        h_ttZ_munu_events_mujdr->SetFillColor(kYellow);
        h_zz_munu_events_mujdr->SetFillColor(kTeal);
//        h_smu_munu_events_mujdr->SetFillColor(kViolet);
//        h_met_munu_events_mujdr->SetFillColor(kPink);

        h_d_munu_events_mujdr->Scale(NWS_F);
        h_ww_munu_events_mujdr->Scale(NWS_F);
        h_wz_munu_events_mujdr->Scale(NWS_F);
        h_ttZ_munu_events_mujdr->Scale(NWS_F);
        h_zz_munu_events_mujdr->Scale(NWS_F);


        h_d_munu_events_mujdr->Draw();
        h_ww_munu_events_mujdr->Draw("SAME");
        h_wz_munu_events_mujdr->Draw("SAME");
        h_ttZ_munu_events_mujdr->Draw("SAME");
        h_zz_munu_events_mujdr->Draw("SAME");
//        h_smu_munu_events_mujdr->Draw("SAME");
//        h_met_munu_events_mujdr->Draw("SAME");

        h_events_mujdr_canvas->BuildLegend();
        h_events_mujdr_canvas->SaveAs("munu_muj_dr.root");
        h_events_mujdr_canvas->SaveAs("munu_muj_dr.pdf");



        auto h_d_munu_events_muzdr = d_munu_z_rec_selection.Histo1D({"MC muzdr_munu_Channel","MC mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
        auto h_ww_munu_events_muzdr = ww_munu_z_rec_selection.Histo1D({"WW muzdr_munu_Channel","WW mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
        auto h_wz_munu_events_muzdr = wz_munu_z_rec_selection.Histo1D({"WZ muzdr_munu_Channel","WZ mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
        auto h_ttZ_munu_events_muzdr = ttZ_munu_z_rec_selection.Histo1D({"ttZ muzdr_munu_Channel","ttZ  mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
        auto h_zz_munu_events_muzdr = zz_munu_z_rec_selection.Histo1D({"zz muzdr_munu_Channel","zz mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
//        auto h_smu_munu_events_muzdr = sm_munu_z_rec_selection.Histo1D({"Single Muon muzdr_munu_Channel","Single Muon mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");
//        auto h_met_munu_events_muzdr = met_munu_z_rec_selection.Histo1D({"MET muzdr_munu_Channel","MET mu and z deltaR in muon-neutrino channel",50,0,10},"z_mu_min_dR");

        auto h_events_muzdr_canvas = new TCanvas("muzdeltar", "muzdeltar",10,10,900,900);

        h_d_munu_events_muzdr->GetXaxis()->SetTitle("muon and z deltaR");
        h_d_munu_events_muzdr->GetYaxis()->SetTitle("Events");


        h_d_munu_events_muzdr->SetFillColor(kBlack);
        h_ww_munu_events_muzdr->SetFillColor(kRed);
        h_wz_munu_events_muzdr->SetFillColor(kOrange);
        h_ttZ_munu_events_muzdr->SetFillColor(kYellow);
        h_zz_munu_events_muzdr->SetFillColor(kTeal);
//        h_smu_munu_events_muzdr->SetFillColor(kViolet);
//        h_met_munu_events_muzdr->SetFillColor(kPink);


        h_d_munu_events_muzdr->Scale(NWS_F);
        h_ww_munu_events_muzdr->Scale(NWS_F);
        h_wz_munu_events_muzdr->Scale(NWS_F);
        h_ttZ_munu_events_muzdr->Scale(NWS_F);
        h_zz_munu_events_muzdr->Scale(NWS_F);


        h_d_munu_events_muzdr->Draw();
        h_ww_munu_events_muzdr->Draw("SAME");
        h_wz_munu_events_muzdr->Draw("SAME");
        h_ttZ_munu_events_muzdr->Draw("SAME");
        h_zz_munu_events_muzdr->Draw("SAME");
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
	electron_pT->SetFillColor(kRed);
        electron_pTCanvas->BuildLegend();
        electron_pT->Draw();
	electron_pTCanvas->SaveAs("Pure_electron_pT.root");


	auto electron_eta = d.Histo1D({"Electron eta","Electron_eta",20,-5,5},"Electron_eta");
        auto electron_etaCanvas = new TCanvas("Electron eta", "Electron_eta",10,10,700,700);
        electron_eta->GetXaxis()->SetTitle("eta");
        electron_eta->GetYaxis()->SetTitle("Events");
        electron_eta->SetFillColor(kRed);
        electron_etaCanvas->BuildLegend();
        electron_eta->Draw();
        electron_etaCanvas->SaveAs("Pure_electron_eta.root");
*/
/////////////////////////////////////////////////////////////////  Z Recon from U Histo2D /////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TO DO
// you need correlation of delta phi vs z mass reconstruction.
//later on W reconstruction : and lepton and MET , delta phi vs w recons, real data : reconstruct trans. mass

///////////////////////////////////////////////////////// THSTACK IN PROGRESS //////////////////////////////////////////////////////////////////////////////////////



	//auto NumDist = new THStack("NumDist","Quark Dist. Per event");//defining a stack that is made up many histograms
   	//auto temp_bnumHist = (TH1D*)&bnumHist.GetValue();//need to copy bnumHist on another hist to make it work with thstack
	//temp_bnumHist->SetFillColor(kRed); // saying to stack use colour red for the histogram
	//bnumHist->SetFillColor(kRed);
	//NumDist->Add(TH1 &bnumHist);
   	//NumDist->Add(temp_bnumHist); //adding the histogram to my stack
   	//auto temp_cnumHist = (TH1D*)&cnumHist.GetValue();
	//temp_cnumHist->SetFillColor(kBlue); //adding another one and so on....
   	//NumDist->Add(temp_cnumHist);
   	//h3->SetFillColor(kGreen);
   	//hs->Add(h3);
   	//auto QDist = new TCanvas("Quark Dist.","Quark Dist",20,20,400,500); //making a canvas for the stack
   	//auto T = new TText(.5,.5,"TStack");
	//T->SetTextFont(42); T->SetTextAlign(21);
   	//QDist->Divide(2,2);
   	//QDist->cd(2); NumDist->Draw("Quark Dist."); T->DrawTextNDC(.5,.95,"Quark Dist");//this is the stack canvas style used in phd
	//QDist->SaveAs("QDist.root");// save it


	//auto hist1_ptr = (TH1D*)&hist1.GetValue();
	//hs->Add(hist1_ptr)



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

