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
double pi = 3.14;

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

   	ROOT::RDataFrame d{"Events", "/data/nanoAOD_2017/tZqlvqq/tZqlvqq/*.root"};
	//ROOT::RDataFrame dc{"Events", "/data/nanoAOD_2017/tZqlvqq/388C7D88-9042-E811-9999-FA163E1CC0EA.root"};
	//auto d = dc.Range(0, 100);
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
        	const bool ele_cut{/*tight_ele_pts.size() == N_E &&*/ tight_ele_pts.size() == loose_ele_pts.size()};
        	//bool lead_pt_cut{false};
        	//lead_pt_cut = tight_ele_pts.empty() ? false : *std::max_element(tight_ele_pts.begin(), tight_ele_pts.end()) > MIN_ELE_PT;
		return ele_cut;
        	//return lead_pt_cut && ele_cut;
    	}};

	auto mu_cut{[](const floats& tight_mu_pts, const floats& loose_mu_pts) {
        	const bool mu_cut{/*tight_ele_pts.size() == N_E &&*/ tight_mu_pts.size() == loose_mu_pts.size()};
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
  		float w_reco_mass{std::numeric_limits<float>::infinity()};
                //size_t l_index_1{std::numeric_limits<size_t>::max()};

		for(int i{0}; i< lep_pt.size();i++)
		{
			const float  reco_mass = sqrt( 2 * lep_pt.at(i) * met_pt * (1 - cos(delta_phi(lep_phi.at(i), met_phi))) );
			if (std::abs(W_MASS - reco_mass) < std::abs(W_MASS - w_reco_mass))
                        {
                        	w_reco_mass = reco_mass;
                                //l_index_1 = i;
                        }
		}
			return w_reco_mass;
	}};


	auto w_mass_cut{[](const float& w_mass) {
        	return std::abs(w_mass - W_MASS) < W_MASS_CUT;
   	}};


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

	/*auto jet_delta_phi_cut{[](const floats& phi1, const floats& phi2){
		return delta_phis >=2;//need to work on this... to say look at each element and pick the ones  >= 2
	}};*/

	auto find_z_pair{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet){
		double z_reco_mass{std::numeric_limits<double>::infinity()};
		size_t jet_index_1{std::numeric_limits<size_t>::max()};
		size_t jet_index_2{std::numeric_limits<size_t>::max()};
		const size_t njets{pts.size()};
		cout<<"in Z Loop"<<endl;
		for (size_t i{0}; i < njets; ++i)
		{
			cout<<"first loop "<<"size "<<njets<<endl;
			for (size_t j{i + 1}; j < njets; ++j)
 			{
				cout<< " i & j before removal of b jet "<<i<<" & "<<j<<endl; 
				if (tight_jets[i] != 0 && tight_jets[j] != 0 && lead_bjet[i] != 1 && lead_bjet[j] != 1)
                		{
                    			continue;
	               		}
                                cout<< " i & j after removal of b jet "<<i<<" & "<<j<<endl;
				cout<<"in second loop and passed the jet phi condition"<<endl;
				auto jet1{TLorentzVector{}};
                		auto jet2{TLorentzVector{}};
                		jet1.SetPtEtaPhiM(pts.at(i), etas.at(i), phis.at(i), ms.at(i));
                		jet2.SetPtEtaPhiM(pts.at(j), etas.at(j), phis.at(j), ms.at(j));

               			if (const double reco_mass{(jet1 + jet2).M()}; std::abs(Z_MASS - reco_mass) < std::abs(Z_MASS - z_reco_mass))
                		{
					cout<<"before reco mass"<<endl;
                    			z_reco_mass = reco_mass;
                    			jet_index_1 = i;
                    			jet_index_2 = j;
					cout<<"after reco mass "<<"i "<<i<<" j "<<j<<endl;
					continue;
       				}
        		}
			cout<<"left loop 2"<<endl;
        	}
		cout<<"left loop one"<<endl;
		ints z_pair(njets, 0);
		cout<<"z pair defined "<<z_pair<<endl;
        	z_pair.at(jet_index_1) = 1; cout<<"z pair line 1 "<<z_pair.at(jet_index_1) <<" z pair "<<z_pair<<endl;
        	z_pair.at(jet_index_2) = 1; cout<<"z pair line 2 "<<z_pair.at(jet_index_2) <<" z pair "<<z_pair<<endl;
		cout<<"z pair assigned"<<endl;
   		return z_pair;
	}};

	auto z_mass_cut{[](const float& z_mass){
   		return abs(z_mass - Z_MASS) < Z_MASS_CUT;
	}};


	auto deltaR_z_l{[](const floats& deltaRzl){
		bool deltaRzl_cut{false};
                deltaRzl_cut = deltaRzl.empty() ? false : *std::max_element(deltaRzl.begin(), deltaRzl.end()) > DELTA_R_ZL;
	return deltaRzl_cut;
        }};

        auto ZW_delta_phi_cut{[](const floats& phis1, const floats& phis2){
                float phi1;
		float phi2;
		float deltaphi;
		for(int i; i < phis1.size(); i++)
		{
			phi1 = phis1.at(i);
			for(int j; j < phis2.size(); j++)
			{
				float deltaphi;
				phi2 = phis2.at(j);
				deltaphi = abs(delta_phi(phi1, phi2));
			}
		}
		return deltaphi > DELTA_PHI_ZW;
		//need to work on this... to say look at each elements and pick the ones with delta  >= 2
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
		const floats& w_pair_pt, const floats& w_pair_eta, const floats& w_pair_phi, const float& w_mass){

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
                               	RecoW.SetPtEtaPhiM(w_pair_pt.at(j), w_pair_eta.at(j), w_pair_phi.at(j), w_mass);
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


///////////////////////////////////////////////////////////////////////////Electron Channel/////////////////////////////////////////////////////////////////////////
/*	auto d_enu_event_selection = d.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
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
							.Define("w_e_mass", transvers_W_mass, {"w_e_pt", "w_e_phi", "MET_electron_pt_Selection", "MET_phi_Selection"})
                 					.Filter(w_mass_cut, {"w_e_mass"}, "W mass cut");


	auto d_enu_jets_selection = d_enu_w_selection.Define("jet_e_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_e_eta", "w_e_phi"})
                   					.Define("tight_jets", tight_jet_id, {"jet_e_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
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
							        .Filter(deltaR_z_l,{"z_e_min_dR"}, "delta R ZL")
                 						.Filter(z_mass_cut, {"z_mass"}, "z mass cut");


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
*/
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
                                                        .Define("w_mu_mass", transvers_W_mass, {"w_mu_pt", "w_mu_phi", "MET_mu_pt_Selection", "MET_phi_Selection"})
                                                        .Filter(w_mass_cut, {"w_mu_mass"}, "W mass cut");

        auto d_munu_jets_selection = d_munu_w_selection.Define("jet_mu_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "w_mu_eta", "w_mu_phi"})
                                                        .Define("tight_jets", tight_jet_id, {"jet_mu_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
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
								.Filter(deltaR_z_l,{"z_mu_min_dR"}, "delta R ZL")
								.Filter(ZW_delta_phi_cut,{"z_pair_phi","w_mu_phi"})
                                                                .Filter(z_mass_cut, {"z_mass"}, "z mass cut");
/*
        auto d_munu_brec_selection = d_munu_z_rec_selection.Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                        .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                        .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                        .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                        .Define("nbjets", numberofbjets, {"bjets"})
                                                        .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "nbjets"});

        auto d_munu_top_selection = d_munu_brec_selection .Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", "w_mu_pt", "w_mu_eta", "w_mu_phi", "w_mu_mass"})
                                                        .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                        .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
                                                        .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
                                                        .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"});
*/

/////////////////////////////////////////////////////////////////////// E-Nu Channel histograms AND Canvases /////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////PT//////////////////////////////////////////////////////////////////////////////////
        //auto h_enu_events_ept = d_enu_top_selection.Histo1D({"electron_pt_enu_Channel","electron pt in electron-neutrino channel",20,0,150},"tight_ele_pt");
	auto h_munu_events_mupt = d_munu_z_rec_selection.Histo1D({"mu_pt_enu_Channel","mu pt in mu-neutrino channel",20,0,150},"tight_mu_pt");

        auto h_events_lept_canvas = new TCanvas("lepton pt", "lepton pt",10,10,900,900);
        /*h_enu_events_ept->GetXaxis()->SetTitle("Pt/GeV");
        h_enu_events_ept->GetYaxis()->SetTitle("Events");
        h_enu_events_ept->SetLineColor(kRed);*/
	h_munu_events_mupt->SetLineColor(kBlue);
        //h_enu_events_ept->Draw();
	h_munu_events_mupt->DrawClone("SAME");
        h_events_lept_canvas->BuildLegend();
	h_events_lept_canvas->SaveAs("leptons_pt.root");


        //auto h_enu_events_wmass = d_enu_top_selection.Histo1D({"enu_w_mass","electron-neutrino w mass",20,0,120},"w_e_mass");
	auto h_munu_events_wmass = d_munu_z_rec_selection.Histo1D({"munu_w_mass","mu-neutrino w mass",20,0,120},"w_mu_mass");

        auto h_events_wmass_canvas = new TCanvas("leptons w mass ", "leptons w mass",10,10,900,900);
        /*h_enu_events_wmass->GetXaxis()->SetTitle("mass/GeV/C^2");
        h_enu_events_wmass->GetYaxis()->SetTitle("Events");
        h_enu_events_wmass->SetLineColor(kRed);*/
        h_munu_events_wmass->SetLineColor(kBlue);
        //h_enu_events_wmass->Draw();
	h_munu_events_wmass->DrawClone("SAME");
        h_events_wmass_canvas->BuildLegend();
        h_events_wmass_canvas->SaveAs("leptons_Wmass.root");



        //auto h_enu_events_jpt = d_enu_top_selection.Histo1D({"jet_pt_enu_Channel","jet pt in electron-neutrino channel",20,0,150},"Jet_pt");
	auto h_munu_events_jpt = d_munu_z_rec_selection.Histo1D({"jet_pt_munu_Channel","jet pt in mu-neutrino channel",20,0,150},"Jet_pt");

	auto h_events_jpt_canvas = new TCanvas("jet pt", "jet pt",10,10,900,900);
        /*h_enu_events_jpt->GetXaxis()->SetTitle("Pt/GeV");
        h_enu_events_jpt->GetYaxis()->SetTitle("Events");
        h_enu_events_jpt->SetLineColor(kRed);*/
	h_munu_events_jpt->SetLineColor(kBlue);
	//h_enu_events_jpt->Draw();
	h_munu_events_jpt->DrawClone("SAME");
	h_events_jpt_canvas->BuildLegend();
        h_events_jpt_canvas->SaveAs("Jet_pts.root");

        //auto h_enu_events_zmass = d_enu_z_rec_selection.Histo1D({"Z_mass_enu_Channel","Z mass in electron-neutrino channel",20,0,150},"z_mass");
	auto h_munu_events_zmass = d_munu_z_rec_selection.Histo1D({"Z_mass_munu_Channel","Z mass in mu-neutrino channel",20,0,150},"z_mass");

        auto h_events_zmass_canvas = new TCanvas("Z mass", "Z mass",10,10,900,900);
        /*h_enu_events_zmass->GetXaxis()->SetTitle("mass/GeVC^2");
        h_enu_events_zmass->GetYaxis()->SetTitle("Events");
        h_enu_events_zmass->SetLineColor(kRed);*/
        h_munu_events_zmass->SetLineColor(kBlue);
        //h_enu_events_zmass->Draw();
	h_munu_events_zmass->DrawClone("SAME");
	h_events_zmass_canvas->BuildLegend();
        h_events_zmass_canvas->SaveAs("Z_mass.root");
/*

        auto h_enu_events_topmass = d_enu_top_selection.Histo1D({"enu_top_mass","electron-neutrino top mass",20,0,250},"Top_Mass");
        auto h_munu_events_topmass = d_munu_top_selection.Histo1D({"munu_top_mass","mu-neutrino top mass",20,0,250},"Top_Mass");

        auto h_events_topmass_canvas = new TCanvas("top mass ", "top mass",10,10,900,900);
        h_enu_events_topmass->GetXaxis()->SetTitle("mass/GeV/C^2");
        h_enu_events_topmass->GetYaxis()->SetTitle("Events");
        h_enu_events_topmass->SetLineColor(kRed);
	h_munu_events_topmass->SetLineColor(kBlue);
        h_enu_events_topmass->Draw();
	h_munu_events_topmass->DrawClone("SAME");
        h_events_topmass_canvas->BuildLegend();
        h_events_topmass_canvas->SaveAs("topmass.root");
*/
/*

///////////////////////////////////////////////////////////////////////////ETA //////////////////////////////////////////////////////////////////////////
        auto h_enu_events_eeta = d_enu_jets_selection.Histo1D({"electron_eta_enu_Channel","electron eta in electron-neutrino channel",20,-5,5},"Electron_eta_Selection");

        auto h_enu_events_eeta_canvas = new TCanvas("Electron eta in enu channel", "electron eta in enu",10,10,900,900);
        h_enu_events_eeta->GetXaxis()->SetTitle("eta");
        h_enu_events_eeta->GetYaxis()->SetTitle("Events");
        h_enu_events_eeta->SetLineColor(kRed);
        h_enu_events_eeta_canvas->BuildLegend();
        h_enu_events_eeta->Draw();
        h_enu_events_eeta_canvas->SaveAs("Electron_eta_enu_eSelection.root");

        auto h_enu_events_jeta = d_enu_jets_selection.Histo1D({"jet_eta_enu_Channel","jet eta in electron-neutrino channel",20,-5,5},"Jet_eta_Selection");

        auto h_enu_events_jeta_canvas = new TCanvas("Jet eta in enu channel", "Jet eta in enu",10,10,900,900);
        h_enu_events_jeta->GetXaxis()->SetTitle("eta");
        h_enu_events_jeta->GetYaxis()->SetTitle("Events");
        h_enu_events_jeta->SetLineColor(kRed);
        h_enu_events_jeta_canvas->BuildLegend();
        h_enu_events_jeta->Draw();
        h_enu_events_jeta_canvas->SaveAs("Jet_eta_enu.root");
/*
        auto h_enu_events_zeta = d_enu_z_rec_selection.Histo1D({"Z_eta_enu_Channel","Z eta in electron-neutrino channel",20,-5,5},"z_eta_pair");

        auto h_enu_events_zeta_canvas = new TCanvas("Z eta in enu channel", "Z eta in enu",10,10,900,900);
        h_enu_events_zeta->GetXaxis()->SetTitle("eta");
        h_enu_events_zeta->GetYaxis()->SetTitle("Events");
        h_enu_events_zeta->SetLineColor(kRed);
        h_enu_events_zeta_canvas->BuildLegend();
        h_enu_events_zeta->Draw();
        h_enu_events_zeta_canvas->SaveAs("Z_eta_enu.root");

/////////////////////////////////////////////////////////////////////// PHI////////////////////////////////////////////////////////////////////////////////////////

	auto h_enu_events_ephi = d_enu_jets_selection.Histo1D({"Electron_phi_enu_Channel","Electron phi in electron-neutrino channel",20,0,8},"Electron_phi_Selection");

        auto h_enu_events_ephi_canvas = new TCanvas("Electron phi in enu channel", "Electrion phi in enu",10,10,900,900);
        h_enu_events_ephi->GetXaxis()->SetTitle("phi/rad");
        h_enu_events_ephi->GetYaxis()->SetTitle("Events");
        h_enu_events_ephi->SetLineColor(kRed);
        h_enu_events_ephi_canvas->BuildLegend();
        h_enu_events_ephi->Draw();
        h_enu_events_ephi_canvas->SaveAs("Electron_phi_enu_eSelection.root");

        auto h_enu_events_jphi = d_enu_jets_selection.Histo1D({"jet_phi_enu_Channel","jet phi in electron-neutrino channel",20,0,8},"Jet_phi_Selection");

        auto h_enu_events_jphi_canvas = new TCanvas("jet phi in enu channel", "jet phi in enu",10,10,900,900);
        h_enu_events_jphi->GetXaxis()->SetTitle("phi/rad");
        h_enu_events_jphi->GetYaxis()->SetTitle("Events");
        h_enu_events_jphi->SetLineColor(kRed);
        h_enu_events_jphi_canvas->BuildLegend();
        h_enu_events_jphi->Draw();
        h_enu_events_jphi_canvas->SaveAs("Jet_phi_enu.root");

/*        auto h_enu_events_zphi = d_enu_z_rec_selection.Histo1D({"z_phi_enu_Channel","z phi in electron-neutrino channel",20,0,8},"z_pair_phi");

        auto h_enu_events_zphi_canvas = new TCanvas("z phi in enu channel", "z phi in enu",10,10,900,900);
        h_enu_events_zphi->GetXaxis()->SetTitle("phi/rad");
        h_enu_events_zphi->GetYaxis()->SetTitle("Events");
        h_enu_events_zphi->SetLineColor(kRed);
        h_enu_events_zphi_canvas->BuildLegend();
        h_enu_events_zphi->Draw();
        h_enu_events_zphi_canvas->SaveAs("Z_phi_enu.root");

///////////////////////////////////////////////////////////////////// MU-NEU CHANNEL Histograms and Canvases ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// PT////////////////////////////////////////////////////////////////////////////////////
        auto h_munu_events_mupt = d_munu_jets_selection.Histo1D({"muon_pt_enu_Channel","muon pt in muon-neutrino channel",20,0,150},"Muon_pt_Selection");

        auto h_munu_events_mupt_canvas = new TCanvas("Muon pt in munu channel", "Muon pt in munu",10,10,900,900);
        h_munu_events_mupt->GetXaxis()->SetTitle("Pt/GeV");
        h_munu_events_mupt->GetYaxis()->SetTitle("Events");
        h_munu_events_mupt->SetLineColor(kRed);
        h_munu_events_mupt_canvas->BuildLegend();
        h_munu_events_mupt->Draw();
        h_munu_events_mupt_canvas->SaveAs("Muon_pt_munu_muSelection.root");


        auto h_munu_events_jpt = d_munu_jets_selection.Histo1D({"muon_pt_enu_Channel","muon pt in electron-neutrino channel",20,0,150},"Jet_pt_Selection");

        auto h_munu_events_jpt_canvas = new TCanvas("jet pt in munu channel", "jet pt in munu",10,10,900,900);
        h_munu_events_jpt->GetXaxis()->SetTitle("Pt/GeV");
        h_munu_events_jpt->GetYaxis()->SetTitle("Events");
        h_munu_events_jpt->SetLineColor(kRed);
        h_munu_events_jpt_canvas->BuildLegend();
        h_munu_events_jpt->Draw();
        h_munu_events_jpt_canvas->SaveAs("Jet_pt_munu.root");
/*
        auto h_munu_events_zpt = d_munu_z_rec_selection.Histo1D({"Z_pt_munu_Channel","Z pt in muon-neutrino channel",20,0,150},"z_pair_pt");

        auto h_munu_events_zpt_canvas = new TCanvas("Z pt in munu channel", "Z pt in munu",10,10,900,900);
        h_munu_events_zpt->GetXaxis()->SetTitle("Pt/GeV");
        h_munu_events_zpt->GetYaxis()->SetTitle("Events");
        h_munu_events_zpt->SetLineColor(kRed);
        h_munu_events_zpt_canvas->BuildLegend();
        h_munu_events_zpt->Draw();
        h_munu_events_zpt_canvas->SaveAs("Z_pt_munu.root");


 ////////////////////////////////////////////////////////////////////////////////////ETA////////////////////////////////////////////////////////////////////////

        auto h_munu_events_mueta = d_munu_jets_selection.Histo1D({"muon_eta_munu_Channel","muon eta in Muon-neutrino channel",20,-5,5},"Muon_eta_Selection");

        auto h_munu_events_mueta_canvas = new TCanvas("Muon eta in enu channel", "Muon eta in enu",10,10,900,900);
        h_munu_events_mueta->GetXaxis()->SetTitle("eta");
        h_munu_events_mueta->GetYaxis()->SetTitle("Events");
        h_munu_events_mueta->SetLineColor(kRed);
        h_munu_events_mueta_canvas->BuildLegend();
        h_munu_events_mueta->Draw();
        h_munu_events_mueta_canvas->SaveAs("Muon_eta_munu_muSelection.root");

        auto h_munu_events_jeta = d_munu_jets_selection.Histo1D({"jet_eta_munu_Channel","jet muta in muon-neutrino channel",20,-5,5},"Jet_eta_Selection");

        auto h_munu_events_jeta_canvas = new TCanvas("Jet eta in munu channel", "Jet eta in munu",10,10,900,900);
        h_munu_events_jeta->GetXaxis()->SetTitle("eta");
        h_munu_events_jeta->GetYaxis()->SetTitle("Events");
        h_munu_events_jeta->SetLineColor(kRed);
        h_munu_events_jeta_canvas->BuildLegend();
        h_munu_events_jeta->Draw();
        h_munu_events_jeta_canvas->SaveAs("Jet_eta_munu.root");
/*
        auto h_munu_events_zeta = d_munu_z_rec_selection.Histo1D({"Z_eta_munu_Channel","Z eta in muon-neutrino channel",20,-5,5},"z_eta_pair");

        auto h_munu_events_zeta_canvas = new TCanvas("Z eta in munu channel", "Z eta in munu",10,10,900,900);
        h_munu_events_zeta->GetXaxis()->SetTitle("eta");
        h_munu_events_zeta->GetYaxis()->SetTitle("Events");
        h_munu_events_zeta->SetLineColor(kRed);
        h_munu_events_zeta_canvas->BuildLegend();
        h_munu_events_zeta->Draw();
        h_munu_events_zeta_canvas->SaveAs("Z_eta_munu.root");

////////////////////////////////////////////////////////////////////////////PHI////////////////////////////////////////////////////////////////////////////////
        auto h_munu_events_muphi = d_munu_jets_selection.Histo1D({"Muon_phi_enu_Channel","Muon phi in Muon-neutrino channel",20,0,8},"Muon_phi_Selection");

        auto h_munu_events_muphi_canvas = new TCanvas("Muon phi in munu channel", "Muon phi in enu",10,10,900,900);
        h_munu_events_muphi->GetXaxis()->SetTitle("phi/rad");
        h_munu_events_muphi->GetYaxis()->SetTitle("Events");
        h_munu_events_muphi->SetLineColor(kRed);
        h_munu_events_muphi_canvas->BuildLegend();
        h_munu_events_muphi->Draw();
        h_munu_events_muphi_canvas->SaveAs("Muon_phi_munu_muSelection.root");

        auto h_munu_events_jphi = d_munu_jets_selection.Histo1D({"jet_phi_munu_Channel","jet phi in muon-neutrino channel",20,0,8},"Jet_phi_Selection");

        auto h_munu_events_jphi_canvas = new TCanvas("jet phi in munu channel", "jet phi in munu",10,10,900,900);
        h_munu_events_jphi->GetXaxis()->SetTitle("phi/rad");
        h_munu_events_jphi->GetYaxis()->SetTitle("Events");
        h_munu_events_jphi->SetLineColor(kRed);
        h_munu_events_jphi_canvas->BuildLegend();
        h_munu_events_jphi->Draw();
        h_munu_events_jphi_canvas->SaveAs("Jet_phi_munu.root");
/*

        auto h_munu_events_zphi = d_munu_z_rec_selection.Histo1D({"z_phi_munu_Channel","z phi in Muon-neutrino channel",20,0,8},"z_pair_phi");

        auto h_munu_events_zphi_canvas = new TCanvas("z phi in munu channel", "z phi in munu",10,10,900,900);
        h_munu_events_zphi->GetXaxis()->SetTitle("phi/rad");
        h_munu_events_zphi->GetYaxis()->SetTitle("Events");
        h_munu_events_zphi->SetLineColor(kRed);
        h_munu_events_zphi_canvas->BuildLegend();
        h_munu_events_zphi->Draw();
        h_munu_events_zphi_canvas->SaveAs("Z_phi_munu.root");

*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* 	auto ReconZ_massCanvas = new TCanvas("Reconstrcuted mass from jets", "Reconstructed mass from jets",10,10,900,900);
        d_munu_z_rec_selection->GetXaxis()->SetTitle("GeV/c^2");
        d_munu_z_rec_selection->GetYaxis()->SetTitle("Events");
        d_munu_z_rec_selection->SetLineColor(kRed);
        d_munu_z_rec_selection->Draw();
        d_enu_z_rec_selection->SetLineColor(kBlue);
	d_enu_z_rec_selection->GetXaxis()->SetTitle("GeV/c^2");
	d_enu_z_rec_selection->GetYaxis()->SetTitle("Events");
        d_enu_z_rec_selection->DrawClone("SAME");
	//ReconZ_massCanvas->SetTitle("Reconstructed Z mass from jets in tZq DataSet");
	gStyle->SetOptTitle(0);
   	//TPaveLabel *title = new TPaveLabel(.11,.95,.35,.99,"new title","Reconstructed Z mass from jets in tZq DataSet");
   	//title->Draw(); 
	//auto legend = new TLegend(0.1,0.7,0.48,0.9);
        //legend->SetHeader("Reconstructed Z mass for e-nu and mu-nu Channel","C"); // option "C" allows to center the header
        //legend->AddEntry("d_enu_z_rec_selection","Reconstructed Z mass from electron-neutrino channel","l");
        //legend->AddEntry("d_munu_z_rec_selection","Reconstructed Z mass fro muon-neutrino channel","l");
        //legend->Draw();
	ReconZ_massCanvas->BuildLegend();
        ReconZ_massCanvas->SaveAs("rec_z_mass_inc_MET_cut.root");
*/

/*
        auto d_enu_event_pt_selection = d.Define("nElectron_Selection",{"nElectron"})
                                        .Define("Electron_pt_Selection", {"Electron_pt"})
                                        .Define("Electron_eta_Selection", {"Electron_eta"})
                                        .Define("Electron_charge_Selection", {"Electron_charge"})
                                        .Define("Electron_IDCut_Selection", {"Electron_cutBased"})
                                        .Define("Electron_isPFcand_Selection", {"Electron_isPFcand"})
                                        //.Filter(enu_selection_function, enu_strings)
                                        .Define("MET_pt_Selection",{"MET_pt"})
                                        //.Filter([](float pts){return pts > 0;},{"MET_pt_Selection"})
                                        //.Filter([](const floats pts){return pts.empty() ? false : std::all_of(pts.cbegin(),pts.cend(), [](float pt){return pt > 0.;})$
                                        .Define("MET_phi_Selection",{"MET_phi"})
                                        .Define("MET_sumEt_Selection",{"MET_sumEt"})
                                        .Histo1D({"PTHist_enu_Channel","electron Pt Dist.",20,0,20},"Electron_pt_Selection");

	auto Electron_pt_Canvas = new TCanvas("Electron PT Canvas", "Electron PT Canvas",10,10,700,700);
        d_enu_event_pt_selection->SetLineColor(kBlue);
        d_enu_event_pt_selection->DrawClone();
        Electron_pt_Canvas->SaveAs("Electron_pt_uncut_check.root");

*/

		//std::cout<<"Z mass is "<< *z_mass<<endl;

				 /*.Define("jets_selection",jet_selection_function,{"Jet_pt","Jet_eta","Jet_phi","nJet","Jet_jetId"})
				.Define("Single_lep",e_or_mu_selection_function,{"nElectron","Electron_eta","Electron_pt","nMuon","Muon_eta","Muon_pt","CaloMET_pt","CaloMET_phi","CaloMET_sumEt"})
				.Define("b_jets",bjet_id,{"Jet_btagCSVV2", "Jet_eta"})
				.Filter(bjet_cut, {"b_jets"}, "b jet cut");
		std::cout<< "Event selection number done "<<std::endl;
	auto Z_Rec_M = event_selection.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass"}) //"tight_jets", "lead_bjet"})
                 			   .Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 	                   .Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 			   .Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 			   .Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 			   .Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
					   //.Count();
					   .Histo1D({"ZReconDMassHist","Recon. Z mass of jets",20,0,200},"z_mass");
                 			//.Filter(w_mass_cut, {"w_mass"}, "W mass cut")};
 		//std::cout<<"Event selection number for Z recon "<< *Z_Rec_M <<std::endl;
	//auto Z_Rec_M_Hist = Z_Rec_M.Histo1D({"ZReconDMassHist","Recon. Z mass of jets",20,0,200},"z_mass");
	auto ReconZ_massCanvas = new TCanvas("Reconstrcuted mass from jets", "Reconstructed mass from jets",10,10,700,700);
        //Z_REC_M_Hist->SetLineColor(kBlue);
	Z_Rec_M->DrawClone();
        Z_Rec_M->SaveAs("rec_z_mass.root");*/
	//ReconZ_massCanvas->SaveAs("rec_z_mass.pdf");
/////////////////////////////////////////////////////////////////////////////// Pt And ETa Test/////////////////////////////////////////////////////////////////

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

