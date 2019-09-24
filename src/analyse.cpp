#include "analyse.hpp"
#include <algorithm>
#include "TLorentzVector.h"
#include "sf.hpp"
#include <TCanvas.h>
#include <TText.h>
#include <THStack.h>
#include <TTreeReaderArray.h>

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


constexpr double MIN_ELE_PT{15};
constexpr float MIN_ELE_LEADING_PT{35.f};
constexpr double MAX_ELE_ETA{2.5};
constexpr double ENDCAP_MIN_ETA{1.566};
constexpr double BARREL_MAX_ETA{1.4442};

constexpr double MIN_MU_PT{20};
constexpr float MIN_MU_LEADING_PT{26.f};
constexpr double MAX_MU_ETA{2.4};
constexpr float MU_LOOSE_ISO{0.15f};
constexpr float MU_TIGHT_ISO{0.25f};

constexpr float Z_MASS{91.1876f};
constexpr float Z_MASS_CUT{20.f};

constexpr float MAX_JET_ETA{4.7f};
constexpr float MIN_JET_PT{30.f};
constexpr float JET_ISO{0.4f};
constexpr unsigned MIN_JETS{4};
constexpr unsigned MAX_JETS{6};

constexpr float MAX_BJET_ETA{2.4f};
constexpr float MIN_BTAG_DISC{0.8838f};
constexpr unsigned MIN_BJETS{1};
constexpr unsigned MAX_BJETS{3};

constexpr float W_MASS{80.385f};
constexpr float W_MASS_CUT{20.f};


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
	std:: cout<<"pt size "<<pts.size()<<std::endl;
    }
	//std::cout<<"z mass being passed is "<< boost::numeric_cast<float>(vec.M())<<std::endl;
    	//std::cout<<"size of vector vec "<< vec.size()<<std::endl;
	//if(boost::numeric_cast<float>(vec.M()) > 0 && boost::numeric_cast<float>(vec.M()) <= 100)
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

   	ROOT::RDataFrame d{"Events", "/data/nanoAOD_2017/tZqlvqq/*.root"};
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
	auto enu_selection_function{[](const unsigned int& nElectron_Selection, 
				     const ints& Electron_IDCut_Selection, 
				     const floats& Electron_eta_Selection, 
				     const ints& Electron_charge_Selection, 
				     const floats& Electron_pt_Selection, 
				     const bools& Electron_isPFcand_Selection, 
				     const float& MET_pt_Selection, 
				     const floats& MET_phi_Selection, 
				     const floats& MET_sumEt_Selection) -> bool 
		{

		return 

			nElectron_Selection == 1 && 
			Electron_IDCut_Selection.at(0) >= 4 && 
			(abs(Electron_eta_Selection.at(0)) < 2.5 && 
			(abs(Electron_eta_Selection.at(0)) < 1.4442 ||
			abs(Electron_eta_Selection.at(0)) > 1.566)) &&
			Electron_charge_Selection.at(0) && 
			Electron_pt_Selection.at(0) > 38 && 
			Electron_isPFcand_Selection.at(0) == 1 &&
			MET_pt_Selection >= 0 &&
                        MET_phi_Selection.at(0) &&
                        MET_sumEt_Selection.at(0) >= 0 ;
		}
				};

	auto munu_selection_function{[](const unsigned int& nMuon_Selection, 
				      const bools& Muon_tightId_Selection, 
                                      const chars& Muon_pfIsoId_Selection, 
                                      const floats& Muon_eta_Selection, 
                                      const ints& Muon_charge_Selection, 
                                      const floats& Muon_pt_Selection, 
                                      const bools& Muon_isPFcand_Selection,
                                      const float& MET_pt_Selection, 
                                      const floats& MET_phi_Selection, 
                                      const floats& MET_sumEt_Selection) -> bool 
		{

		return 

			nMuon_Selection == 1 && 
			Muon_tightId_Selection.at(0)== 1 && 
			Muon_pfIsoId_Selection.at(0) >= 4 && 
			abs(Muon_eta_Selection.at(0)) < 2.4 &&
			Muon_pt_Selection.at(0) > 29 &&
			Muon_isPFcand_Selection.at(0) == 1 &&
                        MET_pt_Selection >= 0 &&
                      	MET_phi_Selection.at(0) &&
                        MET_sumEt_Selection.at(0) >= 0 ;

		}
				};


	auto enu_or_munu_selection_function{[](const unsigned int& nElectron_Selection, 
				   	   const ints& Electron_IDCut_Selection, 
				           const floats& Electron_eta_Selection, 
				           const ints& Electron_charge_Selection, 
				           const floats& Electron_pt_Selection, 
				           const unsigned int& nMuon_Selection, 
			                   const bools& Muon_tightId_Selection, 
				           const chars& Muon_pfIsoId_Selection, 
				           const floats& Muon_eta_Selection, 
				           const ints& Muon_charge_Selection, 
				           const floats& Muon_pt_Selection, 
					   const float& MET_pt_Selection, 
					   const floats& MET_phi_Selection, 
					   const floats& MET_sumEt_Selection, 
				           const bools& Electron_isPFcand_Selection, 
				           const bools& Muon_isPFcand_Selection) -> bool 
						{
							if(MET_pt_Selection>=0)
							{
								return 

									(nElectron_Selection == 1 || 
									 nMuon_Selection == 1) && 
									(Electron_IDCut_Selection.at(0) >= 4 || 
									(Muon_tightId_Selection.at(0) == 1 && 
									Muon_pfIsoId_Selection.at(0) >= 4)) && 
									((abs(Electron_eta_Selection.at(0)) < 2.5 &&
									(abs(Electron_eta_Selection.at(0)) < 1.4442 ||
									abs(Electron_eta_Selection.at(0)) > 1.566)) ||
									abs(Muon_eta_Selection.at(0)) < 2.4) &&  
									(Electron_pt_Selection.at(0) > 25 || 
									Muon_pt_Selection > 25) &&
									(Electron_isPFcand_Selection.at(0) == 1 ||
									Muon_isPFcand_Selection.at(0) == 1) &&
									MET_pt_Selection &&
									MET_phi_Selection.at(0) &&
									MET_sumEt_Selection.at(0) >= 0 ;
							}
						}
					};
	auto deltaphi_e_function{[](const unsigned int& nJet,
				    const floats& Electron_phi_Selection,
				    const floats& Jet_phi_Selection)->bool
					{

						for(int i = 0; i < nJet; i++)
						{

							const floats delta_phi_e = Electron_phi_Selection - Jet_phi_Selection.at(i);
							if(-M_PI <= delta_phi_e.at(i) <= M_PI)
							{
								return delta_phi_e.at(i);

							}

						}

					}
				};

	auto deltaphi_mu_function{[](const unsigned int& nJet,
				     const floats& Muon_phi_Selection,
                                     const floats& Jet_phi_Selection)->bool
					{
						for(int i = 0; i < nJet; i++)
						{

							const floats delta_phi_mu = Muon_phi_Selection - Jet_phi_Selection.at(i);
							if(-M_PI <= delta_phi_mu.at(i) <= M_PI)
							{
								return delta_phi_mu.at(i);
							}

						}

					}
				};



	auto deltaRcheck_function{[&deltaphi_e_function, &deltaphi_mu_function](const unsigned int& nJet,
										const floats& Electron_phi_Selection,
										const floats& Muon_phi_Selection,
										const floats& Electron_eta_Selection,
										const floats& Muon_eta_Selection,
										const floats& Jet_eta_Selection,
										const floats& Jet_phi_Selection)->bool
					{
						for(int i = 0; i < nJet; i++)
						{
							const floats delta_eta_e = Electron_eta_Selection - Jet_eta_Selection.at(i);
							const floats delta_eta_mu = Muon_eta_Selection - Jet_eta_Selection.at(i);
							bool delta_phi_e = deltaphi_e_function(nJet, Electron_phi_Selection, Jet_phi_Selection);
							bool delta_phi_mu = deltaphi_e_function(nJet, Muon_phi_Selection, Jet_phi_Selection);
							floats deltaR_ejet = sqrt(pow(delta_phi_e, 2) + pow(delta_eta_e, 2));
							floats deltaR_mujet = sqrt(pow(delta_phi_mu, 2) + pow(delta_eta_mu, 2));

							cout << "after delta r calculation" << endl;

							bool deltaRcheck_e = any_of(deltaR_ejet.begin(), deltaR_ejet.end(), [](int j){return j > 0.4;});
							cout << "after deltaRcheck_e" << endl;
							bool deltaRcheck_mu = any_of(deltaR_mujet.begin(), deltaR_mujet.end(), [](int k){return k > 0.4;});
							cout << "after deltaRcheck_mu" << endl;
							bool deltaRcheck;

							if(deltaRcheck_e == 1 || deltaRcheck_mu == 1)
							{

        							return

        								Electron_eta_Selection.at(i) &&
									Muon_eta_Selection.at(i) &&
									Jet_eta_Selection.at(i) &&
									Electron_phi_Selection.at(i) &&
									Muon_phi_Selection.at(i) &&
									Jet_phi_Selection.at(i);


							}


						}

					}
				};



	auto jet_selection_function{[&deltaRcheck_function](const floats& Jet_pt_Selection, 
							    const floats& Jet_eta_Selection, 
							    const floats& Jet_phi_Selection, 
							    const unsigned int& nJet, 
							    const floats& Electron_pt_Selection, 
							    const floats& Muon_pt_Selection, 
							    const floats& Electron_phi_Selection, 
							    const floats& Muon_phi_Selection, 
							    const floats& Electron_eta_Selection, 
							    const floats& Muon_eta_Selection, 
							    const ints& Jet_jetId_Selection) ->bool 
					{

						for(int i = 0; i < nJet; i++)
						{
							return
								nJet == 4 &&
								Jet_pt_Selection.at(i) > 30 &&
								Jet_eta_Selection.at(i) < 4.7 &&
								Jet_jetId_Selection.at(i) >= 2; //&&
						}
					}
				};

	auto bjet_id{[](
			const ints& tight_jets,
			const floats& btags,
			const floats& etas) 
			{
				return  tight_jets && (btags > 0.8838f) && (etas < 2.4f);
			}
		    };
	auto bjet_cut{[](const ints& bjets) {const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        					return nbjet >= 1 && nbjet <= 3;
					    }
		     };
	



	auto find_lead_mask{[](const ints& mask, const floats& vals) 
				{
					const auto masked_vals{mask * vals};
					const auto max_idx{boost::numeric_cast<size_t>(std::distance(masked_vals.begin(), max_element(masked_vals.begin(), masked_vals.end())))};
					ints lead_mask(masked_vals.size(), 0); // must be ()
					lead_mask.at(max_idx) = 1;
					return lead_mask;
				}
			   };
	auto find_z_pair{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) 
				{//re-write my own code for reconstructing Z mass
					double z_reco_mass{std::numeric_limits<double>::infinity()};
					size_t jet_index_1{std::numeric_limits<size_t>::max()};
					size_t jet_index_2{std::numeric_limits<size_t>::max()};
					const size_t njets{pts.size()};
					for (size_t i{0}; i < njets; ++i)
					{
						for (size_t j{i + 1}; j < njets; ++j)
            					{
                					if (abs(phis.at(i) - phis.at(j))>=2.0)
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
				}
			};

			auto z_mass_cut{[](const float& z_mass) 
						{
        						return abs(z_mass - Z_MASS) < Z_MASS_CUT;
						}
					};
	vector<string> enu_strings = {"nElectron_Selection", 
				     "Electron_IDCut_Selection", 
                                     "Electron_eta_Selection", 
                                     "Electron_charge_Selection", 
                                     "Electron_pt_Selection", 
                                     "Electron_isPFcand_Selection",
				     "MET_pt_Selection", 
                                     "MET_phi_Selection",
				     "MET_sumEt_Selection"};
	vector<string> munu_strings = {"nMuon_Selection", 
				       "Muon_tightId_Selection", 
				       "Muon_pfIsoId_Selection", 
                                       "Muon_eta_Selection", 
                                       "Muon_charge_Selection", 
                                       "Muon_pt_Selection", 
                                       "Muon_isPFcand_Selection",
				       "MET_pt_Selection", 
                                       "MET_phi_Selection",
                                       "MET_sumEt_Selection"};
	vector<string> enumunu_strings = {"nElectron", 
				      "Electron_IDCut_Selection", 
				      "Electron_eta_Selection", 
                                      "Electron_charge_Selection", 
				      "Electron_pt_Selection", 
				      "nMuon", 
				      "Muon_tightId_Selection", 
                                      "Muon_pfIsoId_Selection", 
                                      "Muon_eta_Selection", 
                                      "Muon_charge_Selection", 
                                      "Muon_pt_Selection",
				      "MET_pt_Selection", 
                                      "MET_phi_Selection",
                                      "MET_sumEt_Selection", 
                                      "Electron_isPFcand_Selection", 
                                      "Muon_isPFcand_Selection"};
	vector<string> jet_strings = {"Jet_pt_Selection",
				      "Jet_eta_Selection", 
				      "Jet_phi_Selection",
				      "nJet",
				      "Electron_pt_Selection",
				      "Muon_pt_Selection",
                                      "Electron_phi_Selection", 
				      "Muon_phi_Selection",
				      "Electron_eta_Selection", 
				      "Muon_eta_Selection",
				      "Jet_jetId_Selection"};

	vector<string> z_pair_strings = {"Jet_pt_Selection",
					 "Jet_eta_Selection",
                                         "Jet_phi_Selection",
                                         "Jet_mass",
                                         "Jet_jetId_Selection",
                                         "lead_bjet"};

	vector<string> b_mass_strings = {"Jet_jetId_Selection", "Jet_CSVv2_Selection", "Jet_eta", "Jet_mass", "nJet"};
	vector<string> b_eta_strings = {"Jet_jetId_Selection", "Jet_CSVv2_Selection", "Jet_eta", "nJet"};
	vector<string> b_pt_strings = {"Jet_jetId_Selection", "Jet_CSVv2_Selection", "Jet_eta", "Jet_pt", "nJet"};
	vector<string> b_phi_strings = {"Jet_jetId_Selection", "Jet_CSVv2_Selection", "Jet_eta", "Jet_phi", "nJet"};



	auto d_enu_event_selection = d.Define("nElectron_Selection",{"nElectron"})
					.Define("Electron_pt_Selection", {"Electron_pt"})
                                	.Define("Electron_eta_Selection", {"Electron_eta"})
                                	.Define("Electron_charge_Selection", {"Electron_charge"})
                                	.Define("Electron_IDCut_Selection", {"Electron_cutBased"})
					.Define("Electron_isPFcand_Selection", {"Electron_isPFcand"})
					.Define("MET_pt_Selection",{"MET_pt"})
					.Define("MET_phi_Selection",{"MET_phi"})
					.Define("MET_sumEt_Selection",{"MET_sumEt"})
                                	.Filter(enu_selection_function, enu_strings);

	auto d_enu_jets_selection = d_enu_event_selection.Define("Jet_pt_Selection", {"Jet_pt"})
					 		.Define("Jet_eta_Selection", {"Jet_eta"})
					 		.Define("Jet_phi_Selection", {"Jet_phi"})
					 		.Define("Muon_pt_Selection", {"Muon_pt"})
                                         		.Define("Muon_eta_Selection", {"Muon_eta"})
					 		.Define("Electron_phi_Selection", {"Electron_phi"})
					 		.Define("Muon_phi_Selection", {"Muon_phi"})
					 		.Define("Jet_jetId_Selection", {"Jet_jetId"})
					 		.Define("deltaRcheck", deltaRcheck_function, {"nJet", "Electron_phi", "Electron_eta", "Muon_phi", "Muon_eta", "Jet_phi", "Jet_eta"})
				         		.Filter(jet_selection_function, jet_strings);

	auto d_enu_jets_bjets_selection = d_enu_jets_selection.Define("bjets", bjet_id, {"Jet_jetId", "Jet_btagCSVV2", "Jet_eta"})
                    						.Filter(bjet_cut, {"bjets"}, "b jet cut");

	auto d_enu_z_rec_selection = d_enu_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                 						.Define("z_reco_jets", find_z_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "Jet_jetId", "lead_bjet"})
                 						.Define("z_pair_pt", select<floats>, {"Jet_pt", "z_reco_jets"})
                 						.Define("z_pair_eta", select<floats>, {"Jet_eta", "z_reco_jets"})
                 						.Define("z_pair_phi", select<floats>, {"Jet_phi", "z_reco_jets"})
                 						.Define("z_pair_mass", select<floats>, {"Jet_mass", "z_reco_jets"})
                 						.Define("z_mass", inv_mass, {"z_pair_pt", "z_pair_eta", "z_pair_phi", "z_pair_mass"})
                 						.Filter(z_mass_cut, {"z_mass"}, "z mass cut")
								.Histo1D({"ZReconMassHist_enu_Channel","Recon. Z mass of jets",20,0,200},"z_mass");

 	auto ReconZ_massCanvas = new TCanvas("Reconstrcuted mass from jets", "Reconstructed mass from jets",10,10,700,700);
        d_enu_z_rec_selection->SetLineColor(kBlue);
        d_enu_z_rec_selection->DrawClone();
        d_enu_z_rec_selection->SaveAs("rec_z_mass.root");
        //ReconZ_massCanvas->SaveAs("rec_z_mass.pdf");

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

