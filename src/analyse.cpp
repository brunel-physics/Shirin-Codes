#include "analyse.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/conversion/cast.hpp> 

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include "TLorentzVector.h"

namespace
{
constexpr double MIN_ELE_PT{15};
constexpr float MIN_ELE_LEADING_PT{35.f};
constexpr double MAX_ELE_ETA{2.5};
constexpr double ENDCAP_MIN_ETA{1.655};
constexpr double BARREL_MAX_ETA{1.4442};

constexpr double MIN_MU_PT{20};
constexpr float MIN_MU_LEADING_PT{26.f};
constexpr double MAX_MU_ETA{2.4};

constexpr float Z_MASS{91.1876f};
constexpr float Z_MASS_CUT{20.f};

constexpr float MAX_JET_ETA{4.7f};
constexpr float MIN_JET_PT{4.7f};
constexpr float JET_ISO{0.3f};
constexpr unsigned MIN_JETS{4};
constexpr unsigned MAX_JETS{6};

constexpr float MAX_BJET_ETA{2.4f};
constexpr float MIN_BTAG_DISC{0.8838f};
constexpr unsigned MIN_BJETS{1};
constexpr unsigned MAX_BJETS{2};

constexpr float W_MASS{80.385f};
constexpr float W_MASS_CUT{20.f};

enum class channels {ee, mumu};

auto deltaR (float eta1, float phi1, float eta2, float phi2)
{
    auto pi = boost::math::constants::pi<float>();
    float dEta = eta1 - eta2;
    float dPhi = phi1 - phi2;
    while (std::abs(dPhi) > pi)
    {
        dPhi += (dPhi > 0 ? -2 * pi : 2 * pi);
    }
    return std::sqrt(dEta* dEta + dPhi * dPhi);
}
}

void analyse(int argc, char* argv[])
{
    // using doubles = ROOT::VecOps::RVec<double>;
    using floats = ROOT::VecOps::RVec<float>;
    using ints = ROOT::VecOps::RVec<int>;
    using bools = ROOT::VecOps::RVec<bool>;
    using chars = ROOT::VecOps::RVec<UChar_t>;  // aka 1 byte ints

    auto channel = channels::ee;

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d("Events", "/scratch/data/tZqSkimsRun2017/tZq_eta/*.root");

    // Trigger cuts
    auto get_triggers = [channel]() {
        switch (channel)
        {
            case channels::ee:
                return "HLT_Ele35_WPTight_Gsf || HLT_Ele32_WPTight_Gsf_L1DoubleEG || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL";
            case channels::mumu:
                return "HLT_IsoMu27 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8";
            default:
                throw std::runtime_error("Unknown channel");
        }
    };

    auto d_trig = d.Filter(get_triggers(), "trigger cut");

    // MET filters
    auto d_met = d_trig.Filter("!(Flag_HBHENoiseFilter <= 0 || Flag_HBHENoiseIsoFilter <= 0 || Flag_globalTightHalo2016Filter <= 0 || Flag_EcalDeadCellTriggerPrimitiveFilter <= 0 || Flag_goodVertices <= 0 || Flag_BadPFMuonFilter <= 0|| Flag_BadChargedCandidateFilter <= 0|| Flag_ecalBadCalibFilter <= 0 || Flag_eeBadScFilter <= 0)", "met filter");

    // Lepton cuts
    auto get_nlep = [channel]() -> std::pair<unsigned, unsigned> {
        switch (channel)
        {
            case channels::ee:
                return {2, 0};
            case channels::mumu:
                return {0, 2};
            default:
                throw std::runtime_error("Unknown channel");
        }
    };

    const auto [N_E, N_MU] = get_nlep();

    auto is_good_ele = [](int target_id, bools& isPFs, floats& pts, floats& etas, ints& ids) {
      etas = abs(etas);
      return (isPFs && pts > MIN_ELE_PT
              && ((etas < MAX_ELE_ETA && etas > ENDCAP_MIN_ETA) || (etas < BARREL_MAX_ETA))
              && ids >= target_id);
    };
    auto is_good_tight_ele = [&is_good_ele](bools& isPFs, floats& pts, floats& etas, ints& ids) {
        return is_good_ele(4, isPFs, pts, etas, ids);
    };
    auto is_good_loose_ele = [&is_good_ele](bools& isPFs, floats& pts, floats& etas, ints& ids) {
        return is_good_ele(1, isPFs, pts, etas, ids);
    };

    // TODO: No muon cut-based ID? Use "MVA ID" for now...
    auto is_good_mu = [](int target_id, int target_iso_id, bools& isPFs, floats& pts, floats& etas, chars& ids, chars& iso_ids) {
      etas = abs(etas);
      return (isPFs && pts > MIN_MU_PT && etas < MAX_MU_ETA && ids >= target_id && iso_ids >= target_iso_id);
    };
    auto is_good_tight_mu = [is_good_mu](bools& isPFs, floats& pts, floats& etas, chars& ids, chars& iso_ids) {
        return is_good_mu(3, 4, isPFs, pts, etas, ids, iso_ids);
    };
    auto is_good_loose_mu = [is_good_mu](bools& isPFs, floats& pts, floats& etas, chars& ids, chars& iso_ids) {
        return is_good_mu(1, 2, isPFs, pts, etas, ids, iso_ids);
    };

    // TODO: this os check makes the N_LEP check ofr e in ee/mu in mumu redundant, shorter way to do it?
    auto is_os = [channel](ints& e_qs, ints& mu_qs)
    {
        switch (channel)
        {
            case channels::ee:
                return e_qs.size() == 2 ? std::signbit(e_qs.at(0)) != std::signbit(e_qs.at(1)) : false;
            case channels::mumu:
                return mu_qs.size() == 2 ? std::signbit(mu_qs.at(0)) != std::signbit(mu_qs.at(1)) : false;
            default:
                throw std::runtime_error("Unknown channel");
        }
    };

    // Because N_E and N_MU were declared in a structured binding we need to change them from names into local variables. C++ is very bad.
    auto lep_cut = [channel, N_E = N_E, N_MU = N_MU](floats& tight_ele_pts, floats& loose_ele_pts, floats& tight_mu_pts, floats& loose_mu_pts, bool os) {
        bool ele_cut = tight_ele_pts.size() == N_E && tight_ele_pts.size() == loose_ele_pts.size();
        bool mu_cut = tight_mu_pts.size() == N_MU && tight_mu_pts.size() == loose_mu_pts.size();
        bool lead_pt_cut{false};
        switch (channel)
        {
            case channels::ee:
                lead_pt_cut = tight_ele_pts.empty() ? false : *std::max_element(tight_ele_pts.begin(), tight_ele_pts.end()) > MIN_ELE_LEADING_PT;
                break;
            case channels::mumu:
                lead_pt_cut = tight_mu_pts.empty() ? false : *std::max_element(tight_mu_pts.begin(), tight_mu_pts.end()) > MIN_MU_LEADING_PT;
                break;
            default:
                throw std::runtime_error("Unknown channel");
        }
        return os && lead_pt_cut && ele_cut && mu_cut;
    };

    auto select = [](ints& mask, floats &a) {return a[mask];};

    auto d_lep = d_met.Define("tight_eles", is_good_tight_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                      .Define("tight_ele_pt", select, {"tight_eles", "Electron_pt"})
                      .Define("tight_ele_charge", "Electron_charge[tight_eles]")
                      .Define("loose_eles", is_good_loose_ele, {"Electron_isPFcand", "Electron_pt", "Electron_eta", "Electron_cutBased"})
                      .Define("loose_ele_pt", select, {"tight_eles", "Electron_pt"})
                      .Define("tight_mus", is_good_tight_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_mvaId", "Muon_pfIsoId"})
                      .Define("tight_mu_pt", select, {"tight_mus", "Muon_pt"})
                      .Define("tight_mu_charge", "Muon_charge[tight_mus]")
                      .Define("loose_mus", is_good_loose_mu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_mvaId", "Muon_pfIsoId"})
                      .Define("loose_mu_pt", select, {"tight_mus", "Muon_pt"})
                      .Define("os", is_os, {"tight_ele_charge", "tight_mu_charge"})
                      .Filter(lep_cut, {"tight_ele_pt", "loose_ele_pt", "tight_mu_pt", "loose_mu_pt", "os"}, "lepton cut");

    // Z mass cut
    // TODO: this function is a bit of a hack, is there a better way to get the quantites for the correct lepton?
    auto get_z_lep_quantity_selector = [channel](const std::string &s) {
        switch (channel)
        {
            case channels::ee:
                return "Electron_" + s + "[tight_eles]";
            case channels::mumu:
                return "Muon_" + s + "[tight_mus]";
            default:
                throw std::runtime_error("Unknown channel");
        }
    };

    auto pair_mass = [](floats& pts, floats& etas, floats& phis, floats& ms)
    {
        auto vec0 = TLorentzVector{};
        auto vec1 = TLorentzVector{};
        vec0.SetPtEtaPhiM(pts.at(0), etas.at(0), phis.at(0), ms.at(0));
        vec1.SetPtEtaPhiM(pts.at(1), etas.at(1), phis.at(1), ms.at(1));
        return boost::numeric_cast<float>((vec0 + vec1).M());
    };

    auto z_mass_cut = [](float& z_mass) {
        return std::abs(z_mass - Z_MASS) < Z_MASS_CUT;
    };

    auto d_z = d_lep.Define("z_lep_eta", get_z_lep_quantity_selector("eta"))
                    .Define("z_lep_phi", get_z_lep_quantity_selector("phi"))
                    .Define("z_lep_mass", get_z_lep_quantity_selector("mass"))
                    .Define("z_lep_pt", get_z_lep_quantity_selector("pt"))
                    .Define("z_mass", pair_mass, {"z_lep_pt", "z_lep_eta", "z_lep_phi", "z_lep_mass"})
                    .Filter(z_mass_cut, {"z_mass"}, "Z mass cut");

    // Jet cuts
    // TODO: Smearing
    auto tight_jet_id = [](floats& jet_lep_min_dRs, floats& pts, floats& etas, ints &ids) {
        return (pts > MIN_JET_PT) && (etas < MAX_JET_ETA) && (jet_lep_min_dRs > JET_ISO) && (ids >= 2);
    };

    auto jet_lep_min_deltaR = [](floats& jet_etas, floats& jet_phis, floats& lep_etas, floats& lep_phis) {
        floats min_dRs{};
        std::transform(jet_etas.begin(), jet_etas.end(), jet_phis.begin(), std::back_inserter(min_dRs), [&](float jet_eta, float jet_phi) {return std::min(deltaR(jet_eta, jet_phi, lep_etas.at(0), lep_phis.at(0)), deltaR(jet_eta, jet_phi, lep_etas.at(1), lep_phis.at(1)));});
        return min_dRs;
    };

    auto jet_cut = [](ints& tight_jets) {
        auto njet = std::count_if(tight_jets.begin(), tight_jets.end(), [](int i){return i;});
        return njet >= MIN_JETS && njet <= MAX_JETS;
    };

    auto d_jet = d_z.Define("jet_lep_min_dR", jet_lep_min_deltaR, {"Jet_eta", "Jet_phi", "z_lep_eta", "z_lep_phi"})
                    .Define("tight_jets", tight_jet_id, {"jet_lep_min_dR", "Jet_pt", "Jet_eta", "Jet_jetId"})
                    .Filter(jet_cut, {"tight_jets"}, "Jet cut   ");

    // B jet cut
    auto bjet_id = [](ints& tight_jets, floats& btags, floats& etas) {
        return tight_jets && (btags > MIN_BTAG_DISC) && (etas < MAX_BJET_ETA);
    };

    auto bjet_cut = [](ints& bjets) {
        auto nbjet = std::count_if(bjets.begin(), bjets.end(), [](int i){return i;});
        return nbjet >= MIN_BJETS && nbjet <= MAX_BJETS;
    };

    auto d_bjet = d_jet.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", "Jet_eta"})
                       .Filter(bjet_cut, {"bjets"}, "b jet cut  ");

    // W mass cut
    auto find_lead_mask = [](ints& mask, floats &vals) {
        auto masked_vals = mask * vals;
        auto max_idx = boost::numeric_cast<size_t>(std::distance(masked_vals.begin(), max_element(masked_vals.begin(), masked_vals.end())));
        ints lead_mask(masked_vals.size(), 0);
        lead_mask.at(max_idx) = 1;
        return lead_mask;
    };

    auto find_w_pair = [](floats &pts, floats &etas, floats &phis, floats &ms, ints& tight_jets, ints& lead_bjet) {
        double w_reco_mass = std::numeric_limits<double>::infinity();
        size_t jet_index_1 = std::numeric_limits<size_t>::max();
        size_t jet_index_2 = std::numeric_limits<size_t>::max();
        size_t njets = pts.size();

        for (size_t i = 0; i < njets; ++i)
        {
            for (size_t j = i + 1; j < njets; ++j)
            {
                if (tight_jets[i] != 0 && tight_jets[j] != 0
                        && lead_bjet[i] != 1 && lead_bjet[j] != 1)
                {
                    continue;
                }

                auto jet1 = TLorentzVector{};
                auto jet2 = TLorentzVector{};
                jet1.SetPtEtaPhiM(pts.at(i), etas.at(i), phis.at(i), ms.at(i));
                jet2.SetPtEtaPhiM(pts.at(j), etas.at(j), phis.at(j), ms.at(j));

                if (double reco_mass =(jet1 + jet2).M(); std::abs(W_MASS - reco_mass) < std::abs(W_MASS - w_reco_mass))
                {
                    w_reco_mass = reco_mass;
                    jet_index_1 = i;
                    jet_index_2 = j;
                }
            }
        }

        ints w_pair(njets, 0);
        w_pair.at(jet_index_1) = 1;
        w_pair.at(jet_index_2) = 1;
        return w_pair;
    };

    auto w_mass_cut = [](float& w_mass) {
        return std::abs(w_mass - W_MASS) < W_MASS_CUT;
    };

    auto d_w = d_bjet.Define("lead_bjet", find_lead_mask, {"bjets", "Jet_pt"})
                     .Define("w_reco_jets", find_w_pair, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass", "tight_jets", "lead_bjet"})
                     .Define("w_pair_pt", "Jet_pt[w_reco_jets]")
                     .Define("w_pair_eta", "Jet_eta[w_reco_jets]")
                     .Define("w_pair_phi", "Jet_phi[w_reco_jets]")
                     .Define("w_pair_mass", "Jet_mass[w_reco_jets]")
                     .Define("w_mass", pair_mass, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_pair_mass"})
                     .Filter(w_mass_cut, {"w_mass"}, "W mass cut");

    auto allCutsReport = d.Report();
    std::cout << "Name\t\tAll\tPass\tEfficiency" << std::endl;
    for (auto&& cutInfo : allCutsReport) {
      std::cout << cutInfo.GetName() << '\t' << cutInfo.GetAll() << '\t' << cutInfo.GetPass() << '\t' << cutInfo.GetEff() << " %" << std::endl;
    }
}
