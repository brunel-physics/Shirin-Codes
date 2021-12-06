#ifndef calchisto_hpp
#define calchisto_hpp

enum      channel      {elnu,munu};
constexpr channel
          channelAll[]={elnu,munu};

enum      dataSource      {tzq,zz,tz1,tz2,ww,wz,met,wjt,wjx,st,stb,stw,stbw,wzll,wjqq,zjt1,zjt2,zjt3,zjqq,cms,ttb,ttl,ttj};//,wjt,met,st,stb,stw,stbw,ttl,ttj,ttb,cms};
constexpr dataSource
          dataSourceAll[]={tzq,zz,tz1,tz2,ww,wz,met,wjt,wjx,st,stb,stw,stbw,wzll,wjqq,zjt1,zjt2,zjt3,zjqq,cms,ttb,ttl,ttj};//,wjt,met,st,stb,stw,stbw,ttl,ttj,ttb,cms};


enum      PtEtaPhiM      {pt,eta,phi,m};
constexpr PtEtaPhiM
          PtEtaPhiMall[]={pt,eta,phi,m};

// The columns name saved by snap shot  (TTreeMaker.cxx && NPL_run.cxx) and read by ttreehisto file
std::initializer_list< std::string >  mc__columns = {"loose_leps","lep__pt","lep_eta","lep_phi","lep_mas"
                  ,"tJ_btagCSVv2","cjer","cTJer","cmet_dpx","cmet_dpy","btagB","btagP"
                  ,"cmet__pt","cmet_phi","cmet_sEt","fin_jets__pt"
                  ,"fin_jets_eta","fin_jets_phi","fin_jets_mas","genW"
                  ,"is_btag_numer","is_btag_denom","no_btag_numer"
                  ,"no_btag_denom","is_btag_numer__pt","is_btag_denom__pt"
                  ,"no_btag_numer__pt","no_btag_denom__pt","is_btag_numer_eta"
                  ,"is_btag_denom_eta","no_btag_numer_eta","no_btag_denom_eta"
                  ,"sfi","sfj","lep_nu_invmass","fin_jets_Dph","bjet__pt","bjet_eta"
                  ,"bjet_phi","bjet_mas","z_reco_jets","z_pair__pt","z_pair_eta","z_pair_phi"
                  ,"z_pair_mas","z_mas","z__pt","z_eta","z_phi","z_lep_min_dR","zw_Dph","zmet_Dph"
                  ,"z_pair_Dph","WZ_deltaR","recoTtop","ttop__pt","ttop_eta","ttop_phi","ttop_mas"
                  ,"NoEffBTagged","P___ej","P_sfej","P_sf_i","btag_w","ttbSF","lepSF","mostSF"
                  ,"sf","EvtWeight","roccorSF","LHEPdfWeight","MET_sumEt","MET_pt","tw_lep_mas"
                  ,"nw_lep__pt","nw_lep_eta","nw_lep_phi","nw_lep_mas","nw_cmet__pt","nw_cmet_phi"
                  ,"nw_cmet_sEt","nw_fin_jets__pt","nw_fin_jets_eta","nw_fin_jets_phi","nw_fin_jets_mas"
                  ,"nw_bjet__pt","nw_bjet_eta","nw_bjet_phi","nw_bjet_mas","nw_jet_lep_min_dR"
                  ,"nw_z_lep_min_dR","nw_WZ_deltaR","nw_tw_lep_mas","nw_z_mas","nw_fin_jets_Dph"
                  ,"nw_z_pair_Dph","nw_zmet_Dph","nw_zw_Dph","nw_lep_nu_invmass","nw_ttop__pt"
                  ,"nw_ttop_mas"
                  };

std::initializer_list< std::string >  cms_columns = {"cmet__pt","cmet_phi" ,"fin_jets__pt","fin_jets_eta","fin_jets_phi","fin_jets_mas"
                        ,"btag_w","LHEPdfWeight","lepSF","mostSF","lep__pt","lep_eta","lep_phi","lep_mas","tw_lep_mas"
                        ,"lep_nu_invmass","fin_jets_Dph","lead_bjet","bjet__pt","bjet_eta","bjet_phi"
                        ,"bjet_mas","z_reco_jets","z_pair__pt","z_pair_eta","z_pair_phi","z_pair_mas"
                        ,"z_mas","z__pt","z_eta","z_phi","z_lep_min_dR","zw_Dph","zmet_Dph","z_pair_Dph"
                        ,"WZ_deltaR","recoTtop","ttop__pt","ttop_eta","ttop_phi","ttop_mas"
                        ,"MET_sumEt","MET_pt","sf","EvtWeight","btagB","btagP"
                        ,"nw_lep__pt","nw_lep_eta","nw_lep_phi","nw_lep_mas","nw_cmet__pt","nw_cmet_phi"
                        ,"nw_cmet_sEt","nw_fin_jets__pt","nw_fin_jets_eta","nw_fin_jets_phi","nw_fin_jets_mas"
                        ,"nw_bjet__pt","nw_bjet_eta","nw_bjet_phi","nw_bjet_mas","nw_jet_lep_min_dR"
                        ,"nw_z_lep_min_dR","nw_WZ_deltaR","nw_tw_lep_mas","nw_z_mas","nw_fin_jets_Dph"
                        ,"nw_z_pair_Dph","nw_zmet_Dph","nw_zw_Dph","nw_lep_nu_invmass","nw_ttop__pt"
                        ,"nw_ttop_mas"
                        };


void calchisto(const channel, const dataSource);

#endif /* calchisto_hpp */
