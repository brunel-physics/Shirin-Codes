#ifndef calchisto_hpp
#define calchisto_hpp

enum      channel      {elnu,munu};
constexpr channel
          channelAll[]={elnu,munu};

enum      dataSource      {zz,ttz,tzq,ww,wz,met,ttb,cms};
constexpr dataSource
          dataSourceAll[]={zz,ttz,tzq,ww,wz,met,ttb,cms};


enum      PtEtaPhiM      {pt,eta,phi,m};
constexpr PtEtaPhiM
          PtEtaPhiMall[]={pt,eta,phi,m};

void calchisto(const channel, const dataSource);

#endif /* calchisto_hpp */
