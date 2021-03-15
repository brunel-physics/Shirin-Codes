#ifndef calchisto_hpp
#define calchisto_hpp

enum      channel      {elnu,munu};
constexpr channel
          channelAll[]={elnu,munu};

enum      dataSource      {tzq,zz,tz1,tz2,ww,wz,met,wjt,st,stb,stw,stbw,wzll,wjqq,cms,ttb,ttl,ttj};//,wjt,met,st,stb,stw,stbw,ttl,ttj,ttb,cms};
constexpr dataSource
          dataSourceAll[]={tzq,zz,tz1,tz2,ww,wz,met,wjt,st,stb,stw,stbw,wzll,wjqq,cms,ttb,ttl,ttj};//,wjt,met,st,stb,stw,stbw,ttl,ttj,ttb,cms};


enum      PtEtaPhiM      {pt,eta,phi,m};
constexpr PtEtaPhiM
          PtEtaPhiMall[]={pt,eta,phi,m};

void calchisto(const channel, const dataSource);

#endif /* calchisto_hpp */
