#ifndef calchisto_hpp
#define calchisto_hpp

enum  channel  {elnu,munu};
const channel
channelAll[] = {elnu,munu};

enum  dataSource  {tzq,ww,wz,zz,ttz,met,cms};
const dataSource
dataSourceAll[] = {tzq,ww,wz,zz,ttz,met,cms};

void calchisto(const channel, const dataSource); 

#endif /* calchisto_hpp */
