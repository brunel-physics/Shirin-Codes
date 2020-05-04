#ifndef calchisto_hpp
#define calchisto_hpp

enum  channels  {elnu,munu};
const channels
channelsAll[] = {elnu,munu};

enum  dataSource  {tzq,ww,wz,zz,ttz};
const dataSource
dataSourceAll[] = {tzq,ww,wz,zz,ttz};
/*
enum  exptData  {CMS,MET};
const exptData
exptDataAll[] = {CMS,MET};
*/
void calchisto(const dataSource); 

#endif
