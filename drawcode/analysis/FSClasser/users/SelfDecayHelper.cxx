#include <EVENT/MCParticle.h>
#include "SelfDecayHelper.h"
#include "lcio.h"
#include <iostream>
using namespace lcio;
using namespace std;
lcio::MCParticle * SelfDecay_RecurseDown(lcio::MCParticle * p){
/*
     if((p->getDaughters().size() == 1) && p->getDaughters().at(0)->getPDG() == 94) // Higgs decay product
     {
         for(unsigned int i =0; i< p->getDaughters().at(0)->getDaughters().size(); i++){
              if(p->getDaughters().at(0)->getDaughters().at(i)->getPDG() == p->getPDG() )
              { return SelfDecay_RecurseDown(p->getDaughters().at(0)->getDaughters().at(i));}    
         } 
     }

     if(p->getDaughters().size() > 1) return p;
     else if( p->getDaughters().size() == 1 &&  p->getDaughters().at(0)->getPDG() == p->getPDG() ){
         MCParticle * p_tmp = p;
         while(p_tmp->getDaughters().size() == 1 && (p_tmp->getDaughters().at(0)->getPDG() == p_tmp->getPDG()  ) ){
               p_tmp = p_tmp->getDaughters().at(0);
         }
         return p_tmp;
     }
     else return 0;
*/
     if( p->getDaughters().size() == 1 &&  p->getDaughters().at(0)->getPDG() == p->getPDG() && (p->getDaughters().at(0)->getParents()).size() == 1)
         return SelfDecay_RecurseDown(p->getDaughters().at(0));
     else return p;
}

lcio::MCParticle * SelfDecay_RecurseUp(lcio::MCParticle * p){
//    cout<<"SelfDecay_RecurseUp: input type: "<<p->getPDG()<<", with energy "<<p->getEnergy()<<endl;
//    cout<<"SelfDecay_RecurseUp: p->getParents().size() = "<<p->getParents().size()<<endl;
    if( p->getParents().size() == 1 &&  p->getParents().at(0)->getPDG() == p->getPDG() &&  p->getParents().at(0)->getDaughters().size() == 1 ){
//         cout<<"SelfDecay_RecurseUp: Enterloop"<<endl;
/*
         MCParticle * p_tmp = p;
         while(p_tmp->getParents().size() == 1 && (p_tmp->getParents().at(0)->getPDG() == p_tmp->getPDG()  ) && ( p->getParents().at(0)->getDaughters().size() == 1)  ){
               p_tmp = p_tmp->getParents().at(0);
         }
         return p_tmp;
*/
         return SelfDecay_RecurseUp(p->getParents().at(0));
     }
     else {/*cout<<"SelfDecay_RecurseUp: stop! will return input pointer"<<endl;*/return p;}
}

bool OnOneBranch(MCParticle* p1, MCParticle* p2){ //check if p2 is from p1, after several radiation
     if(SelfDecay_RecurseUp(p1) == SelfDecay_RecurseUp(p2)) return true; //Same parton
     if(p1->getPDG() != p2->getPDG()) return false;
     bool result =false;
     while(SelfDecay_RecurseUp(p2)->getParents().size()>0 && SelfDecay_RecurseUp(p2)->getParents().at(0)->getPDG()!=92 && SelfDecay_RecurseUp(p2)->getParents().at(0)->getPDG()!=94){
          MCParticle* p2_tmp = SelfDecay_RecurseUp(p2);
          if(p2_tmp->getParents().size() == 1 && p2_tmp->getParents().at(0)->getDaughters().size() ==2 ){
               if(p2_tmp == p2_tmp->getParents().at(0)->getDaughters().at(0)){
                    if(p2_tmp->getParents().at(0)->getDaughters().at(1)->getPDG() == 21 || p2_tmp->getParents().at(0)->getDaughters().at(1)->getPDG() == 22){
                        p2_tmp = p2_tmp->getParents().at(0);return OnOneBranch(p1,p2_tmp);
                    }
               }else if(p2_tmp == p2_tmp->getParents().at(0)->getDaughters().at(1)){
                    if(p2_tmp->getParents().at(0)->getDaughters().at(0)->getPDG() == 21 || p2_tmp->getParents().at(0)->getDaughters().at(0)->getPDG() == 22){
                        p2_tmp = p2_tmp->getParents().at(0);return OnOneBranch(p1,p2_tmp);  
                    }
               }else return false;//Unknown vertex
          }  
     }
     return result;  
}
