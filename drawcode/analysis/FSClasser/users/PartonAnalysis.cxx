#include <EVENT/MCParticle.h>
#include "PartonAnalysis.h"
#include "SelfDecayHelper.h"
#include "lcio.h"
#include "FSClasserProcessor.h"
using namespace lcio;

vector<MCParticle *> FSClasserProcessor::MakePartonList(MCTruthHelper* truthhelper){
        vector<MCParticleImpl*>  mc_list = truthhelper->AllParticleList();

        vector<MCParticle*>  mc_before_hadronization; //Parents of pdg = 92 particle
        mc_before_hadronization.clear();


        for(unsigned int i=0; i<mc_list.size(); i++){
             if(mc_list.at(i)->getPDG() == 92){
                 for(unsigned int j = 0; j< mc_list.at(i)->getParents().size();j++){//check if duplicate ones?
                        mc_before_hadronization.push_back(mc_list.at(i)->getParents().at(j));                          
                 }
             }  //EndOf Loop in 92's parents
        }//Endof Loop of mc particle

//        cout<<"PartonListMaker : Parton number just before hadronization = "<<mc_before_hadronization.size()<<endl;        
        map<MCParticle*, vector<MCParticle*> > parton_map;
        if(mc_before_hadronization.size() > 0)
             parton_map = GetPrimaryParton(mc_before_hadronization);

//        cout<<"PartonListMaker : PrimaryParton number = "<<parton_map.size()<<endl;
        map<MCParticle*, vector<MCParticle*> >::iterator parton_map_itr;
        
        vector<MCParticle*> result;          

        for(parton_map_itr = parton_map.begin();parton_map_itr!=parton_map.end();++parton_map_itr){
               if(parton_map_itr->second.size() ==0){
                   std::cout<<"Empty parton list for parton "<<parton_map_itr->first->getPDG()<<std::endl; 
                   continue;
               }
               if(parton_map_itr->second.size() ==1){
                   result.push_back(SelfDecay_RecurseUp(parton_map_itr->second.at(0))); continue;
               }
               if(abs(parton_map_itr->first->getPDG()) <6){  //Gluon Radiation Finding
//                  cout<<"This is a quark with pdgid = "<<parton_map_itr->first->getPDG()<<endl;
                  vector<MCParticle*> v_tmp = parton_map_itr->second; 
                  vector<MCParticle*> gluon_brothers;
                  vector<MCParticle*> gluon;
                  for(unsigned int i = 0; i< parton_map_itr->second.size();i++){                  
                        bool quarkfromgluon = false;
//                        cout<<"Mark A"<<endl;        
                        if(abs(v_tmp.at(i)->getPDG())<6 && SelfDecay_RecurseUp(v_tmp.at(i)) != SelfDecay_RecurseUp(parton_map_itr->first)) {
                             MCParticle * quark_tmp = SelfDecay_RecurseUp(v_tmp.at(i));
                             while(SelfDecay_RecurseUp(SelfDecay_RecurseUp(quark_tmp)->getParents().at(0)) != SelfDecay_RecurseUp(parton_map_itr->first) ){
                                 if(SelfDecay_RecurseUp(quark_tmp)->getParents().at(0)->getPDG() == 21){
//                                        cout<<"PartonAnalysis: parimary parton with pdg "<<parton_map_itr->first->getPDG()<<", secondary quark with pdg "<<quark_tmp->getPDG()<<endl;
                                        quarkfromgluon=true;break;
                                 }
                                 quark_tmp = SelfDecay_RecurseUp(SelfDecay_RecurseUp(quark_tmp)->getParents().at(0));
                             }
                        }
//                        cout<<"Mark B"<<endl; 
                        if(v_tmp.at(i)->getPDG() == 21 || (abs(v_tmp.at(i)->getPDG())<6 && quarkfromgluon)){
                             MCParticle * gluon_tmp;
//                             cout<<"PartonAnalysis: v_tmp["<< i<<"] is a "<<v_tmp.at(i)->getPDG()<<" with energy "<<v_tmp.at(i)->getEnergy()<<endl;
                             if(v_tmp.at(i)->getPDG() == 21) gluon_tmp = v_tmp.at(i);
//                             cout<<"Mark C"<<endl;
                             if(abs(v_tmp.at(i)->getPDG())<6){
//                                  cout<<"Mark C1"<<endl;
                                  
                                  MCParticle * quark_tmp = SelfDecay_RecurseUp(v_tmp.at(i));
//                                  cout<<"Mark C2"<<endl;
                                  while(SelfDecay_RecurseUp(SelfDecay_RecurseUp(quark_tmp)->getParents().at(0)) != SelfDecay_RecurseUp(parton_map_itr->first) ){
//                                      cout<<"Mark C3"<<endl;
//                                      cout<<"PartonAnalysis: quark "<<quark_tmp->getPDG()<<", with energy "<<quark_tmp->getEnergy()<<endl;
//                                      cout<<"PartonAnalysis: gradefather "<<SelfDecay_RecurseUp(SelfDecay_RecurseUp(SelfDecay_RecurseUp(quark_tmp)->getParents().at(0))->getParents().at(0))<<endl;
//                                      cout<<"PartonAnalysis: parton_map_itr->first : "<<SelfDecay_RecurseUp(parton_map_itr->first)<<endl;
                                      if(SelfDecay_RecurseUp(quark_tmp)->getParents().at(0)->getPDG() == 21 && OnOneBranch(parton_map_itr->first,SelfDecay_RecurseUp(SelfDecay_RecurseUp(quark_tmp)->getParents().at(0))->getParents().at(0)) ){
                                            gluon_tmp = SelfDecay_RecurseUp(quark_tmp)->getParents().at(0);break;
                                      } 
                                      quark_tmp = SelfDecay_RecurseUp(quark_tmp)->getParents().at(0);
//                                      cout<<"Mark C4"<<endl;
                                  }
                             }
//                             cout<<"gluon energy before: "<<SelfDecay_RecurseUp(gluon_tmp)->getEnergy()<<", parents pdg"<<SelfDecay_RecurseUp(gluon_tmp)->getPDG()<<endl;
                             //cout<<"gluon energy before: "<<gluon_tmp->getEnergy()<<", parents pdg"<<gluon_tmp->getParents().at(0)->getPDG()<<endl;
                             //cout<<"PartonAnalysis: gluon energy before: "<<gluon_tmp->getEnergy()<<", parents pdg"<<SelfDecay_RecurseUp(gluon_tmp)->getParents().at(0)->getPDG()<<endl;
//                             cout<<"PartonAnalysis: gluon energy before: "<<gluon_tmp->getEnergy()<<endl;
                             while(SelfDecay_RecurseUp(gluon_tmp)->getParents().size()>0&&SelfDecay_RecurseUp(gluon_tmp)->getParents().at(0)->getPDG() == 21){ //Make sure only initial gluon considered
                                 gluon_tmp = SelfDecay_RecurseUp(gluon_tmp)->getParents().at(0);
//                                 cout<<"gluon energy : "<<gluon_tmp->getEnergy()<<endl;
                             }
//                             cout<<"PartonAnalysis: gluon energy after: "<<gluon_tmp->getEnergy()<<endl;
                             TLorentzVector p_gluon,p_gluon_mother;
                             p_gluon.SetPxPyPzE(gluon_tmp->getMomentum()[0],gluon_tmp->getMomentum()[1],gluon_tmp->getMomentum()[2],gluon_tmp->getEnergy());
                             MCParticle *gluon_mother = SelfDecay_RecurseUp(gluon_tmp)->getParents().at(0); 
                             p_gluon_mother.SetPxPyPzE(gluon_mother->getMomentum()[0],gluon_mother->getMomentum()[1],gluon_mother->getMomentum()[2],gluon_mother->getEnergy()); 
                             TLorentzVector p_gluon_brothers = p_gluon_mother-p_gluon;
//                             cout<<"PartonAnalysis: DeltaR between gluon with energy "<<p_gluon.E()<<" and quark with energy "<<p_gluon_brothers.E()<<" is "<<p_gluon.DeltaR(p_gluon_brothers)<<endl;
                             if(p_gluon.DeltaR(p_gluon_brothers) > 0.7 && p_gluon.Pt()>15.0 && (p_gluon.E()/p_gluon_mother.E())>0.2) {
/*
                                    if(find(result.begin(),result.end(),gluon_tmp) == result.end() ){
                                       result.push_back(gluon_tmp);
                                    }
*/
                                    if(find(gluon.begin(),gluon.end(),SelfDecay_RecurseUp(gluon_tmp)) == gluon.end() ){  
                                       gluon.push_back(SelfDecay_RecurseUp(gluon_tmp));
                                       if(SelfDecay_RecurseDown(gluon_tmp)->getDaughters().size() == 2 && abs(SelfDecay_RecurseDown(gluon_tmp)->getDaughters().at(0)->getPDG())<6){ //Radiated Gluon Splitting
                                            MCParticle* quark1 = SelfDecay_RecurseDown(gluon_tmp)->getDaughters().at(0);
                                            MCParticle* quark2 = SelfDecay_RecurseDown(gluon_tmp)->getDaughters().at(1);
                                            TLorentzVector p_quark1,p_quark2;
                                            p_quark1.SetPxPyPzE(quark1->getMomentum()[0],quark1->getMomentum()[1],quark1->getMomentum()[2],quark1->getEnergy());
                                            p_quark2.SetPxPyPzE(quark2->getMomentum()[0],quark2->getMomentum()[1],quark2->getMomentum()[2],quark1->getEnergy());
                                            cout<<"Delta R between quark1,2"<<p_quark1.DeltaR(p_quark2)<<endl;
                                            if(p_quark1.E()/(p_quark1.E()+p_quark2.E())>0.25 && p_quark2.E()/(p_quark1.E()+p_quark2.E())>0.25 && ( (p_quark1.DeltaR(p_quark2)>0.7 && p_quark1.E()>10.0 && p_quark2.E()>10.0)||(p_quark1.DeltaR(p_quark2)>0.5 && p_quark1.E()>30.0 && p_quark2.E()>30.0)) ){
                                                 if( find(result.begin(),result.end(),SelfDecay_RecurseUp(quark1)) == result.end()) 
                                                     result.push_back(SelfDecay_RecurseUp(quark1));
                                                 if( find(result.begin(),result.end(),SelfDecay_RecurseUp(quark2)) == result.end())
                                                     result.push_back(SelfDecay_RecurseUp(quark2));
                                            }else if(abs(quark1->getPDG())>=3 && (p_quark1.E()/(p_quark1.E()+p_quark2.E())>0.75 || p_quark2.E()/(p_quark1.E()+p_quark2.E())>0.75) )  {
                                                  if(quark1->getEnergy() > quark2->getEnergy()){
                                                      if( find(result.begin(),result.end(),SelfDecay_RecurseUp(quark1)) == result.end())
                                                        result.push_back(SelfDecay_RecurseUp(quark1)); 
                                                  }else{
                                                      if( find(result.begin(),result.end(),SelfDecay_RecurseUp(quark2)) == result.end())
                                                        result.push_back(SelfDecay_RecurseUp(quark2));
                                                  }
                                            }else{
                                                  if(find(result.begin(),result.end(),SelfDecay_RecurseUp(gluon_tmp)) == result.end())
                                                    result.push_back(SelfDecay_RecurseUp(gluon_tmp));       
                                            }         
                                         }
                                         else if(SelfDecay_RecurseDown(gluon_tmp)->getDaughters().size() == 2 && abs(SelfDecay_RecurseDown(gluon_tmp)->getDaughters().at(0)->getPDG())==21){
                                            MCParticle * gluon1 = SelfDecay_RecurseDown(gluon_tmp)->getDaughters().at(0);
                                            MCParticle * gluon2 = SelfDecay_RecurseDown(gluon_tmp)->getDaughters().at(1);
                                            TLorentzVector p_gluon1,p_gluon2;
                                            p_gluon1.SetPxPyPzE(gluon1->getMomentum()[0],gluon1->getMomentum()[1],gluon1->getMomentum()[2],gluon1->getEnergy());
                                            p_gluon2.SetPxPyPzE(gluon2->getMomentum()[0],gluon2->getMomentum()[1],gluon2->getMomentum()[2],gluon2->getEnergy());  
                                            if(p_gluon1.E()/(p_gluon1.E() + p_gluon2.E())>0.3 && p_gluon2.E()/(p_gluon1.E() + p_gluon2.E())>0.3 && ((p_gluon1.E() > 15.0 && p_gluon2.E() > 15.0 && p_gluon1.DeltaR(p_gluon2)>1.0) || (p_gluon1.E() > 30.0 && p_gluon2.E() > 30.0 && p_gluon1.DeltaR(p_gluon2)>0.7))  )  {
                                              if( find(result.begin(),result.end(),SelfDecay_RecurseUp(gluon1)) == result.end())
                                                 result.push_back(SelfDecay_RecurseUp(gluon1));
                                              if( find(result.begin(),result.end(),SelfDecay_RecurseUp(gluon2)) == result.end()) 
                                                 result.push_back(SelfDecay_RecurseUp(gluon2));
                                            }else{
                                                if( find(result.begin(),result.end(),SelfDecay_RecurseUp(gluon_tmp)) == result.end())
                                                    result.push_back(SelfDecay_RecurseUp(gluon_tmp));
                                            }
                                         }  
                                         else{    
                                            if(find(result.begin(),result.end(),SelfDecay_RecurseUp(gluon_tmp)) == result.end()) 
                                              result.push_back(SelfDecay_RecurseUp(gluon_tmp)); 
                                         }
                                    }
                                    for(unsigned int j=0; j < gluon_mother->getDaughters().size();j++){
                                       if(SelfDecay_RecurseUp(gluon_tmp) != gluon_mother->getDaughters().at(j))
                                          if(find(gluon_brothers.begin(),gluon_brothers.end(),gluon_mother->getDaughters().at(j)) == gluon_brothers.end())
                                            gluon_brothers.push_back(gluon_mother->getDaughters().at(j));  
                                    }
                             }
                        }   
                  } //Endof the gluon finding
//                  cout<<"Mark B"<<endl;
                  vector<MCParticle*>::iterator gb_itr;//Overlap removal between gluons and gluon brothers
                  if(gluon.size() > 0){
                     for(gb_itr=gluon_brothers.begin();gb_itr!=gluon_brothers.end();){
                           bool todelete = false;
                           for(unsigned int i=0; i<gluon.size();i++){
                               if(gluon.at(i) == (*gb_itr)){todelete=true;break;}
                           }
                           if(todelete) gb_itr = gluon_brothers.erase(gb_itr);
                           else ++gb_itr;
                     }
                  }  
//                  cout<<"Mark C"<<endl;
                  sort(gluon_brothers.begin(),gluon_brothers.end(),SortDecayChain);
//                  cout<<"Mark D"<<endl;
                  for(gb_itr=gluon_brothers.begin();gb_itr!=gluon_brothers.end();++gb_itr){
                     if(abs((*gb_itr)->getPDG())>6) cout<<"There is gluon in gluon brothers list!"<<endl;
                  } 
//                  cout<<"Mark E : gluon brothers size "<<gluon_brothers.size()<<endl; 
                  if(gluon_brothers.size()) //Hard gluon radiation found
                     result.push_back(gluon_brothers.at(gluon_brothers.size()-1));
                  else 
                     result.push_back(SelfDecay_RecurseUp(parton_map_itr->first)); 
//                  cout<<"Mark F"<<endl;
               } //Endof quark ancestor

               if (abs(parton_map_itr->first->getPDG()) == 21){  //Gluon Splitting(to quark pairs) Finding, Consider only the primary effect
                   MCParticle * g_tmp = SelfDecay_RecurseDown(parton_map_itr->first);
                   bool gluon_splitting = false;
                   for(unsigned int i =0; i<g_tmp->getDaughters().size();i++){
                        if(abs(g_tmp->getDaughters().at(i)->getPDG())<6 ){gluon_splitting = true;
                        break;}
                   } 

                   if(gluon_splitting){
                      TLorentzVector quark1,quark2;
                      quark1.SetPxPyPzE(g_tmp->getDaughters().at(0)->getMomentum()[0],g_tmp->getDaughters().at(0)->getMomentum()[1],g_tmp->getDaughters().at(0)->getMomentum()[2],g_tmp->getDaughters().at(0)->getEnergy());
                      quark2.SetPxPyPzE(g_tmp->getDaughters().at(1)->getMomentum()[0],g_tmp->getDaughters().at(1)->getMomentum()[1],g_tmp->getDaughters().at(1)->getMomentum()[2],g_tmp->getDaughters().at(1)->getEnergy()); 
                      if(quark1.E()/(quark1.E()+quark2.E())>0.25 && quark2.E()/(quark1.E()+quark2.E())>0.25&&( (quark1.E()>10.0 && quark2.E()>10.0 && quark1.DeltaR(quark2)>0.7)||(quark1.E()>30.0 && quark2.E()>30.0 && quark1.DeltaR(quark2)>0.5) )){
                        if(find(result.begin(),result.end(),g_tmp->getDaughters().at(0)) == result.end() )
                          result.push_back(g_tmp->getDaughters().at(0));
                        if(find(result.begin(),result.end(),g_tmp->getDaughters().at(1)) == result.end() )
                          result.push_back(g_tmp->getDaughters().at(1));
                      }
                      else if(abs(g_tmp->getDaughters().at(0)->getPDG() >=3) && (quark1.E()/(quark1.E()+quark2.E())>0.75 ||quark2.E()/(quark1.E()+quark2.E())>0.75) ){
                          if(g_tmp->getDaughters().at(0)->getEnergy()>g_tmp->getDaughters().at(1)->getEnergy()){
                              if(find(result.begin(),result.end(),g_tmp->getDaughters().at(0)) == result.end()){
                                  result.push_back(g_tmp->getDaughters().at(0));
                              }  
                          }else{
                              if(find(result.begin(),result.end(),g_tmp->getDaughters().at(1)) == result.end()){
                                  result.push_back(g_tmp->getDaughters().at(1));
                              }
                          }
                      }else{
                          if(find(result.begin(),result.end(),SelfDecay_RecurseUp(g_tmp)) == result.end())  
                              result.push_back(SelfDecay_RecurseUp(g_tmp));  
                      }
                   }else if(g_tmp->getDaughters().size() == 2 && g_tmp->getDaughters().at(0)->getPDG() == 21 && g_tmp->getDaughters().at(1)->getPDG() == 21)
                   {
                          TLorentzVector gmomentum1,gmomentum2; 
                          double E1_gluon = g_tmp->getDaughters().at(0)->getEnergy();
                          double E2_gluon = g_tmp->getDaughters().at(1)->getEnergy();
                          gmomentum1.SetPxPyPzE(g_tmp->getDaughters().at(0)->getMomentum()[0],g_tmp->getDaughters().at(0)->getMomentum()[1],g_tmp->getDaughters().at(0)->getMomentum()[2],E1_gluon);
                          gmomentum2.SetPxPyPzE(g_tmp->getDaughters().at(1)->getMomentum()[0],g_tmp->getDaughters().at(1)->getMomentum()[1],g_tmp->getDaughters().at(1)->getMomentum()[2],E2_gluon);
                          if(E1_gluon/(E1_gluon+E2_gluon)>0.3&&E2_gluon/(E1_gluon+E2_gluon)>0.3 && ((gmomentum1.DeltaR(gmomentum2)>1.0 && E1_gluon> 15.0 && E2_gluon>15.0) ||(gmomentum1.DeltaR(gmomentum2)>0.7 && E1_gluon> 30.0 && E2_gluon>30.0) )) {
                              if(find(result.begin(),result.end(),g_tmp->getDaughters().at(0)) == result.end())  
                                 result.push_back(g_tmp->getDaughters().at(0)); 
                              if(find(result.begin(),result.end(),g_tmp->getDaughters().at(1)) == result.end())  
                                 result.push_back(g_tmp->getDaughters().at(1));
                          }
                          else{ 
                              if(find(result.begin(),result.end(),SelfDecay_RecurseUp(g_tmp)) == result.end()) 
                                 result.push_back(g_tmp); 
                          }
                   }
                   else{  
                          if(find(result.begin(),result.end(),SelfDecay_RecurseUp(g_tmp)) == result.end())   
                            result.push_back(g_tmp);
                   }
               }
        }

        return result;
}

bool SortDecayChain(MCParticle *p1, MCParticle *p2){
     if(SelfDecay_RecurseUp(p1) == SelfDecay_RecurseUp(p2)) return true;
     
     bool result=true;
     while(p1->getParents().size()>0 ){
           if(p2 == (p1->getParents().at(0)) ) {result = false; break;}
           else p1 = p1->getParents().at(0);
     } 

     bool result2 =true;
     if(!result){
         while(p2->getParents().size()>0 ){
             if(p1 == (p2->getParents().at(0)) ) {result2 = false; break;}
             else p2 = p2->getParents().at(0);
         }
     }
     if(result||result2) std::cout<<"Sorting Parton: they are not in the same decay chain! Please check!"<<std::endl;
     return result;     
}

map<MCParticle*,vector<MCParticle*> > GetPrimaryParton(vector<MCParticle*> inputPartonList){
//        if(inputPartonList.size() == 0) return NULL;

        map<MCParticle*, vector<MCParticle*> > result;
        map<MCParticle*, vector<MCParticle*> >::iterator itr_result;

        for(unsigned int i=0; i<inputPartonList.size();i++ ){

            MCParticle * tmp  = inputPartonList.at(i);

            if(tmp->getParents().size()>0){
//                tmp = inputPartonList.at(i)->getParents().at(0);
                while(tmp->getParents().size()>0 && tmp->getParents().at(0)->getPDG()!=94 && tmp->getParents().at(0)->getPDG()!=23 && abs(tmp->getParents().at(0)->getPDG())!=24 && tmp->getParents().at(0)->getPDG()!=25){
                    tmp = tmp->getParents().at(0);
                }
            }
//a            else tmp = inputPartonList.at(i);

            bool findthis_ancesoter = false;

            for(itr_result = result.begin();itr_result!=result.end();++itr_result){
                  if( tmp == itr_result->first){
                        (itr_result->second).push_back(inputPartonList.at(i));
                        findthis_ancesoter = true; break;
                  }
            }

            if(!findthis_ancesoter){
                  result[tmp].push_back(inputPartonList.at(i));
            }  
        }

        return result;
}
