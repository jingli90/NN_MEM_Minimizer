#ifndef simpleNN_h
#define simpleNN_h

#include <iostream>
using namespace std;
#include <vector>
#include <time.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>

#include "base_function.h"

class simpleNN{
	public:
		int nLayer;
		int nNodes[20];
		TString NeuronType;
		vector<double>* var;
		vector<double>* net[20];
		vector<double>* o[20];
		vector<vector<double>*>* weight[20];
		vector<double>* discriminant;

		simpleNN(TString myNeuronType, int mynLayer, int mynNodes[20]);
		void SetWeight(const double* par);
		vector<double>* Eval(const vector<double>*) const;
};

simpleNN::simpleNN(TString myNeuronType, int mynLayer, int mynNodes[20]){
	//cout<<endl;
	cout<<"init the neural network ...";

	for(int i=0;i<20;i++){
		nNodes[i]=0;
		net[i] = new vector<double>;
		o[i] = new vector<double>;
		weight[i] = new vector<vector<double>*>;
	}
	var = new vector<double>;
	discriminant = new vector<double>;
	cout<<endl;

	NeuronType=myNeuronType;

	nLayer=mynLayer;
	if(mynLayer>18){
		cout<<"Too much inner layers. Set number of inner layers = 18"<<endl;
		nLayer=18;
	}
	for(int i=0;i<=mynLayer;i++){
		nNodes[i]=mynNodes[i];
	}
	nNodes[nLayer+1]=mynNodes[nLayer+1];
	//cout<<"Number of inner layers = "<<nLayer<<endl;

}

void simpleNN::SetWeight(const double* par){
	cout<<"Set weights..."<<endl;

	int index=0;
	for(int i=1;i<=nLayer+1;i++){
		cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			vector<double>* tmp = new vector<double>;
			weight[i]->push_back(tmp);
			weight[i]->at(j)->clear();
			for(int k=0;k<=nNodes[i-1];k++){
				if(j==0){
					weight[i]->at(j)->push_back(0);
				}
				else{
					weight[i]->at(j)->push_back(par[index]);
					index++;
				}
				if(j!=0)
				cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" ";
			}
			cout<<endl;
		}
	}

}

vector<double>* simpleNN::Eval(const std::vector<double>* par) const{

	var->clear();
	net[0]->clear();
	o[0]->clear();

	net[0]->push_back(1);
	o[0]->push_back(1);

	for(int ivar=0;ivar<nNodes[0];ivar++){
		var->push_back(par->at(ivar));
		o[0]->push_back(var->at(ivar));
		//cout<<var_input->at(ivar)<<"_nrm="<<var->at(ivar)<<endl;
	}

	for(int i=1;i<=nLayer+1;i++){
		//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		net[i]->clear();
		o[i]->clear();
		net[i]->push_back(1);
		o[i]->push_back(1);
		for(int j=1;j<=nNodes[i];j++){
			double net_tmp=0;
			double o_tmp=0;
			for(int k=0;k<=nNodes[i-1];k++){
				net_tmp=net_tmp + weight[i]->at(j)->at(k) * o[i-1]->at(k);
				//cout<<"w"<<j<<k<<endl;
			}
			o_tmp=activation_function(net_tmp, NeuronType);
			net[i]->push_back(net_tmp);
			o[i]->push_back(o_tmp);
		}
	}

	discriminant->clear();
	for(int i=1;i<=nNodes[nLayer+1];i++){
		discriminant->push_back(o[nLayer+1]->at(i));
		//cout<<"discriminant["<<i-1<<"]="<<discriminant->at(i-1)<<endl;
	}

	return discriminant;
}

#endif
