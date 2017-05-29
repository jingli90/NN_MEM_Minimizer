#ifndef NN_MEM_h
#define NN_MEM_h

#include <iostream>
using namespace std;
#include <vector>
#include <time.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>

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

#include "simpleNN.h"
#include "NN_MEM_Eval.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MinosError.h"
#include "Minuit2/ContoursError.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Integrator.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

class NN_MEM{
	public:
		TString Filename;
		TString label;
		TString Filename_other[10];
		TString Filelabel_other[10];
		TString outFilename;
		int nOther;

		Long64_t nentries;
		TString cut;
		TString weightExpression;

		vector<vector<double>*>* weight[20];
		vector<vector<double>*>* weight_init[20];
		std::vector<double> init_par;
		std::vector<double> par;
		std::vector<double> par_old;
		int total;

		int nEpochs;
		int nLayer;
		int nNodes[20];
		TString NeuronType;
		int num_outputNN;

		NN_MEM_Eval * theEval;

		vector<TString>* var_input_str;
		vector<TString>* spectator_var_input_str;
		vector<double>* var_max_int;
		vector<double>* var_min_int;
		vector<double>* xL;
		vector<double>* xU;
		int kMode;
		int kCatJets;

		//double discriminant;
		double duration;

		TString configFile;

		ROOT::Minuit2::Minuit2Minimizer* minimizer;

		NN_MEM(TString File, TString outfilename, TString label_, TString configFile_);
		virtual void SetNLayer(int nlayer);
		virtual void SetNNodes(int i, int nnodes);
		void SetNum_outputNN(int N){num_outputNN=N; nNodes[nLayer+1]=N; cout<<"number of input to MEM: "<< N <<endl;}
		void SetNEpochs(int N){nEpochs=N;cout<<"nEpochs = "<< N <<endl;}
		void SetNeuronType(TString myNeuronType){NeuronType=myNeuronType; cout<<"NeuronType="<<NeuronType<<endl;}
		void SetCut(const TString& cut_){cut = cut_; cout<<"Cut: "<<cut<<endl;}
		void SetWeightExpression(TString expression){weightExpression=expression; cout<<"Weight expression is: \""<< weightExpression<<"\""<<endl;}

		void AddVariable(TString var){cout<<"Add input variable: "<<var<<endl; var_input_str->push_back(var);}
		void AddSpectator(TString var){cout<<"Add spectator: "<<var<<endl; spectator_var_input_str->push_back(var);}
		void SetXLXU(vector<double>* xL_, vector<double>* xU_){xL=xL_; xU=xU_;}
		void SetMinMax(vector<double>* min, vector<double>* max){var_max_int=max; var_min_int=min;}
		void SetModeCat(int kMode_, int kCatJets_){kMode=kMode_; kCatJets=kCatJets_;}

		void AddOtherFiles(TString other, TString label){Filename_other[nOther]=other; Filelabel_other[nOther]=label; nOther++; cout<<"will eval also on "<<other<<endl;}

		virtual void myana();
		virtual void InitParameters();
		virtual void InitFromFile(TString weightInit);
		virtual void doTraining();
		virtual void Eval_Multi(int N);
		virtual void Eval_Tree();

};

NN_MEM::NN_MEM(TString File, TString outfilename, TString label_, TString configFile_){
	cout<<endl;
	cout<<"Welcome to the training ..."<<endl;
	cout<<"Filename: "<<File<<endl;
	Filename=File;
	outFilename = outfilename;
	label=label_;
	configFile=configFile_;
	nOther=0;

	for(int i=0;i<20;i++){
		nNodes[i]=0;
		weight[i] = new vector<vector<double>*>;
		weight_init[i] = new vector<vector<double>*>;
	}
	NeuronType="tanh";
	nEpochs=-1;

	var_input_str = new vector<TString>;
	spectator_var_input_str = new vector<TString>;
	var_max_int = new vector<double>;
	var_min_int = new vector<double>;
	xU = new vector<double>;
	xL = new vector<double>;
	kMode=kMEM_TTLL_TopAntitopDecay;
	kCatJets=kCat_3l_2b_2j;

	minimizer = new ROOT::Minuit2::Minuit2Minimizer( ROOT::Minuit2::kMigrad );
	minimizer->SetPrintLevel(0);
}

void NN_MEM::SetNLayer(int nlayer){
	cout<<endl;
	nLayer=nlayer;
	if(nlayer>18){
		cout<<"Too much inner layers. Set number of inner layers = 18"<<endl;
		nLayer=18;
	}
	cout<<"Number of inner layers = "<<nLayer<<endl;
} 

void NN_MEM::SetNNodes(int i, int nnodes){
	if(i>=nLayer){ 
		cout<<"Inner layer "<<i<<" does not exit"<<endl;
	}
	nNodes[i]=nnodes;
	cout<<"Number of nodes in inner layer "<<i<<" is "<<nNodes[i]<<endl;
}

void NN_MEM::InitParameters(){
	nNodes[0]=var_input_str->size();
	cout<<endl;

	cout<<"Init weights:"<<endl;
	gRandom = new TRandom3();
	gRandom->SetSeed(0);
	int index=0.;
	for(int i=1;i<=nLayer+1;i++){
		cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			vector<double>* tmp = new vector<double>;
			vector<double>* tmp_init = new vector<double>;
			weight[i]->push_back(tmp);
			weight_init[i]->push_back(tmp_init);
			for(int k=0;k<=nNodes[i-1];k++){
				if(j==0){
					weight[i]->at(j)->push_back(0);
					weight_init[i]->at(j)->push_back(0);
				}
				else{
					weight[i]->at(j)->push_back(2*gRandom->Rndm()-1);
					weight_init[i]->at(j)->push_back(weight[i]->at(j)->at(k));
					init_par.push_back(weight_init[i]->at(j)->at(k));
					par.push_back(weight_init[i]->at(j)->at(k));
					par_old.push_back(weight_init[i]->at(j)->at(k));
					index++;
				}
				cout<<"w"<<j<<k<<" = "<<weight[i]->at(j)->at(k)<<" ";
			}
			cout<<endl;
		}
	}
	total=index;

	theEval = new NN_MEM_Eval(Filename, configFile, weight_init, NeuronType, nLayer, nNodes);
	theEval->SetCut(cut);
	theEval->SetWeightExpression(weightExpression);
	theEval->SetVariables(var_input_str, spectator_var_input_str);
	theEval->SetXLXU(xL, xU);
	theEval->SetMinMax(var_max_int, var_min_int);
	theEval->SetModeCat(kMode, kCatJets);
	theEval->Init();
}

void NN_MEM::InitFromFile(TString weightInit){

	cout<<"Init weights from "<<weightInit<<endl;

	ifstream fweightInit;
	fweightInit.open(weightInit.Data());
	string line;

	init_par.clear();
	par.clear();
	par_old.clear();
	double weight_tmp=0;
	for(int i=1;i<=nLayer+1;i++){
		fweightInit>>line;
		//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		cout<<line<<endl;
		weight[i]->clear();
		weight_init[i]->clear();
		for(int j=0;j<=nNodes[i];j++){
			vector<double>* tmp = new vector<double>;
			vector<double>* tmp_init = new vector<double>;
			weight[i]->push_back(tmp);
			weight_init[i]->push_back(tmp_init);
			for(int k=0;k<=nNodes[i-1];k++){
				if(j==0){
					weight[i]->at(j)->push_back(0);
					weight_init[i]->at(j)->push_back(0);
				}
				else{
					fweightInit>>line;
					fweightInit>>line;
					fweightInit>>line;
					weight_tmp=strtod(line.c_str(), NULL);
					cout<<"w"<<j<<k<<" = "<<weight_tmp<<" ";

					weight[i]->at(j)->push_back(weight_tmp);
					weight_init[i]->at(j)->push_back(weight[i]->at(j)->at(k));
					init_par.push_back(weight_init[i]->at(j)->at(k));
					par.push_back(weight_init[i]->at(j)->at(k));
					par_old.push_back(weight_init[i]->at(j)->at(k));
				}
			}
			if(j!=0) cout<<endl;
		}
	}


}

void NN_MEM::doTraining(){
	cout<<endl;
	cout<<"start training ..."<<endl;

	clock_t start, finish;
	start = clock();

	//double par_test[total];
	//for(int i=0;i<total;i++)
		//par_test[i]=par[i];
	//double average_discriminant_nlog=theEval->Eval_logAvg(par_test);
	//cout<<"average_discriminant_nlog="<<average_discriminant_nlog<<endl;

	ROOT::Math::Functor FunctorHyp(theEval, &NN_MEM_Eval::Eval_logAvg, total);
	minimizer->SetFunction(FunctorHyp);
	if(nEpochs>0)minimizer->SetMaxFunctionCalls(nEpochs);
	minimizer->SetPrintLevel(1);
	for(int i=0; i<total; i++)
		//minimizer->SetVariable(i,Form("%d", i), init_par[i], 0.002);
		//minimizer->SetVariable(i,Form("%d", i), init_par[i], 0.005);
		minimizer->SetVariable(i,Form("%d", i), init_par[i], 0.01);
	bool isSuccess=minimizer->Minimize();
	cout<<"isSuccess="<<isSuccess<<endl;

	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC;
	//cout<<"CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
	cout<<"duration="<<duration<<" sec"<<endl;
	cout<<endl;

	const double *xs = minimizer->X();
	par.clear();
	for(int l=0; l<total; l++){
		par.push_back(xs[l]);
	}
	cout<<"Final weights: "<<endl;
	int i_xs=0;
	for(int i=1;i<=nLayer+1;i++){
		cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		weight[i]->clear();
		for(int j=0;j<=nNodes[i];j++){
			vector<double>* tmp = new vector<double>;
			weight[i]->push_back(tmp);
			for(int k=0;k<=nNodes[i-1];k++){
				if(j==0){
					weight[i]->at(j)->push_back(0);
				}
				else{
					weight[i]->at(j)->push_back(xs[i_xs]);
					cout<<"w"<<j<<k<<" = "<<xs[i_xs]<<" ";
					i_xs++;
				}
			}
			if(j!=0) cout<<endl;
		}
	}

}

void NN_MEM::Eval_Multi(int N){
	cout<<endl;
	cout<<"start init many times ..."<<endl;

	clock_t start, finish;
	start = clock();

	double par_test[total];
	for(int n_init=0; n_init<N; n_init++){
		cout<<"try "<<n_init<<endl;
		for(int i=0;i<total;i++)
			par_test[i]=2*gRandom->Rndm()-1;
		double average_discriminant_nlog=theEval->Eval_logAvg(par_test);
		cout<<"average_discriminant_nlog="<<average_discriminant_nlog<<endl;
	}

	finish = clock();
	duration = (double)(finish-start)/CLOCKS_PER_SEC;
	//cout<<"CLOCKS_PER_SEC="<<CLOCKS_PER_SEC<<endl;
	cout<<"duration="<<duration<<" sec"<<endl;
	cout<<endl;

}

void NN_MEM::Eval_Tree(){
	cout<<endl;
	cout<<"start evaluation ..."<<endl;

	double par_test[total];
	for(int i=0;i<total;i++)
		par_test[i]=par[i];

	theEval->Eval_Tree(par_test, Filename, Form("%s_%s.root",outFilename.Data(), label.Data()));


	for(int i=0; i<nOther; i++){
		cout<<"start evaluation on "<< Filename_other[i] <<endl;
		theEval->Eval_Tree(par_test, Filename_other[i], Form("%s_%s.root",outFilename.Data(), Filelabel_other[i].Data()));
	}

}

void NN_MEM::myana(){
	InitParameters();
	InitFromFile("weightInit.txt");
	doTraining();
	//Eval_Multi(1000);
	Eval_Tree();
}


#endif
