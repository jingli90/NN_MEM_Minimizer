#ifndef NN_MEM_Eval_H
#define NN_MEM_Eval_H

class NN_MEM_Eval{
	public:
		TString Filename;
		TString label;

		Long64_t nentries;
		TString cut;
		TString weightExpression;

		int nLayer;
		int nNodes[20];
		TString NeuronType;
		vector<vector<double>*>* weight_init[20];

		vector<TString>* var_input_str;
		vector<TString>* spectator_var_input_str;
		vector<double>* var_max_int;
		vector<double>* var_min_int;
		vector<double>* xL;
		vector<double>* xU;
		int kMode;
		int kCatJets;

		TString configFile;
		mutable ConfigParser cfgParser;
		int nhyp;
		string* shyp;
		int* hyp;
		int* nPointsHyp;
		int* index_hyp;
		int doMinimization;
		string JetChoice;
		int nPermutationJetSyst;

		mutable MultiLepton multiLepton;
		mutable HypIntegrator hypIntegrator;
		mutable Permutations* MEMpermutations;
		mutable ReadGenFlatTree tree;

		mutable simpleNN * nn;

		mutable int num_Calls;
		mutable int num_notZero;

		NN_MEM_Eval(TString File, TString configFile, vector<vector<double>*>* myweight[20], TString myNeuronType, int mynLayer, int mynNodes[20]);
		void SetCut(const TString& cut_){cut = cut_; cout<<"Cut: "<<cut<<endl;}
		void SetWeightExpression(TString expression){weightExpression=expression; cout<<"Weight expression is: \""<< weightExpression<<"\""<<endl;}
		void SetVariables(vector<TString>* var, vector<TString>* spectator){var_input_str=var; spectator_var_input_str=spectator;}
		void SetXLXU(vector<double>* xL_, vector<double>* xU_){xL=xL_; xU=xU_;}
		void SetMinMax(vector<double>* min, vector<double>* max){var_max_int=max; var_min_int=min;}
		void SetModeCat(int kMode_, int kCatJets_){kMode=kMode_; kCatJets=kCatJets_;}
		void Init();
		virtual double Eval_logAvg(const double *) const;
		void Eval_Tree(const double* par, TString filename_tmp, TString foutName);
};

NN_MEM_Eval::NN_MEM_Eval(TString File, TString configFile_, vector<vector<double>*>* myweight[20], TString myNeuronType, int mynLayer, int mynNodes[20]){
	Filename=File;
	//cfgParser.GetConfigFromFile(configFile.Data());
	configFile=configFile_;

	for(int i=0;i<20;i++){
		nNodes[i]=0;
		weight_init[i] = new vector<vector<double>*>;
	}

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

	//cout<<"Init weights:"<<endl;
	for(int i=1;i<=nLayer+1;i++){
		//cout<<"Inner layer "<<i-1<<" to "<<i<<endl;
		for(int j=0;j<=nNodes[i];j++){
			vector<double>* tmp = new vector<double>;
			weight_init[i]->push_back(tmp);
			for(int k=0;k<=nNodes[i-1];k++){
				weight_init[i]->at(j)->push_back(myweight[i]->at(j)->at(k));
				//cout<<"w"<<j<<k<<" = "<<weight_init[i]->at(j)->at(k)<<" ";
			}
			//cout<<endl;
		}
	}

	//nn = new simpleNN(weight_init, NeuronType, nLayer, nNodes);

	var_input_str=new vector<TString>;
	spectator_var_input_str=new vector<TString>;
	var_max_int = new vector<double>;
	var_min_int = new vector<double>;
	xU = new vector<double>;
	xL = new vector<double>;
	kMode=kMEM_TTLL_TopAntitopDecay;
	kCatJets=kCat_3l_2b_2j;

	num_Calls=0;
}

void NN_MEM_Eval::Init(){
	cfgParser.GetConfigFromFile(configFile.Data());
	cfgParser.LoadHypotheses(&nhyp, &shyp, &hyp, &nPointsHyp, &index_hyp);
	doMinimization = cfgParser.valDoMinimization;
	cfgParser.LoadJetChoice(&JetChoice);
	nPermutationJetSyst = cfgParser.nJetSyst;

	cfgParser.LoadIntegrationRange(&multiLepton.JetTFfracmin, &multiLepton.JetTFfracmax, &multiLepton.NeutMaxE);
	hypIntegrator.InitializeIntegrator(&cfgParser);
	MEMpermutations = new Permutations[nhyp];

	tree.InitializeMEMRun(Filename.Data());

	int nNodes_tmp[20];
	for(int i=0;i<20;i++){
		nNodes_tmp[i]=nNodes[i];
	}
	nn = new simpleNN(NeuronType, nLayer, nNodes_tmp);

	num_Calls=0;
}

double NN_MEM_Eval::Eval_logAvg(const double* par) const {
	num_Calls++;
	//if(num_Calls>500){
	//	cout<<"hehe ... "<<endl;
	//	return -log(1e-300);
	//}
	//cout<<"num_Calls="<<num_Calls<<endl;

	nn->SetWeight(par);

	double evt_weight=1;
	//double evt_weight_max=1;
	double sum_weight=0.;
	//evt_weight_max=tree.tInput->GetMaximum(weightExpression.Data());
	double discriminant=0;
	double discriminant_nlog;
	double sum_discriminant_nlog=0;
	vector<double>* var_input_NN;
	vector<double>* var_input_NN_beforeNorm;
	vector<double>* var_output_NN;
	var_input_NN = new vector<double>;
	var_input_NN_beforeNorm = new vector<double>;
	var_output_NN = new vector<double>;
	double * x;
	x=new double[nNodes[nLayer+1]];

	int nEvents=tree.tInput->GetEntries();
	int ih=0;
	int Hypothesis=hyp[ih];
	//bool isStop=false;
	num_notZero=0;
	int num_total=0;
	int num_total1=0;
	for(Long64_t jentry=0; jentry<nEvents;jentry++)
		//for(Long64_t jentry=0; jentry<10;jentry++)
		//for(Long64_t jentry=0; jentry<nEvents && isStop==false;jentry++)
	{
		//if(jentry%1000==0)cout<<"entry: "<<jentry<<endl;
		//cout<<"entry: "<<jentry<<endl;
		tree.ReadMultilepton(jentry, &multiLepton);
		//cout << "Event number: " << tree.mc_event << endl;
		if((tree.catJets==kCatJets) && (Hypothesis==kMode)){
			int initresult = 0;
			MEMpermutations[ih].SetMultiLepton(&multiLepton, &hypIntegrator);
			initresult = MEMpermutations[ih].InitializeHyp(&hypIntegrator, hyp[ih], nPointsHyp[ih], shyp[ih], doMinimization, JetChoice, nPermutationJetSyst); // outside the loop only once			--Jing Li @ 2017-04-26 
			if (initresult==1){
				num_total++;

				int combcheck = 1;
				int permreslep = 1;
				int permresjet = 1;
				int permresbjet = 1;
				int iperm=0;
				bool isCombcheck=false;

				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet1_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet1_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet1_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet2_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet2_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet2_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw1_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw1_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw1_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw2_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw2_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw2_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton1_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton1_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton1_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton2_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton2_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton2_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton3_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton3_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton3_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_mET_ptr->Pt());
				var_input_NN_beforeNorm->push_back(tree.multilepton_mET_ptr->Phi());
				//var_input_NN->clear();
				for(unsigned int ivar=0;ivar<var_input_str->size();ivar++){
					//cout<<"var_input_NN_beforeNorm->at(ivar)="<<var_input_NN_beforeNorm->at(ivar)<<endl;
					double var_tmp=var_input_NN_beforeNorm->at(ivar);
					var_input_NN->push_back(2 * (var_tmp - var_min_int->at(ivar)) / (var_max_int->at(ivar) - var_min_int->at(ivar)) - 1);
				}
				//for(unsigned int ivar=0; ivar<var_input_NN->size(); ivar++)
				//	cout<<"var_input_NN["<<ivar<<"]="<<var_input_NN->at(ivar)<<endl;

				//var_output_NN->clear();
				var_output_NN=nn->Eval(var_input_NN);
				for(unsigned int ivar=0; ivar<var_output_NN->size();ivar++){
					//cout<<"var_output_NN->at("<<ivar<<")="<<var_output_NN->at(ivar)<<endl;
					double var_tmp = (var_output_NN->at(ivar)+1)*0.5 * (xU->at(ivar)-xL->at(ivar)) + xL->at(ivar);
					x[ivar]=var_tmp;
					//cout<<"x["<<ivar<<"]="<<x[ivar]<<endl;
				}

				double discriminant_tmp = 0;
				double discriminant_max = 0;

				//cout << "Start looping" << endl;

				MEMpermutations[ih].multiLepton.DoSort(&MEMpermutations[ih].multiLepton.Bjets);
				do{
					MEMpermutations[ih].multiLepton.DoSort(&MEMpermutations[ih].multiLepton.Jets);
					do{
						MEMpermutations[ih].multiLepton.DoSort(&MEMpermutations[ih].multiLepton.Leptons);
						do
						{
							combcheck = MEMpermutations[ih].multiLepton.CheckPermutationHyp(Hypothesis);
							//cout<<"combcheck="<<combcheck<<endl;
							//cout<<"iperm="<<iperm<<" combcheck="<<combcheck<<endl;
							if (combcheck) {
								isCombcheck=true;
								//cout<<"leptonSize="<<MEMpermutations[ih].multiLepton.Leptons.size()<<endl;
								//for (unsigned int il=0; il<MEMpermutations[ih].multiLepton.Leptons.size(); il++) cout << " Lepton"<< il<<"Id="<<MEMpermutations[ih].multiLepton.Leptons.at(il).Id; cout<<endl;
								//for (unsigned int ib=0; ib<MEMpermutations[ih].multiLepton.Bjets.size(); ib++) cout << " Bjet"<< ib<<"Pt="<<MEMpermutations[ih].multiLepton.Bjets.at(ib).P4.Pt(); cout<<endl;
								//for (unsigned int ij=0; ij<MEMpermutations[ih].multiLepton.Jets.size(); ij++) cout << " Jet"<< ij<<"Pt="<<MEMpermutations[ih].multiLepton.Jets.at(ij).P4.Pt(); cout<<endl;
								MEMpermutations[ih].multiLepton.FillParticlesHypothesis(Hypothesis, &(hypIntegrator.meIntegrator));
								MEMpermutations[ih].multiLepton.SwitchJetSyst(0);

								discriminant_tmp=hypIntegrator.meIntegrator->Eval(x);
								//cout<<"hypIntegrator.meIntegrator->Eval(x)="<<discriminant_tmp<<endl;
								if(discriminant_tmp>discriminant_max)
									discriminant_max=discriminant_tmp;
							}

							iperm++;

							if(MEMpermutations[ih].doPermutationLep){
								if (Hypothesis!=kMEM_TTbar_TopAntitopSemiLepDecay)
									permreslep = MEMpermutations[ih].multiLepton.DoPermutation(&MEMpermutations[ih].multiLepton.Leptons);
								else permreslep = MEMpermutations[ih].multiLepton.DoPermutationLinear(&MEMpermutations[ih].multiLepton.Leptons);
							}
							//cout<<"permreslep="<<permreslep<<endl;

						}while (MEMpermutations[ih].doPermutationLep && permreslep);

						if(MEMpermutations[ih].doPermutationJet){
							if (MEMpermutations[ih].multiLepton.kCatJets!=kCat_3l_2b_1j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_3l_1b_1j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_2b_3j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_1b_3j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_2b_2j)
								permresjet = MEMpermutations[ih].multiLepton.DoPermutation(&MEMpermutations[ih].multiLepton.Jets);
							else
								permresjet = MEMpermutations[ih].multiLepton.DoPermutationMissingJet("jet");
						}
						//cout<<"permresjet="<<permresjet<<endl;

					}while(MEMpermutations[ih].doPermutationJet && permresjet);

					if(MEMpermutations[ih].doPermutationBjet){
						if(Hypothesis==kMEM_TLLJ_TopLepDecay)
							permresbjet = MEMpermutations[ih].multiLepton.DoPermutationLinear(&MEMpermutations[ih].multiLepton.Bjets);
						else if(MEMpermutations[ih].multiLepton.kCatJets!=kCat_3l_1b_2j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_3l_1b_1j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_4l_1b && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_1b_4j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_1b_3j)
							permresbjet = MEMpermutations[ih].multiLepton.DoPermutation(&MEMpermutations[ih].multiLepton.Bjets);
						else
							permresbjet = MEMpermutations[ih].multiLepton.DoPermutationMissingJet("bjet");
					}
					//cout<<"permresbjet"<<permresbjet<<endl;

				}while (MEMpermutations[ih].doPermutationBjet && permresbjet);

				discriminant=discriminant_max;

				if(isCombcheck==true){
					num_total1++;

					if(discriminant>0){
						discriminant_nlog = -log(discriminant);
						num_notZero++;
					}
					else
						discriminant_nlog = -log(1e-300);

					evt_weight=tree.weight;
					sum_discriminant_nlog = sum_discriminant_nlog + discriminant_nlog * evt_weight;
					sum_weight=sum_weight+evt_weight;
				}

				//isStop=true;

				var_input_NN_beforeNorm->clear();
				var_input_NN->clear();
				var_output_NN->clear();
			}
		}
	}

	sum_discriminant_nlog=sum_discriminant_nlog/sum_weight;
	cout<<"num_Calls="<<setw(10)<<num_Calls<<" ";
	cout<<"num_total="<<setw(10)<<num_total<<" ";
	cout<<"num_total1="<<setw(10)<<num_total1<<" ";
	cout<<"num_notZero="<<setw(10)<<num_notZero<<"  ";
	cout<<"average_discriminant_nlog="<<sum_discriminant_nlog<<endl;
	return sum_discriminant_nlog;
}

void NN_MEM_Eval::Eval_Tree(const double* par, TString filename_tmp, TString foutName) {
	TFile * fout = new TFile(foutName,"RECREATE");
	TTree * tout = new TTree();

	tree.InitializeMEMRun(filename_tmp.Data());
	//tout = tree.tInput->CloneTree(0);

	nn->SetWeight(par);

	double evt_weight=1;
	double discriminant=0;
	double discriminant_nlog;
	vector<double>* var_input_NN;
	vector<double>* var_input_NN_beforeNorm;
	vector<double>* var_output_NN;
	var_input_NN = new vector<double>;
	var_input_NN_beforeNorm = new vector<double>;
	var_output_NN = new vector<double>;
	double * x;
	x=new double[nNodes[nLayer+1]];

	tout->Branch("discriminant",&discriminant,"discriminant/D");
	tout->Branch("discriminant_nlog",&discriminant_nlog,"discriminant_nlog/D");
	tout->Branch("weight",&evt_weight,"weight/D");

	vector<double>* o[20];
	for(int i=0;i<20;i++){
		o[i] = new vector<double>;
	}
	for(int i=0;i<=nLayer+1;i++){
		o[i]->clear();
		for(int j=0;j<=nNodes[i];j++){
			o[i]->push_back(-2.);
		}
	}
	for(int i=0;i<=nLayer+1;i++){
		for(int j=0;j<=nNodes[i];j++){
			tout->Branch(Form("l%dn%d",i,j), &(o[i]->at(j)), Form("l%dn%d/D",i,j) );
		}
	}


	int nEvents=tree.tInput->GetEntries();
	int ih=0;
	int Hypothesis=hyp[ih];
	for(Long64_t jentry=0; jentry<nEvents;jentry++)
		//for(Long64_t jentry=0; jentry<100;jentry++)
	{
		//if(jentry%1000==0)cout<<"entry: "<<jentry<<endl;
		//cout<<"entry: "<<jentry<<endl;
		tree.ReadMultilepton(jentry, &multiLepton);
		//cout << "Event number: " << tree.mc_event << endl;
		if((tree.catJets==kCatJets) && (Hypothesis==kMode)){
			int initresult = 0;
			MEMpermutations[ih].SetMultiLepton(&multiLepton, &hypIntegrator);
			initresult = MEMpermutations[ih].InitializeHyp(&hypIntegrator, hyp[ih], nPointsHyp[ih], shyp[ih], doMinimization, JetChoice, nPermutationJetSyst); // outside the loop only once			--Jing Li @ 2017-04-26 
			if (initresult==1){

				int combcheck = 1;
				int permreslep = 1;
				int permresjet = 1;
				int permresbjet = 1;
				int iperm=0;
				bool isCombcheck=false;

				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet1_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet1_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet1_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet2_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet2_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Bjet2_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw1_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw1_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw1_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw2_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw2_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_JetClosestMw2_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton1_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton1_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton1_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton2_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton2_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton2_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton3_P4_ptr->E());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton3_P4_ptr->Theta());
				var_input_NN_beforeNorm->push_back(tree.multilepton_Lepton3_P4_ptr->Phi());
				var_input_NN_beforeNorm->push_back(tree.multilepton_mET_ptr->Pt());
				var_input_NN_beforeNorm->push_back(tree.multilepton_mET_ptr->Phi());
				//var_input_NN->clear();
				for(unsigned int ivar=0;ivar<var_input_str->size();ivar++){
					//cout<<"var_input_NN_beforeNorm->at(ivar)="<<var_input_NN_beforeNorm->at(ivar)<<endl;
					double var_tmp=var_input_NN_beforeNorm->at(ivar);
					var_input_NN->push_back(2 * (var_tmp - var_min_int->at(ivar)) / (var_max_int->at(ivar) - var_min_int->at(ivar)) - 1);
				}
				//for(unsigned int ivar=0; ivar<var_input_NN->size(); ivar++)
				//	cout<<"var_input_NN["<<ivar<<"]="<<var_input_NN->at(ivar)<<endl;

				//var_output_NN->clear();
				var_output_NN=nn->Eval(var_input_NN);
				for(unsigned int ivar=0; ivar<var_output_NN->size();ivar++){
					//cout<<"var_output_NN->at("<<ivar<<")="<<var_output_NN->at(ivar)<<endl;
					double var_tmp = (var_output_NN->at(ivar)+1)*0.5 * (xU->at(ivar)-xL->at(ivar)) + xL->at(ivar);
					x[ivar]=var_tmp;
					//cout<<"x["<<ivar<<"]="<<x[ivar]<<endl;
				}

				double discriminant_tmp = 0;
				double discriminant_max = 0;

				//cout << "Start looping" << endl;

				MEMpermutations[ih].multiLepton.DoSort(&MEMpermutations[ih].multiLepton.Bjets);
				do{
					MEMpermutations[ih].multiLepton.DoSort(&MEMpermutations[ih].multiLepton.Jets);
					do{
						MEMpermutations[ih].multiLepton.DoSort(&MEMpermutations[ih].multiLepton.Leptons);
						do
						{
							combcheck = MEMpermutations[ih].multiLepton.CheckPermutationHyp(Hypothesis);
							//cout<<"combcheck="<<combcheck<<endl;
							//cout<<"iperm="<<iperm<<" combcheck="<<combcheck<<endl;
							if (combcheck) {
								isCombcheck=true;
								//cout<<"leptonSize="<<MEMpermutations[ih].multiLepton.Leptons.size()<<endl;
								//for (unsigned int il=0; il<MEMpermutations[ih].multiLepton.Leptons.size(); il++) cout << " Lepton"<< il<<"Id="<<MEMpermutations[ih].multiLepton.Leptons.at(il).Id; cout<<endl;
								//for (unsigned int ib=0; ib<MEMpermutations[ih].multiLepton.Bjets.size(); ib++) cout << " Bjet"<< ib<<"Pt="<<MEMpermutations[ih].multiLepton.Bjets.at(ib).P4.Pt(); cout<<endl;
								//for (unsigned int ij=0; ij<MEMpermutations[ih].multiLepton.Jets.size(); ij++) cout << " Jet"<< ij<<"Pt="<<MEMpermutations[ih].multiLepton.Jets.at(ij).P4.Pt(); cout<<endl;
								MEMpermutations[ih].multiLepton.FillParticlesHypothesis(Hypothesis, &(hypIntegrator.meIntegrator));
								MEMpermutations[ih].multiLepton.SwitchJetSyst(0);

								discriminant_tmp=hypIntegrator.meIntegrator->Eval(x);
								//cout<<"hypIntegrator.meIntegrator->Eval(x)="<<discriminant_tmp<<endl;
								if(discriminant_tmp>discriminant_max){
									discriminant_max=discriminant_tmp;

									for(int i=0;i<=nLayer+1;i++){
										o[i]->clear();
										for(int j=0;j<=nNodes[i];j++){
											o[i]->push_back(nn->o[i]->at(j));
										}
									}
								}
							}

							iperm++;

							if(MEMpermutations[ih].doPermutationLep){
								if (Hypothesis!=kMEM_TTbar_TopAntitopSemiLepDecay)
									permreslep = MEMpermutations[ih].multiLepton.DoPermutation(&MEMpermutations[ih].multiLepton.Leptons);
								else permreslep = MEMpermutations[ih].multiLepton.DoPermutationLinear(&MEMpermutations[ih].multiLepton.Leptons);
							}
							//cout<<"permreslep="<<permreslep<<endl;

						}while (MEMpermutations[ih].doPermutationLep && permreslep);

						if(MEMpermutations[ih].doPermutationJet){
							if (MEMpermutations[ih].multiLepton.kCatJets!=kCat_3l_2b_1j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_3l_1b_1j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_2b_3j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_1b_3j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_2b_2j)
								permresjet = MEMpermutations[ih].multiLepton.DoPermutation(&MEMpermutations[ih].multiLepton.Jets);
							else
								permresjet = MEMpermutations[ih].multiLepton.DoPermutationMissingJet("jet");
						}
						//cout<<"permresjet="<<permresjet<<endl;

					}while(MEMpermutations[ih].doPermutationJet && permresjet);

					if(MEMpermutations[ih].doPermutationBjet){
						if(Hypothesis==kMEM_TLLJ_TopLepDecay)
							permresbjet = MEMpermutations[ih].multiLepton.DoPermutationLinear(&MEMpermutations[ih].multiLepton.Bjets);
						else if(MEMpermutations[ih].multiLepton.kCatJets!=kCat_3l_1b_2j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_3l_1b_1j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_4l_1b && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_1b_4j && MEMpermutations[ih].multiLepton.kCatJets!=kCat_2lss_1b_3j)
							permresbjet = MEMpermutations[ih].multiLepton.DoPermutation(&MEMpermutations[ih].multiLepton.Bjets);
						else
							permresbjet = MEMpermutations[ih].multiLepton.DoPermutationMissingJet("bjet");
					}
					//cout<<"permresbjet"<<permresbjet<<endl;

				}while (MEMpermutations[ih].doPermutationBjet && permresbjet);

				discriminant=discriminant_max;

				if(isCombcheck==true){
					if(discriminant>0){
						discriminant_nlog = -log(discriminant);
					}
					else
						discriminant_nlog = -log(1e-300);

					evt_weight=tree.weight;

					tout->Fill();
				}

				//isStop=true;

				var_input_NN_beforeNorm->clear();
				var_input_NN->clear();
				var_output_NN->clear();

			}
		}
	}
	//tout->Print();
	fout->cd();
	tout->Write("Tree");

}

#endif
