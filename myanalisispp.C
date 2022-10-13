#include "TList.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TColor.h"
#include "TMath.h"
#include <typeinfo>
#include "TStyle.h"
#include "THistPainter.h"
#include "TLatex.h"
#include "my_tools.C"
#include "my_settings.C"


/*
   root
   .L myanalisispp.C
   MCforpp13("13dNdpT_output.root","arbolespp")
   MCforpp502("502dNdpT_output.root","arbolespp")	

*/
Double_t pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852;

void LoadLibs();
void SetStyle(Bool_t graypalette=kFALSE);
Float_t meanpt, cent, multiplicity;
TTree * fTree;

void MCforpp13(const Char_t * inFilepp, const Char_t * outDirpp){

	cout<<" Creating directory for dNdptppMc..."<<endl;
	CreateDir(outDirpp); 

	cout<<" Opening root file..."<<endl;
	//Archivo de salida
	TFile * outFiletree = 0;
	outFiletree = new TFile(Form("%s/newcutTreeV0Mpp13.root",outDirpp),"RECREATE");

	//<dN/deta>*deltaeta for SPD and V0M at 13 Tev
	Double_t meandNSPD[10]={86.56,71.36,62.24,54.56,46.88,39.2,31.2,23.04,14.44,4.65};
	Double_t truemeandNSPD[10]={4.65,14.44,23.04,31.2,39.2,46.88,54.56,62.24,71.36,86.56};
	Double_t meandNV0M[10]={42.56,32.80,26.72,22.88,20.16,16.96,13.54,10.91,7.90,4.06};
	Double_t truemeandNV0M[10]={4.06,7.90,10.91,13.54,16.96,20.26,22.88,26.72,32.80,42.56};
	//Numero de entradas para las dist de poisson
        Double_t NentriesSPD[10]={5,47,108,306,827,2082,4894,10694,22356,50000};
	Double_t trueNentriesSPD[10]={50000,22536,10694,4894,2082,827,306,108,47,5};
	Double_t NentriesV0M[10]={166,666,833,833,833,1666,1666,1666,3333,5000};
	//Porcentajes eventos
	Double_t percentsSPD[10] = {0.00006,0.00052,0.00119,0.00336,0.00906,0.0228,0.0536,0.11711,0.2448,0.5475};
	Double_t truepercentsSPD[10]={0.5475,0.2448,0.11711,0.0536,0.0228,0.00906,0.00336,0.00119,0.00052,0.00006};
        Double_t percentsV0M[10] = {0.01,0.04,0.05,0.05,0.05,0.1,0.1,0.1,0.2,0.3};
	Double_t truepercentsV0M[10] = {0.3,0.2,0.1,0.1,0.1,0.5,0.5,0.5,0.04,0.01};
	Int_t numberevents = 1700000;//numero de eventos 
        gRandom->SetSeed(0);

	Double_t N13dataerror[10]={0.26,0.28,0.34,0.40,0.5,0.5,0.6,0.7,0.8,1.1};
	
	//Creamos el árbol
	fTree = 0;
	fTree = new TTree("tree","Event data");
	fTree->Branch("meanpt",&meanpt,"meanpt/F");
	fTree->Branch("multiplicity",&multiplicity,"multiplicity/F");
	fTree->Branch("cent",&cent,"cent/F");
	//archivoentrada
	TFile * fpp = 0;
	fpp = TFile::Open(inFilepp);

	//Iniciamos histogramas
	TH1 * dNdptrandom[10];
	TH1 * hpoissonpp[10];
	TF1 * fhpp[10];
	TH1 * yrandompp[10];
	TH1 * hRatpp[10];
        TH1 * hRat[10];
	TH1 * hMCpp[10];
	TH1 * hdatapp[10];
	TH1 * hMC[10];
	TH1 * hdata[10];
	//Simulamos <dN/deta>_pp como distribuciones de Poisson
	for (Int_t s=0 ; s<10; s++)
    	   {  
       
      		hpoissonpp[s] = 0;
      		fhpp[s] = 0;
          	hpoissonpp[s] = new TH1F(Form("h%d",s),"Poisson",100,0,100);
          	fhpp[s] = new TF1(Form("fh%d",s),"TMath::Poisson(x,[0])",0,100);
          	fhpp[s]->SetParameter(0,truemeandNV0M[s]);
          	hpoissonpp[s]->Sumw2();
          	hpoissonpp[s]->FillRandom(Form("fh%d",s), truepercentsV0M[s]*numberevents);
          }

	// Load necessary libraries
	LoadLibs();
	// Set the default style
	SetStyle();
        //Graficando las dist de poisson para pp
	TCanvas * poissonpp = new TCanvas();
	poissonpp->SetTitle("Poisson");
	poissonpp->SetLogy(1);
	auto fpoisspp = new TLegend(0.60,0.65,0.90,0.90);
	fpoisspp->SetFillStyle(0);
	fpoisspp->SetTextSize(0.04);
	fpoisspp->SetHeader("V0M pp (#surd s_{NN} = 13 TeV)","C");
	fpoisspp->AddEntry(hpoissonpp[0],"X (#lambda = 4.06)" ,"p");
	fpoisspp->AddEntry(hpoissonpp[1],"IX (#lambda = 7.90)","p");
	fpoisspp->AddEntry(hpoissonpp[2],"VIII (#lambda = 10.91)" ,"p");
	fpoisspp->AddEntry(hpoissonpp[3],"VII (#lambda = 13.54)","p");
	fpoisspp->AddEntry(hpoissonpp[4],"VI (#lambda = 16.96)" ,"p");
	fpoisspp->AddEntry(hpoissonpp[5],"V (#lambda = 20.26)","p");
	fpoisspp->AddEntry(hpoissonpp[6],"IV (#lambda = 22.88)" ,"p");
	fpoisspp->AddEntry(hpoissonpp[7],"III (#lambda = 26.72)" ,"p");
	fpoisspp->AddEntry(hpoissonpp[8],"II (#lambda = 32.80)" ,"p");
	fpoisspp->AddEntry(hpoissonpp[9],"I (#lambda = 42.56)" ,"p");
	
	//Titulos a los ejes para las distribuciones de Poisson	
	hpoissonpp[0]->GetXaxis()->SetTitle("#it{N}_{ch}");
	hpoissonpp[0]->GetYaxis()->SetTitle("Number of events");
	hpoissonpp[0]->GetXaxis()->SetTitleSize(0.08);
	hpoissonpp[0]->GetYaxis()->SetTitleSize(0.08);
	hpoissonpp[0]->GetXaxis()->SetTitleOffset(.6);
	hpoissonpp[0]->GetYaxis()->SetTitleOffset(.6);

	hpoissonpp[0]->SetAxisRange(1, 800, "Y");
	//Estilo de las distribuciones
	hpoissonpp[0]->SetMarkerStyle(kFullDiamond);
	hpoissonpp[1]->SetMarkerStyle(kFullDoubleDiamond);
	hpoissonpp[2]->SetMarkerStyle(kFullCrossX);
	hpoissonpp[3]->SetMarkerStyle(kFullCircle);
	hpoissonpp[4]->SetMarkerStyle(kFullStar);
	hpoissonpp[5]->SetMarkerStyle(kFullSquare);
	hpoissonpp[6]->SetMarkerStyle(kFullTriangleUp);
	hpoissonpp[7]->SetMarkerStyle(kOpenDiamondCross);
	hpoissonpp[8]->SetMarkerStyle(kFullFourTrianglesPlus);
	hpoissonpp[9]->SetMarkerStyle(kOpenCrossX);

	hpoissonpp[0]->SetMarkerSize(2);
	hpoissonpp[1]->SetMarkerSize(2);
	hpoissonpp[2]->SetMarkerSize(2);
	hpoissonpp[3]->SetMarkerSize(2);
	hpoissonpp[4]->SetMarkerSize(2);
	hpoissonpp[5]->SetMarkerSize(2);
	hpoissonpp[6]->SetMarkerSize(2);
	hpoissonpp[7]->SetMarkerSize(2);
	hpoissonpp[8]->SetMarkerSize(2);
	hpoissonpp[9]->SetMarkerSize(2);
	
	for(Int_t j=0;j<9;j++){
		//hpoissonpp[j]=0;
	   hpoissonpp[j]->SetMarkerColor(j+1);
	   hpoissonpp[j]->SetLineColor(j+1);
           hpoissonpp[j]->Draw("SAMEP EX0");
        }
	hpoissonpp[9]->SetMarkerColor(49);
        hpoissonpp[9]->SetLineColor(49);
        hpoissonpp[9]->Draw("SAMEP EX0");
        fpoisspp->Draw("SAME");
	
	
	//Obtenemos las distribuciones dN/dpt pp que queremos simular y asignamos errores
	TH1 * dNdpt[10];
        TH1 * dNdpterrtotsys[10];
        TH1 * dNdpterrunsys[10];
	for(Int_t l=0 ; l < 10 ; l++)
		     {
			dNdpt[l] = 0;
        		dNdpterrtotsys[l] = 0;
        		dNdpterrunsys[l] = 0;
                        //Ref for SPD and V0M for V0M
        		fpp->GetObject(Form("Table 3/Hist1D_y%d",l+1),dNdpt[l]);        
        		fpp->GetObject(Form("Table 3/Hist1D_y%d_e1",l+1),dNdpterrtotsys[l]);
    			fpp->GetObject(Form("Table 3/Hist1D_y%d_e2",l+1),dNdpterrunsys[l]);
                        
			Int_t binesstapp = dNdpterrunsys[l]->GetNbinsX();
        
				
				/*
				//quitamos el factor 1/2pip_T, para obtener d2N/detadpT 
				Double_t twopi = 2*pi;
				for (Int_t bin=1 ; bin<=binesstapp ; bin++)
				{
					Double_t pTerrsys = dNdpterrtotsys[l]->GetBinCenter(bin);
					dNdpterrtotsys[l]->SetBinContent(bin,(dNdpterrtotsys[l]->GetBinError(bin))*twopi*pTerrsys);

				}

				for (Int_t bin=1 ; bin<=binesstapp ; bin++)
				{
					Double_t pTerrunsys = dNdpterrunsys[l]->GetBinCenter(bin);
					dNdpterrunsys[l]->SetBinContent(bin,(dNdpterrunsys[l]->GetBinError(bin))*twopi*pTerrunsys);

				}

				for (Int_t bin=1 ; bin<=binesstapp ; bin++)
				{
					Double_t pT = dNdpt[l]->GetBinCenter(bin);
					dNdpt[l]->SetBinContent(bin,(dNdpt[l]->GetBinContent(bin))*twopi*pT);
					//dNdpt[l]->SetBinError(bin,(dNdpt[l]->GetBinError(bin))*twopi*pT);
				}

				*/

				for (Int_t j = 0; j< binesstapp; j++)
           		           {
             				Double_t binerrorstapp = dNdpterrtotsys[l]->GetBinError(j+1);
             				Double_t binerrorsyspp = dNdpterrunsys[l]->GetBinError(j+1);
             				Double_t binerrtotpp = sqrt(pow(binerrorstapp,2)+pow(binerrorsyspp,2));
             				dNdpt[l]->SetBinError(j+1,binerrtotpp);
            
            			   }

				for (Int_t binppe=1; binppe<=binesstapp; binppe++){
						Double_t ptforppe = dNdpterrunsys[l]->GetBinCenter(binppe);
						if((ptforppe<0.15) || (ptforppe>10.) ){	

		                 		dNdpterrunsys[l]->SetBinError(binppe,0.);
						}
			   		}
		     }
	
	TCanvas * Qppb =0;
	Qppb = new TCanvas();
	Qppb->SetTitle("QpPb Spectra");
	Qppb->SetLogy(1);
	for (Int_t j = 0;j<9;j++){
		
		dNdpt[j]->SetMarkerStyle(kFullCircle);
                dNdpt[j]->SetMarkerColor(j+1);
	        dNdpt[j]->SetLineColor(j+1);
		dNdpt[j]->Draw("SAMEP EX0");
	}
	dNdpt[9]->SetMarkerStyle(kFullCircle);
	dNdpt[9]->SetMarkerColor(49);
	dNdpt[9]->SetLineColor(49);
        dNdpt[9]->Draw("SAMEP EX0");
	
	Double_t meanptdata[10];
	for (Int_t j = 0;j<10;j++){
		meanptdata[j]=0;
                meanptdata[j] = dNdpterrunsys[j]->GetMeanError();
	        cout << meanptdata[j] << endl;
	}
	
	
	
	
	//Quitamos delta pt y delta eta para obtener N_char de partículas
	
	Float_t sumpT;
	Float_t pt;
	Float_t averagepT;
	Int_t ranmult;
	Int_t countpt;
	Float_t sumdeltaN[10];
	for(Int_t k=0;k<10;k++){

		sumdeltaN[k]=0;
		Int_t ppbinessta = dNdpt[k]->GetNbinsX();
		
		for (Int_t binpt=1 ; binpt <= ppbinessta ; binpt++)
             		{
                 		Double_t deltapt = dNdpt[k]->GetBinWidth(binpt); 
                 		dNdpt[k]->SetBinContent(binpt,(dNdpt[k]->GetBinContent(binpt))*deltapt);
                 		dNdpt[k]->SetBinError(binpt,(dNdpt[k]->GetBinError(binpt))*deltapt); 
                        }
			
			dNdpt[k]->Scale(1.6);
		
		for (Int_t binpp=1; binpp<=ppbinessta; binpp++){
				Double_t ptforpp = dNdpt[k]->GetBinCenter(binpp);
				if((ptforpp<0.15) || (ptforpp>20.) ){	
				dNdpt[k]->SetBinContent(binpp,0.);
                 		dNdpt[k]->SetBinError(binpp,0.);
				}
			} 

			
			
			
		       //simulando de manera aleatoria las distribuciones dN/dp_t			
			dNdptrandom[k]=0;
			TAxis * x = dNdpt[k]->GetXaxis();
			const TArrayD * arx = x->GetXbins();
			dNdptrandom[k] = new TH1D(Form("dNdpt%d",k),Form("dNdpt%d",k),dNdpt[k]->GetNbinsX(),arx->GetArray());

			
			
			//llenamos de manera aleatoria
				
			 for (Int_t n = 0; n < Int_t(truepercentsV0M[k]*numberevents); n++)//loop over events
    	    			{
					
					ranmult = 0;
        				ranmult= hpoissonpp[k]->GetRandom();
					                 
        				dNdptrandom[k]->FillRandom(dNdpt[k],ranmult);
						
				        sumpT=0;
					pt=0;
					countpt=0;
					for (Int_t np = 0; np<ranmult; np++)//loop over particles
						{		
					             pt = dNdptrandom[k]->GetRandom();
						     if ((pt > 0.15) && (pt > 0)){
						                sumpT += pt;
								countpt++;
							 	}
						}
					averagepT = 0;
					averagepT = sumpT/(1.0*countpt);
					meanpt = averagepT;
					if(ranmult>0){	
					multiplicity = ranmult*1.0;
					}
				        cent=k;
					fTree->Fill();									
                                 }

			


	}

	


	//Comparando MC con los datos
 	TCanvas * MCvsdata1 = new TCanvas();
	MCvsdata1->SetTitle("MCvsdata1");
	MCvsdata1->SetLogy();
        auto legend1 = new TLegend(0.58,0.73,0.95,0.95);
	legend1->SetFillStyle(0);
	legend1->SetTextSize(0.03);
	legend1->SetHeader("V0M pp (#surd s_{NN} = 13 TeV)","C");
	legend1->AddEntry(dNdpt[0],"X" ,"p");
	legend1->AddEntry(dNdpt[1],"IX","p");
	legend1->AddEntry(dNdpt[2],"VIII" ,"p");
	legend1->AddEntry(dNdpt[3],"VII","p");
	legend1->AddEntry(dNdpt[4],"VI" ,"p");
	legend1->AddEntry(dNdpt[5],"V","p");
	legend1->AddEntry(dNdpt[6],"IV" ,"p");
	legend1->AddEntry(dNdpt[7],"III" ,"p");
	legend1->AddEntry(dNdpt[8],"II","p");
	legend1->AddEntry(dNdpt[9],"I" ,"p");

	auto legpoiss1 = new TLegend(0.58,0.50,0.95,0.70);
	legpoiss1->SetFillStyle(0);
	legpoiss1->SetTextSize(0.03);
	legpoiss1->SetHeader(" ","C");

	legpoiss1->AddEntry(dNdptrandom[0],"X Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[1],"IX Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[2],"VIII Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[3],"VII Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[4],"IV Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[5],"V Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[6],"IV Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[7],"III Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[8],"II Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[9],"I Toy MC" ,"p");
	//titulo delos ejes
	//dNdpt[1]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	//dNdpt[1]->GetYaxis()->SetTitle("N_{char}");
	for (Int_t i=0;i<9;i++){
		dNdpt[i]->SetMarkerStyle(kDiamond);
		dNdpt[i]->SetMarkerColor(i+1);
		dNdpt[i]->SetLineColor(i+1);

		dNdptrandom[i]->Scale(1.0/(truepercentsV0M[i]*numberevents));
		dNdptrandom[i]->SetMarkerStyle(kFullCircle);
		dNdptrandom[i]->SetMarkerColor(i+1);
		dNdptrandom[i]->SetLineColor(i+1);

		dNdptrandom[i]->Draw("SAME P");
		dNdpt[i]->Draw("SAME P ");
		
	}
	dNdpt[9]->SetMarkerStyle(kDiamond);
	dNdpt[9]->SetMarkerColor(49);
	dNdpt[9]->SetLineColor(49);
	dNdptrandom[9]->Scale(1.0/(truepercentsV0M[9]*numberevents));
	dNdptrandom[9]->SetMarkerStyle(kFullCircle);
	dNdptrandom[9]->SetMarkerColor(49);
	dNdptrandom[9]->SetLineColor(49);
	dNdptrandom[9]->Draw("SAME P ");
	dNdpt[9]->Draw("SAME P ");
	legend1->Draw("SAME");
	legpoiss1->Draw("SAME");
	
	//Ratio toyMC/data
	TCanvas *ratioMC = 0;
	ratioMC = new TCanvas();
	ratioMC->SetTitle("data/MC");
	TH1 * hr = 0;
	hr = ratioMC->DrawFrame(0.15,0.5,50.5,1.5);
	hr->SetXTitle("#it{p}_{T} (GeV/#it{c})");
	hr->SetYTitle("ToyMC/Data");
	hr->GetXaxis()->SetTitleSize(0.08);
	hr->GetYaxis()->SetTitleSize(0.08);
	hr->GetXaxis()->SetTitleOffset(.6);
	hr->GetYaxis()->SetTitleOffset(.6);
	auto rat = new TLegend(0.3,0.3,0.3,0.3);
	rat->SetFillStyle(0);
	rat->SetTextSize(0.04);
	rat->SetHeader("V0M pp (#surd s_{NN} = 13 TeV)","C");

        for (Int_t i=0;i<9;i++){
	hRat[i]=0;
	hRat[i] = (TH1 *)dNdptrandom[i]->Clone("hratio");
	hRat[i]->Divide(dNdpt[i]);
	hRat[i]->SetMarkerStyle(kFullCircle);
	hRat[i]->SetMarkerSize(2);
	hRat[i]->SetMarkerColor(i+1);
	hRat[i]->Draw("SAMEP EX0");
        }
	hRat[9] = (TH1 *)dNdptrandom[9]->Clone("hratio");
	hRat[9]->Divide(dNdpt[9]);
	hRat[9]->SetMarkerStyle(kFullCircle);
	hRat[9]->SetMarkerSize(2);
	hRat[9]->SetMarkerColor(49);
	hRat[9]->Draw("SAMEP EX0");



	rat->AddEntry(hRat[0],"X " ,"p");
	rat->AddEntry(hRat[1],"IX " ,"p");
	rat->AddEntry(hRat[2],"VIII " ,"p");
	rat->AddEntry(hRat[3],"VII " ,"p");
	rat->AddEntry(hRat[4],"VI " ,"p");
	rat->AddEntry(hRat[5],"V " ,"p");
	rat->AddEntry(hRat[6],"IV " ,"p");
	rat->AddEntry(hRat[7],"III " ,"p");
	rat->AddEntry(hRat[8],"II " ,"p");
	rat->AddEntry(hRat[9],"I" ,"p");
	rat->Draw("SAME");



	//Lo que reportamos es d²N/detadpt
	TCanvas * MCvsdata = 0;
	MCvsdata = new TCanvas();
	MCvsdata->SetTitle("FinalMCvsdata");
	MCvsdata->SetLogy(1);
	for(Int_t l =0;l<10;l++){
		hMC[l] = 0;
		hdata[l] = 0;
      		hMC[l] = (TH1 *)dNdptrandom[l]->Clone("hMC");
		hdata[l] = (TH1 *)dNdpt[l]->Clone("hdata");
		//agregamos deltaeta y deltapt
		Int_t ppbin = hdata[l]->GetNbinsX();
		for (Int_t ptbin = 1  ; ptbin <= ppbin ; ptbin++)
             		{
                 		Double_t detadpt = hdata[l]->GetBinWidth(ptbin); 
                 		hdata[l]->SetBinContent(ptbin,(hdata[l]->GetBinContent(ptbin))/(detadpt*1.6));
                 		hdata[l]->SetBinError(ptbin,(hdata[l]->GetBinError(ptbin))/(detadpt*1.6));

 		 		Double_t detadptMC = hMC[l]->GetBinWidth(ptbin); 
                 		hMC[l]->SetBinContent(ptbin,(hMC[l]->GetBinContent(ptbin))/(detadptMC*1.6));
                 		hMC[l]->SetBinError(ptbin,(hMC[l]->GetBinError(ptbin))/(detadptMC*1.6));
             		}
	}

	
        auto legend = new TLegend(0.60,0.65,0.90,0.90);
	legend->SetFillStyle(0);
	legend->SetTextSize(0.03);
	legend->SetHeader("V0M pp (#surd s_{NN} = 13 TeV)","C");
	legend->AddEntry(hdata[0],"X" ,"p");
	legend->AddEntry(hdata[1],"IX","p");
	legend->AddEntry(hdata[2],"VIII" ,"p");
	legend->AddEntry(hdata[3],"VII","p");
	legend->AddEntry(hdata[4],"VI" ,"p");
	legend->AddEntry(hdata[5],"V","p");
	legend->AddEntry(hdata[6],"IV" ,"p");
	legend->AddEntry(hdata[7],"III" ,"p");
	legend->AddEntry(hdata[8],"II","p");
	legend->AddEntry(hdata[9],"I" ,"p");

	auto legpoiss = new TLegend(0.58,0.50,0.95,0.70);
	legpoiss->SetFillStyle(0);
	legpoiss->SetTextSize(0.04);
	legpoiss->SetHeader(" ","C");

	legpoiss->AddEntry(hMC[0],"X Toy MC" ,"p");
	legpoiss->AddEntry(hMC[1],"IX Toy MC" ,"p");
	legpoiss->AddEntry(hMC[2],"VIII Toy MC" ,"p");
	legpoiss->AddEntry(hMC[3],"VII Toy MC" ,"p");
	legpoiss->AddEntry(hMC[4],"VI Toy MC" ,"p");
	legpoiss->AddEntry(hMC[5],"V Toy MC" ,"p");
	legpoiss->AddEntry(hMC[6],"IV Toy MC" ,"p");
	legpoiss->AddEntry(hMC[7],"III Toy MC" ,"p");
	legpoiss->AddEntry(hMC[8],"II Toy MC" ,"p");
	legpoiss->AddEntry(hMC[9],"I Toy MC" ,"p");
	//titulo de los ejes
	hMC[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	hMC[0]->GetYaxis()->SetTitle("d^{2}N/d#etad#it{p}_{T} (GeV/#it{c})^{-1}");
	hMC[0]->GetXaxis()->SetTitleSize(0.08);
	hMC[0]->GetYaxis()->SetTitleSize(0.08);
	hMC[0]->GetXaxis()->SetTitleOffset(.6);
	hMC[0]->GetYaxis()->SetTitleOffset(.6);
	for (Int_t i=0;i<9;i++){
		hMC[i]->SetMarkerStyle(kDiamond);
		hMC[i]->SetMarkerColor(i+1);
		hMC[i]->SetLineColor(i+1);

		hMC[i]->SetMarkerSize(2);
		hdata[i]->SetMarkerStyle(kFullCircle);
		hdata[i]->SetMarkerColor(i+1);
		hdata[i]->SetLineColor(i+1);
		hdata[i]->SetMarkerSize(2);

		hMC[i]->Draw("SAME P");
		hdata[i]->Draw("SAME P");
	}
	hMC[9]->SetMarkerStyle(kDiamond);
	hMC[9]->SetMarkerColor(49);
	hMC[9]->SetLineColor(49);
	hMC[9]->SetMarkerSize(2);
	hdata[9]->SetMarkerStyle(kFullCircle);
	hdata[9]->SetMarkerColor(49);
	hdata[9]->SetLineColor(49);
	hdata[9]->SetMarkerSize(2);

	hMC[9]->Draw("SAME P");
	hdata[9]->Draw("SAME P");

	legend->Draw("SAME");
	legpoiss->Draw("SAME");
	
	
	//Guardamos el árbol
	outFiletree->cd();
	fTree->Write();
	outFiletree->Close();
	
	




}


void MCforpp502(const Char_t * inFilepp, const Char_t * outDirpp){

	cout<<" Creating directory for dNdptppMc..."<<endl;
	CreateDir(outDirpp); 

	cout<<" Opening root file..."<<endl;
	//Archivo de salida
	TFile * outFiletree = 0;
	outFiletree = new TFile(Form("%s/newcutTreeV0Mpp502.root",outDirpp),"RECREATE");

	//<dN/deta>*deltaeta for SPD and V0M at 5.02 Tev
	Double_t meandNSPD[9]={55.36,47.84,41.92,35.84,29.60,23.36,16.96,10.53,3.54};
	Double_t truemeandNSPD[9]={3.54,10.53,15.96,23.36,29.60,35.84,41.92,47.84,55.36};
	Double_t meandNV0M[10]={30.72,24.16,19.84,17.12,15.15,12.86,10.50,8.62,6.48,3.63};
	Double_t truemeandNV0M[10]={3.63,6.48,8.62,10.50,12.86,15.15,17.12,19.84,24.16,30.71};
	//Numero de entradas para las dist de poisson
        Double_t NentriesSPD[9]={5,16,44,112,270,605,1267,2571,50000};
	Double_t NentriesV0M[10]={166,666,833,833,833,1666,1666,1666,3333,5000};
	//Porcentajes eventos
	Double_t percentsSPD[9] = {0.00079,0.00165,0.00447,0.0114,0.02733,0.0612,0.1281,0.2598,0.5052};
	Double_t truepercentsSPD[9]={0.5052,0.2598,0.1281,0.0612,0.02733,0.0114,0.00447,0.00165,0.00079}; 
        Double_t percentsV0M[10] = {0.01,0.04,0.05,0.05,0.05,0.1,0.1,0.1,0.2,0.3};
	Double_t truepercentsV0M[10]={0.3,0.2,0.1,0.1,0.1,0.05,0.05,0.05,0.04,0.01};
	Int_t numberevents = 1700000;//numero de eventos 
        gRandom->SetSeed(0);
	
	//Creamos el árbol
	fTree = 0;
	fTree = new TTree("tree","Event data");
	fTree->Branch("meanpt",&meanpt,"meanpt/F");
	fTree->Branch("multiplicity",&multiplicity,"multiplicity/F");
	fTree->Branch("cent",&cent,"cent/F");
	//archivoentrada
	TFile * fpp = 0;
	fpp = TFile::Open(inFilepp);

	//Iniciamos histogramas
	TH1 * dNdptrandom[10];
	TH1 * hpoissonpp[10];
	TF1 * fhpp[10];
	TH1 * yrandompp[10];
	TH1 * hRatpp[10];
	TH1 * hMC[10];
	TH1 * hdata[10];

	//Simulamos <dN/deta>_pp como distribuciones de Poisson
	for (Int_t s=0 ; s<10; s++)
    	   {  
       
      		hpoissonpp[s] = 0;
      		fhpp[s] = 0;
          	hpoissonpp[s] = new TH1F(Form("h%d",s),"Poisson",100,0,100);
          	fhpp[s] = new TF1(Form("fh%d",s),"TMath::Poisson(x,[0])",0,100);
          	fhpp[s]->SetParameter(0,truemeandNV0M[s]);
          	hpoissonpp[s]->Sumw2();
          	hpoissonpp[s]->FillRandom(Form("fh%d",s),truepercentsV0M[s]*numberevents);
          }

	// Load necessary libraries
	LoadLibs();
	// Set the default style
	SetStyle();
        //Graficando las dist de poisson para pp
	TCanvas * poissonpp = new TCanvas();
	poissonpp->SetTitle("Poisson");
	auto fpoisspp = new TLegend(0.58,0.73,0.95,0.95);
	fpoisspp->SetFillStyle(0);
	fpoisspp->SetTextSize(0.03);
	fpoisspp->SetHeader("V0M","C");
	fpoisspp->AddEntry(hpoissonpp[0],"I" ,"p");
	fpoisspp->AddEntry(hpoissonpp[1],"II ","p");
	fpoisspp->AddEntry(hpoissonpp[2],"III " ,"p");
	fpoisspp->AddEntry(hpoissonpp[3],"IV ","p");
	fpoisspp->AddEntry(hpoissonpp[4],"V" ,"p");
	fpoisspp->AddEntry(hpoissonpp[5],"VI","p");
	fpoisspp->AddEntry(hpoissonpp[6],"VII" ,"p");
	fpoisspp->AddEntry(hpoissonpp[7],"VIII" ,"p");
	fpoisspp->AddEntry(hpoissonpp[8],"IX" ,"p");
	fpoisspp->AddEntry(hpoissonpp[9],"X" ,"p");
	
	//Titulos a los ejes para las distribuciones de Poisson	
	hpoissonpp[0]->GetXaxis()->SetTitle("N_{ch}");
	hpoissonpp[0]->GetYaxis()->SetTitle("P(N_{ch}) ");
	hpoissonpp[0]->SetAxisRange(0, 800, "Y");
	//Estilo de las distribuciones
	hpoissonpp[0]->SetMarkerStyle(kFullDiamond);
	hpoissonpp[1]->SetMarkerStyle(kFullDoubleDiamond);
	hpoissonpp[2]->SetMarkerStyle(kFullCrossX);
	hpoissonpp[3]->SetMarkerStyle(kFullCircle);
	hpoissonpp[4]->SetMarkerStyle(kFullStar);
	hpoissonpp[5]->SetMarkerStyle(kFullSquare);
	hpoissonpp[6]->SetMarkerStyle(kFullTriangleUp);
	hpoissonpp[7]->SetMarkerStyle(kOpenDiamondCross);
	hpoissonpp[8]->SetMarkerStyle(kFullFourTrianglesPlus);
	hpoissonpp[9]->SetMarkerStyle(kOpenCrossX);
	
	for(Int_t j=0;j<9;j++){
	   hpoissonpp[j]->SetMarkerColor(j+1);
	   hpoissonpp[j]->SetLineColor(j+1);
           hpoissonpp[j]->Draw("SAMEP EX0");
        }
	hpoissonpp[9]->SetMarkerColor(49);
        hpoissonpp[9]->SetLineColor(49);
        hpoissonpp[9]->Draw("SAMEP EX0");
        fpoisspp->Draw("SAME");
        fpoisspp->Draw("SAME");
	
	//Obtenemos las distribuciones dN/dpt pp que queremos simular y asignamos errores
	TH1 * dNdpt[10];
        TH1 * dNdpterrtotsys[10];
        TH1 * dNdpterrunsys[10];
	for(Int_t l=0 ; l < 10 ; l++)
		     {
			dNdpt[l] = 0;
        		dNdpterrtotsys[l] = 0;
        		dNdpterrunsys[l] = 0;
                        //Ref for SPD and V0M for V0M
        		fpp->GetObject(Form("Table 4/Hist1D_y%d",l+1),dNdpt[l]);        
        		fpp->GetObject(Form("Table 4/Hist1D_y%d_e1",l+1),dNdpterrtotsys[l]);
    			fpp->GetObject(Form("Table 4/Hist1D_y%d_e2",l+1),dNdpterrunsys[l]);
                        
			Int_t binesstapp = dNdpterrunsys[l]->GetNbinsX();
        			/*
				//quitamos el factor 1/2pip_T, para obtener d2N/detadpT 
				Double_t twopi = 2*pi;
				for (Int_t bin=1 ; bin<=binesstapp ; bin++)
				{
					Double_t pTerrsys = dNdpterrtotsys[l]->GetBinCenter(bin);
					dNdpterrtotsys[l]->SetBinContent(bin,(dNdpterrtotsys[l]->GetBinError(bin))*twopi*pTerrsys);

				}

				for (Int_t bin=1 ; bin<=binesstapp ; bin++)
				{
					Double_t pTerrunsys = dNdpterrunsys[l]->GetBinCenter(bin);
					dNdpterrunsys[l]->SetBinContent(bin,(dNdpterrunsys[l]->GetBinError(bin))*twopi*pTerrunsys);

				}

				for (Int_t bin=1 ; bin<=binesstapp ; bin++)
				{
					Double_t pT = dNdpt[l]->GetBinCenter(bin);
					dNdpt[l]->SetBinContent(bin,(dNdpt[l]->GetBinContent(bin))*twopi*pT);
					dNdpt[l]->SetBinError(bin,(dNdpt[l]->GetBinError(bin))*twopi*pT);
				}

			
				*/

				for (Int_t j = 0; j< binesstapp; j++)
           		           {
             				Double_t binerrorstapp = dNdpterrtotsys[l]->GetBinError(j+1);
             				Double_t binerrorsyspp = dNdpterrunsys[l]->GetBinError(j+1);
             				Double_t binerrtotpp = sqrt(pow(binerrorstapp,2)+pow(binerrorsyspp,2));
             				dNdpt[l]->SetBinError(j+1,binerrtotpp);
            
            			   }

					for (Int_t binppe=1; binppe<=binesstapp; binppe++){
						Double_t ptforppe = dNdpterrunsys[l]->GetBinCenter(binppe);
						if((ptforppe<0.15) || (ptforppe>10.) ){	

		                 		dNdpterrunsys[l]->SetBinError(binppe,0.);
						}
			   		}
		     
		     }

	TCanvas * Qppb =0;
	Qppb = new TCanvas();
	Qppb->SetTitle("QpPb Spectra");
	Qppb->SetLogy(1);
	for (Int_t j = 0;j<9;j++){
		
		dNdpt[j]->SetMarkerStyle(kFullCircle);
                dNdpt[j]->SetMarkerColor(j+1);
	        dNdpt[j]->SetLineColor(j+1);
		dNdpt[j]->Draw("SAMEP EX0");
	}
	dNdpt[9]->SetMarkerStyle(kFullCircle);
	dNdpt[9]->SetMarkerColor(49);
	dNdpt[9]->SetLineColor(49);
        dNdpt[9]->Draw("SAMEP EX0");

	Double_t meanptdataerror[10];
	for (Int_t j = 0;j<10;j++){
		meanptdataerror[j]=0;
                meanptdataerror[j] = dNdpterrunsys[j]->GetMeanError();
	        cout << meanptdataerror[j] << endl;
	}

	//Quitamos delta pt para obtener N_char de partículas
	
	Float_t sumpT;
	Float_t pt;
	Float_t averagepT;
	Int_t ranmult;
	Int_t countpt;
	for(Int_t k=0;k<10;k++){
		Int_t ppbinessta = dNdpt[k]->GetNbinsX();
		for (Int_t binpt=1 ; binpt <= ppbinessta ; binpt++)
             		{
                 		Double_t deltapt = dNdpt[k]->GetBinWidth(binpt); 
                 		dNdpt[k]->SetBinContent(binpt,(dNdpt[k]->GetBinContent(binpt))*deltapt);
                 		dNdpt[k]->SetBinError(binpt,(dNdpt[k]->GetBinError(binpt))*deltapt); 
                        }
			dNdpt[k]->Scale(1.6);
		       //simulando de manera aleatoria las distribuciones dN/dp_t
			TAxis * x = dNdpt[k]->GetXaxis();
			const TArrayD * arx = x->GetXbins();
			dNdptrandom[k] = new TH1D(Form("dNdpt%d",k),Form("dNdpt%d",k),dNdpt[k]->GetNbinsX(),arx->GetArray());

			for (Int_t binpp=1; binpp<=ppbinessta; binpp++){
				Double_t ptforpp = dNdpt[k]->GetBinCenter(binpp);
				if((ptforpp<0.15) || (ptforpp>20.) ){	
				dNdpt[k]->SetBinContent(binpp,0.);
                 		dNdpt[k]->SetBinError(binpp,0.);
				}
			} 
			
			//llenamos de manera aleatoria
				
			 for (Int_t n = 0; n < Int_t(truepercentsV0M[k]*numberevents); n++)//loop over events
    	    			{
					
					ranmult = 0;
        				ranmult= hpoissonpp[k]->GetRandom();
					                 
        				dNdptrandom[k]->FillRandom(dNdpt[k],ranmult);
						
				        sumpT=0;
					pt=0;
					countpt=0;
					for (Int_t np = 0; np<ranmult; np++)//loop over particles
						{		
					             pt = dNdptrandom[k]->GetRandom();
						     if ((pt > 0.15) && (pt > 0)){
						                sumpT += pt;
								countpt++;
							 	}
						}
					averagepT = 0;
					averagepT = sumpT/(1.0*countpt);
					meanpt = averagepT;
					if(ranmult>0){	
					multiplicity = ranmult*1.0;
					}
				        cent=k;
					fTree->Fill();									
                                 }
					


	}
	
	

	//Comparando MC con los datos
 	TCanvas * MCvsdata1 = new TCanvas();
	MCvsdata1->SetTitle("MCvsdata1");
	MCvsdata1->SetLogy();
        auto legend1 = new TLegend(0.58,0.73,0.95,0.95);
	legend1->SetFillStyle(0);
	legend1->SetTextSize(0.03);
	legend1->SetHeader("V0M (#surd s_{NN} = 5.02 TeV)","C");
	legend1->AddEntry(dNdpt[0],"I" ,"p");
	legend1->AddEntry(dNdpt[1],"II","p");
	legend1->AddEntry(dNdpt[2],"III" ,"p");
	legend1->AddEntry(dNdpt[3],"IV","p");
	legend1->AddEntry(dNdpt[4],"V" ,"p");
	legend1->AddEntry(dNdpt[5],"VI","p");
	legend1->AddEntry(dNdpt[6],"VII" ,"p");
	legend1->AddEntry(dNdpt[7],"VIII" ,"p");
	legend1->AddEntry(dNdpt[8],"IX","p");
	legend1->AddEntry(dNdpt[9],"X" ,"p");

	auto legpoiss1 = new TLegend(0.58,0.50,0.95,0.70);
	legpoiss1->SetFillStyle(0);
	legpoiss1->SetTextSize(0.03);
	legpoiss1->SetHeader(" ","C");

	legpoiss1->AddEntry(dNdptrandom[0],"I Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[1],"II Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[2],"III Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[3],"IV Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[4],"V Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[5],"VI Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[6],"VII Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[7],"VIII Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[8],"IX Toy MC" ,"p");
	legpoiss1->AddEntry(dNdptrandom[9],"X Toy MC" ,"p");
	//titulo delos ejes
	//dNdpt[1]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	//dNdpt[1]->GetYaxis()->SetTitle("N_{char}");
	for (Int_t i=0;i<9;i++){
		dNdpt[i]->SetMarkerStyle(kDiamond);
		dNdpt[i]->SetMarkerColor(i+1);
		dNdpt[i]->SetLineColor(i+1);

		dNdptrandom[i]->Scale(1.0/(truepercentsV0M[i]*numberevents));
		dNdptrandom[i]->SetMarkerStyle(kFullCircle);
		dNdptrandom[i]->SetMarkerColor(i+1);
		dNdptrandom[i]->SetLineColor(i+1);

		dNdptrandom[i]->Draw("SAME P");
		dNdpt[i]->Draw("SAME P ");
		
	}
	dNdpt[9]->SetMarkerStyle(kDiamond);
	dNdpt[9]->SetMarkerColor(49);
	dNdpt[9]->SetLineColor(49);
	dNdptrandom[9]->Scale(1.0/(truepercentsV0M[8]*numberevents));
	dNdptrandom[9]->SetMarkerStyle(kFullCircle);
	dNdptrandom[9]->SetMarkerColor(49);
	dNdptrandom[9]->SetLineColor(49);
	dNdptrandom[9]->Draw("SAME P ");
	dNdpt[9]->Draw("SAME P ");
	legend1->Draw("SAME");
	legpoiss1->Draw("SAME");
	
	//Lo que reportamos es d²N/detadpt
	TCanvas * MCvsdata = 0;
	MCvsdata = new TCanvas();
	MCvsdata->SetTitle("FinalMCvsdata");
	MCvsdata->SetLogy(1);
	for(Int_t l =0;l<10;l++){
		hMC[l] = 0;
		hdata[l] = 0;
      		hMC[l] = (TH1 *)dNdptrandom[l]->Clone("hMC");
		hdata[l] = (TH1 *)dNdpt[l]->Clone("hdata");
		//agregamos deltaeta y deltapt
		Int_t ppbin = hdata[l]->GetNbinsX();
		for (Int_t ptbin = 1  ; ptbin <= ppbin ; ptbin++)
             		{
                 		Double_t detadpt = hdata[l]->GetBinWidth(ptbin); 
                 		hdata[l]->SetBinContent(ptbin,(hdata[l]->GetBinContent(ptbin))/(detadpt*1.6));
                 		hdata[l]->SetBinError(ptbin,(hdata[l]->GetBinError(ptbin))/(detadpt*1.6));

 		 		Double_t detadptMC = hMC[l]->GetBinWidth(ptbin); 
                 		hMC[l]->SetBinContent(ptbin,(hMC[l]->GetBinContent(ptbin))/(detadptMC*1.6));
                 		hMC[l]->SetBinError(ptbin,(hMC[l]->GetBinError(ptbin))/(detadptMC*1.6));
             		}
	}


 	
        auto legend = new TLegend(0.60,0.65,0.90,0.90);
	legend->SetFillStyle(0);
	legend->SetTextSize(0.03);
	legend->SetHeader("V0M (#surd s_{NN} = 5.02 TeV)","C");
	legend->AddEntry(hdata[0],"I" ,"p");
	legend->AddEntry(hdata[1],"II","p");
	legend->AddEntry(hdata[2],"III" ,"p");
	legend->AddEntry(hdata[3],"IV","p");
	legend->AddEntry(hdata[4],"V" ,"p");
	legend->AddEntry(hdata[5],"VI","p");
	legend->AddEntry(hdata[6],"VII" ,"p");
	legend->AddEntry(hdata[7],"VIII" ,"p");
	legend->AddEntry(hdata[8],"IX","p");
	legend->AddEntry(hdata[9],"X" ,"p");

	auto legpoiss = new TLegend(0.58,0.50,0.95,0.70);
	legpoiss->SetFillStyle(0);
	legpoiss->SetTextSize(0.03);
	legpoiss->SetHeader(" ","C");

	legpoiss->AddEntry(hMC[0],"I Toy MC" ,"p");
	legpoiss->AddEntry(hMC[1],"II Toy MC" ,"p");
	legpoiss->AddEntry(hMC[2],"III Toy MC" ,"p");
	legpoiss->AddEntry(hMC[3],"IV Toy MC" ,"p");
	legpoiss->AddEntry(hMC[4],"V Toy MC" ,"p");
	legpoiss->AddEntry(hMC[5],"VI Toy MC" ,"p");
	legpoiss->AddEntry(hMC[6],"VII Toy MC" ,"p");
	legpoiss->AddEntry(hMC[7],"VIII Toy MC" ,"p");
	legpoiss->AddEntry(hMC[8],"IX Toy MC" ,"p");
	legpoiss->AddEntry(hMC[9],"X Toy MC" ,"p");
	//titulo de los ejes
	hMC[0]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	hMC[0]->GetYaxis()->SetTitle("d^{2}N/d#etad#it{p}_{T} (GeV/#it{c})^{-1}");
	for (Int_t i=0;i<9;i++){
		hMC[i]->SetMarkerStyle(kDiamond);
		hMC[i]->SetMarkerColor(i+1);
		hMC[i]->SetLineColor(i+1);

		
		hdata[i]->SetMarkerStyle(kFullCircle);
		hdata[i]->SetMarkerColor(i+1);
		hdata[i]->SetLineColor(i+1);

		hMC[i]->Draw("SAME P");
		hdata[i]->Draw("SAME P");
	}
	hMC[9]->SetMarkerStyle(kDiamond);
	hMC[9]->SetMarkerColor(49);
	hMC[9]->SetLineColor(49);
	
	hdata[9]->SetMarkerStyle(kFullCircle);
	hdata[9]->SetMarkerColor(49);
	hdata[9]->SetLineColor(49);

	hMC[9]->Draw("SAME P");
	hdata[9]->Draw("SAME P");

	legend->Draw("SAME");
	legpoiss->Draw("SAME");

	/*
	//Guardamos el árbol
	outFiletree->cd();
	fTree->Write();
	outFiletree->Close();
	*/
	




}




void SetStyle(Bool_t graypalette) {
	cout << "Setting style!" << endl;

	gStyle->Reset("Plain");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	if(graypalette) gStyle->SetPalette(8,0);
	else gStyle->SetPalette(1);
	gStyle->SetCanvasColor(10);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetFrameLineWidth(1);
	gStyle->SetFrameFillColor(kWhite);
	gStyle->SetPadColor(10);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadLeftMargin(0.12);

	gStyle->SetPadTopMargin(0.03);
	gStyle->SetPadRightMargin(0.03);

	gStyle->SetHistLineWidth(1);
	gStyle->SetHistLineColor(kRed);
	gStyle->SetFuncWidth(2);
	gStyle->SetFuncColor(kGreen);
	gStyle->SetLineWidth(2);
	gStyle->SetLabelSize(0.045,"xyz");
	gStyle->SetLabelOffset(0.002,"y");
	gStyle->SetLabelOffset(0.001,"x");
	gStyle->SetLabelColor(kBlack,"xyz");
	gStyle->SetTitleSize(0.05,"xyz");
	gStyle->SetTitleOffset(1.12,"y");
	gStyle->SetTitleOffset(0.95,"x");
	gStyle->SetTitleFillColor(kWhite);
	gStyle->SetTextSizePixels(26);
	gStyle->SetTextFont(42);

	gStyle->SetLegendBorderSize(0);
	gStyle->SetLegendFillColor(kWhite);
	//  gStyle->SetFillColor(kWhite);
	gStyle->SetLegendFont(42);
}

void LoadLibs() {
	gSystem->Load("libCore.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libVMC");
	gSystem->Load("libTree");
	gSystem->Load("libMinuit");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libCORRFW");
	gSystem->Load("libPWGTools");
}
