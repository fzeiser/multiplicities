//////////////////////////////////////////////////////
//  Script to plot fission gamma spectra and multip.//
//  Fabio, Jan 2016                                 //
//////////////////////////////////////////////////////
//#include <string>
{
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetOptTitle(0);
	// gStyle.SetOptStat(0);
	// gStyle.SetFillColor(0);
	gStyle.SetPadBorderMode(0);
	gStyle->SetOptStat("nemri");
    
	m = (TH1F*)gROOT->FindObject("h");
	if (m) m->Delete();
	TCanvas *c1_1 = new TCanvas("c1","Fission Gamma Multiplicities");
	// c1->Divide(1,1,0,0);

    // Threshold Gamma energy (currently for the integral & Normalization
    // and Multiplicities! )
    int Egammamin = 450;    //diskusjon med Cec: i alle fall 400?! 300 er "original" threshold
    int Egammamax = 12e3;
    
        
    //================== ====== ==================
    // Neutron/gamma correction factor (~20% of all counts are neutrons. Correct for this. Ask Sunniva or Fabio)
    // float neutronCorrectionFactor = 0.78;
    float neutronCorrectionFactor = 1.0;   //max 
    // float neutronCorrectionFactor = 0.84;   //min 
    float oneMinusNeutronCorrectionFactor = 1.0 - neutronCorrectionFactor;
    // Based on Häusser NIM 1983, we extracted a exponential fitting function:
    // f(E_g) = A * exp(-2*E_g)
    // The total area is (assuming upper limit = infinity), a being the lower limit
    // \int_a^\infty f(E_g) dE_g = 0.5*A*exp(-2.0*a)
    // We solve for A so that total area is oneMinusNeutronCorrectionFactor of total count which gives:
    // float neutronCorrectionFunctionScaling = 2.0 * oneMinusNeutronCorrectionFactor * 1.072 * 1e6 / exp(-2.0*Egammamin/1000.); // divided by 1000 because of units
    float neutronCorrectionFunctionScaling = 2.0 * oneMinusNeutronCorrectionFactor * 1.072 * 1e6; // We numerically find the scaling further down
    
    

    // How many projection shall be made
	int nIntervals = 7;
	// Region of Interest
	int Emin = 4.5e3; // Excitation energy
    int Emax = 8.8e3;
	int Erange =  Emax - Emin; 

    // Efficiency of the NaI detectors
    // denne avhenger av antallet detektorer
    // I Sunnivas setup: 25(?) detektorer
    double efficiency = 0.141;   //originalt fra Fabio - altså for 26 detektorer. For 25 stykk (Sunnivas): 0.136
    // oups - read the comment!
    // !!!!! set to 0 --> There is an error, but it's systematic! and we ignore it!
    double efficiencyErr = 0.01; 
    if (efficiencyErr!=0){
        cout << "WARNING: Review the error propagation on this part! - maybe somewhere: assuming this is 0" << endl;
    }
    


    float nGammaCoinc[nIntervals],nGammaCoincErr[nIntervals],
          nFissions[nIntervals],  nFissionsErr[nIntervals],
          multiplicity[nIntervals],  multiplicityErr[nIntervals],  multiplicityLow[nIntervals],  multiplicityUp[nIntervals],
          AverageEnergy[nIntervals], AverageEnergyErr[nIntervals], AverageEnergyLow[nIntervals], AverageEnergyUp[nIntervals],
          TotalEnergy[nIntervals],   TotalEnergyErr[nIntervals],   TotalEnergyLow[nIntervals],   TotalEnergyUp[nIntervals],
          meanEnergies[nIntervals];
          


    //================== Reading a Mama Matrix ==================
    //================== ======                ==================


    m = (TH2F*)gROOT->FindObject("matrix");
    if (m) { 
        m->Delete(); 
    }
    // m2 = (TH2F*)gROOT->FindObject("matrix2");
    // if (m2) {m2->Delete();}                    

    // declarations and stuff
    char filename[10000];
    sprintf(filename, "../m_alfna_fiss_un");
    // sprintf(filename, "../fiss_n_corrected_mama_un");
    ifstream  ifile1(filename);
    string line;
    string cal_dummy;
    string dim_dummy;
    char pdf_filename[512];
    int dim;
    int dim_start;
    int dim_stop;
    int dim_size;
    int position;
    int file_length;
    line.resize(200);   // need long enough line to read MAMA headers
    double x_cal[3] = {0.,1.,0.};   // calibration coeffs. on x axis: a0, a1, a2
    double y_cal[3] = {0.,1.,0.};   // calibration coeffs. on y axis: a0, a1, a2
    int dx, dy; // dimension on x and y axis
    int ix, iy;
    double value;
    double x,y , d1, d2;
    double number_of_counts = 0.;
    double new_y1, new_y2; 
    int sign_ycal;

    
    // open file to read
    if(!ifile1){
        cout << "\n Could not open file " << filename << endl;
        exit(1);
    } else {
        cout << "\n Successful opening of file"  << endl;
    }

    // read MAMA header (fixed format). The 10 first lines are info text 
    if(!getline(ifile1,line) || line.substr(0,10) != "!FILE=Disk"){ // check correct format
        printf("\n This is not a MAMA file!!!\n ");
        exit(2);
    }   
    getline(ifile1,line);   // skip !KIND=Spectrum
    getline(ifile1,line);   // skip !LABORATORY=Oslo Cyclotron Laboratory (OCL) 
    getline(ifile1,line);   // skip !EXPERIMENT=mama
    getline(ifile1,line);   // skip !COMMENT=Sorted simulated data
    getline(ifile1,line);   // skip !TIME=DATE:    19/11/09 11:47:26
    getline(ifile1,line);   // get line with calibration
    cout << "\n Reading calibration coeffs.:" << endl;
    // calibration on x axis
    cal_dummy = line.substr(20,13); // position 20, length 13 characters
    if(!(istringstream(cal_dummy) >> x_cal[0])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a0 on x axis is: " << x_cal[0] << " keV." << endl;
    cal_dummy = line.substr(34,13); 
    if(!(istringstream(cal_dummy) >> x_cal[1])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a1 on x axis is: " << x_cal[1] << " keV/ch." << endl;
    cal_dummy = line.substr(48,13); 
    if(!(istringstream(cal_dummy) >> x_cal[2])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a2 on x axis is: " << x_cal[2] << " (keV/ch)^2." << endl;
    // calibration on y axis
    cal_dummy = line.substr(62,13); 
    if(!(istringstream(cal_dummy) >> y_cal[0])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a0 on y axis is: " << y_cal[0] << " keV." << endl;
    cal_dummy = line.substr(76,13); 
    if(!(istringstream(cal_dummy) >> y_cal[1])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a1 on y axis is: " << y_cal[1] << " keV/ch." << endl;
    cal_dummy = line.substr(90,13); 
    if(!(istringstream(cal_dummy) >> y_cal[2])) cout << "Could not convert string to number." << endl;
    else cout << " Calibration coeff. a2 on y axis is: " << y_cal[2] << " (keV/ch)^2.\n" << endl;
    getline(ifile1,line);   // skip !PRECISION=16
    getline(ifile1,line);   // get dimension
    // dimension of matrix
    dim_start = line.find_first_of("=") + 1;
    dim_dummy = line.substr(dim_start,1);
    if(!(istringstream(dim_dummy) >> dim)) cout << "Could not convert string to number." << endl;
    else cout << " Dimension of matrix is: " << dim << endl;    
    getline(ifile1,line);   // get channels
    // dimension on x axis
    dim_start = line.find_first_of(":") + 1;
    dim_stop = line.find_last_of(",");
    dim_size = dim_stop - dim_start;
    dim_dummy = line.substr(dim_start,dim_size);
    if(!(istringstream(dim_dummy) >> dx)) cout << "Could not convert string to number." << endl;
    else cout << " Dimension on x axis is: " << dx+1 << " ch." << endl; 
    dx = dx+1;
    // dimension on y axis
    dim_start = line.find_last_of(":");
    dim_stop = line.find_last_of(")");
    dim_size = dim_stop - dim_start;
    dim_dummy = line.substr(dim_start+1,dim_size-1);
    if(!(istringstream(dim_dummy) >> dy)) cout << "Could not convert string to number." << endl;
    else cout << " Dimension on y axis is: " << dy+1 << " ch." << endl; 
    dy = dy+1;

    // Test if negative calibration coeff. on Ex, then invert axis:
    if(y_cal[1] < 0.){
        sign_ycal = -1;
        new_y1 = y_cal[0] + (y_cal[1]*(double) dy);
        new_y2 = y_cal[0];
        
        y_cal[0] = new_y1;
        y_cal[1] = (-1.)*y_cal[1];
        cout << " New calibration on y axis: y_cal[0] = " << y_cal[0] << ", y_cal[1] = " << y_cal[1] << endl;

    }
        
    // x_cal[0] /= 1000.; // to get MeV instead of keV
    // x_cal[1] /= 1000.;
    // y_cal[0] /= 1000.; // to get MeV instead of keV
    // y_cal[1] /= 1000.;

    
    // Make histogram 
    cout << "Creating ALFNA matrix with number of bins: " << dx << ", " << dy << endl;
    cout << "Histogram lowest bin center x: " << x_cal[0] << " and highest " << dx*x_cal[1]+x_cal[0] << endl;
    cout << "Histogram lowest bin center y: " << y_cal[0] << " and highest " << dy*y_cal[1]+y_cal[0] << endl;
    TH2D *matrix = new TH2D("matrix"," ",dx,x_cal[0],dx*x_cal[1]+x_cal[0],dy,y_cal[0],dy*y_cal[1]+y_cal[0]);
    matrix->SetOption("colz");
    gStyle->SetPalette(1);
    

    if(sign_ycal < 0.){ // if negative calibration coeff. on y axis
        for(iy=dy;iy>0;iy--){
            for(ix=0;ix<dx;ix++){
                ifile1 >> value;
                number_of_counts += value;
                matrix->SetBinContent(ix,iy,value);
            }
        }
    } else{   // if positive calibration coeff. on y axis
        for(iy=0;iy<dy;iy++){
            for(ix=0;ix<dx;ix++){
                ifile1 >> value;
                number_of_counts += value;
                matrix->SetBinContent(ix,iy,value);
            }            
        }
    }



///////////////////////////
// READ THE DATA : Verbinski
///////////////////////////
// double x,y;
// string string_x;
// char delim;

// opening the data_file
// char *myfile = "Verbinski_239Pu.csv".c_str();

ifstream thefile;
thefile.open("Verbinski_239Pu.csv");
int nheaderlines = 1;

// throw error if file does not exists
if(thefile.fail()){
  cout << "\n Could not open file!!!: "  << "\n ";
  exit(1);
}
else cout << "\n Successful opening of file" << endl;


int nLines=0;
//std::ifstream myfile("main.cpp");
std::string lines;
while (std::getline(thefile, lines))
       ++nLines;
thefile.close();

thefile.open("Verbinski_239Pu.csv");
string line; // ignore the header line
line.resize(500); //

for(int i=0;i<nheaderlines;i++){ // (for(i=0;i<2;i++) --> 3 header lines)
getline(thefile,line); // 
} // ignore the header line

nPoints = nLines - nheaderlines;

// shift_filenames.resize(nPoints);
// shift_channels.resize(nPoints);
double energy[nPoints], spectrum[nPoints], lower_y[nPoints], higher_y[nPoints];

  //int i = 0;
  double d1_1, d2_2;

double percentage_error = 0.105; // to be assumed if we didn't retreive the error
  while(thefile){
    thefile >> x >> y >> d1_1 >> d2_2;            // readfrist columns into "x"
    if(i<nPoints){
            energy[i]   = x*1000; // to get to keV
            spectrum[i] = y;
            d1 =d1_1;
            d2 =d2_2;
            cout << d1 << d2 << endl; 
            if (d1>=0){
            higher_y[i] = d1 - y;
            lower_y[i]  = y - d2;
            }
            else{
            higher_y[i] = percentage_error * y;
            lower_y[i]  = percentage_error * y;   
            }

            // cout << energy[i] << "\t" << spectrum[i] << endl;
    }
    i++;
  }

  // for(int i=0;i<nPoints;i++)
  // {
  //   cout << energy[i] << "\t" << spectrum[i] << "\t" << lower_y[i] << "\t" << higher_y[i] << endl;
  // }

cout << "end of reading the file " << endl;

    // TGraphAsymmErrors (Int_t n, const Double_t *x, const Double_t *y, const Double_t *exl=0, const Double_t *exh=0, const Double_t *eyl=0, const Double_t *eyh=0)
TGraphAsymmErrors *Verbinski = new TGraphAsymmErrors(nPoints,energy,spectrum,0,0,lower_y,higher_y);

//TGraph *Verbinski = new TGraph(nPoints,energy,spectrum);

//////////////////////////////////////////
///////////////////////////////////////////



    // Correct for neutron contribution in CACTUS
    // First, find scaling factor for neutron correction function so it contributes to 22% of total count
    // double totalNeutronCount = 0;
    // for(ix=0; ix<dx; ix++){ // Loop over all gamma energies
    //     double binCenterEnergy = (x_cal[0] + x_cal[1]*(ix+0.5)); // (lowest bincenter + deltaX * binNumber)
    //     if(binCenterEnergy >= Egammamin && binCenterEnergy < Egammamax) { // we only work with energies within this range
    //         totalNeutronCount += exp(-2.0*binCenterEnergy*0.001); // See where neutronCorrectionFunctionScaling is declared in this file . Also scale to MeV
    //     }
    // }

    // // We want that neutronCorrectionFunctionScaling*totalNeutronCount = oneMinusNeutronCorrectionFactor*number_of_counts (i.e. 22%)
    // neutronCorrectionFunctionScaling = oneMinusNeutronCorrectionFactor*number_of_counts / totalNeutronCount;

    // for(ix=0; ix<dx; ix++){ // Loop over all gamma energies
    //     // Gamma energy of this bin
    //     double binCenterEnergy = (x_cal[0] + x_cal[1]*(ix+0.5)); // (lowest bincenter + deltaX * binNumber) * 0.001 to convert from keV to MeV
    //     double neutronCount = neutronCorrectionFunctionScaling*exp(-2.0*binCenterEnergy*0.001);
    //     if(binCenterEnergy >= Egammamin && binCenterEnergy < Egammamax) {
    //         double totalCount = 0;

    //         for(iy=0;iy<dy;iy++) { // Loop over all excitation energies
    //             double countInBin = matrix->GetBinContent(ix,iy);
    //             totalCount += countInBin;
    //         }

    //         double fraction = neutronCount / totalCount;
    //         // cout << "E_γ = " << binCenterEnergy << ", totalCount = " << totalCount << " and neutronCount = " << neutronCount << ", fraction = " << fraction << endl;

    //         for(iy=0;iy<dy;iy++) { // Loop over all excitation energies
    //             double countInBin = matrix->GetBinContent(ix,iy);
    //             // chose on of the two follwing options for neutron correction ("countInBin"):
    //             countInBin *= (1.0 - fraction); // Per gamma ray neutron scaling
    //             // countInBin *= neutronCorrectionFactor; // Constant neutron scaling

    //             // matrix->SetBinContent(ix,iy,countInBin);
    //         }
    //     }
    // }

    
    cout << "\n F.g. matrix for XXX1 is now filled.\n" << endl;

    //================== Offline Sort and                ==================
    //================== Creating of the first Multiplicites===============



    TFile SortingFile("../offline_Pu239_dp_039.root");
    // TFile SortingFile("../offline_233U_deutrons_060616_averagetimegates.root");
        
 
    TH1F *h1 = (TH1F*)SortingFile.Get("h_ex_fiss");
    // TH1I *h1cp = h1-> Clone(); 
    TH2F *h2 = (TH2F*)SortingFile.Get("m_alfna_fiss_neutron_corrected"); // denne må forandres ti den neutron corrected
    // TH2F *h3a = h2-> Clone();

    // Neutron/gamma discrimination
    // for(int binIndex=1; binIndex<=numberOfBins; binIndex++) { // Loop through all bins
    //     double value = h4[i]->GetBinContent(binIndex); 

    /////////////////
    // Warnings that will be created due to the assumptions made about the
    // spectra size
    if( h2->GetXaxis()->GetXmin() != matrix->GetXaxis()->GetXmin() ) {
        cout<< "warning: 
                Xmin of the spectra is not the same -> revise code"
            <<  endl;
        }
    
     if( h2->GetYaxis()->GetXmin() != matrix->GetYaxis()->GetXmin() ) {
        cout<< "warning: 
               Ymin of the spectra is not the same -> revise code"
            <<  endl; 
        }

    if( h2->GetSize() != matrix->GetSize() ) {
       cout<< "warning: 
               nBin on of the matrixes is not the same-> revise code"
           <<  endl; 
        } 
    /////////////////    

    int binmin, binmax;
    TAxis *yaxis_h2 = h2->GetYaxis();
    TAxis *xaxis_h2 = h2->GetXaxis();
    TAxis *xaxis_h1 = h1->GetXaxis();
    binmin = xaxis_h2->FindBin(Egammamin);
    // binmax = xaxis_h1->FindBin(Egammamax);

    // TH1D *h3 = h2->ProjectionY("Fission coincidences (_py)",binmin,-1);
    // TH2F *h4a = h2-> Clone();

    TH1D *h3_unf = matrix->ProjectionY("Fission coincidences (_py), unf.",binmin,-1);
    // TH2F *h4a = h2-> Clone();

    TH1F *hdiv=new TH1F("hdiv","Multiplicity (raw Spectrum)",h1->GetSize(), 
                        h1->GetXaxis()->GetXmin(), 
                        h1->GetXaxis()->GetXmax() );
    // TH1I *hdiv = h1-> Clone(); 
    // hdiv->Divide(h3,h1,1.,1.,"");

    TH1F *hdiv_unf=new TH1F("hdiv_unf","Multiplicity (unfolded spectrum)",h1->GetSize(), 
                        h1->GetXaxis()->GetXmin(), 
                        h1->GetXaxis()->GetXmax() );

    // TH1*  h_multi_1 = new TH1F("h_multi_1","Multiplicities Histogram (raw)",
    //                             nIntervals,Emin,Emax);

    /////////////////////////////////////////////
    // Finding the Std Err from m_alfna_fiss
    /////////////////////////////////////////////
    TH2F *m_alfna_fiss = (TH2F*)SortingFile.Get("m_alfna_fiss"); // denne må forandres ti den neutron corrected
    int gammaBins = m_alfna_fiss->GetXaxis()->GetLast();
    double std_in_Interval[nIntervals][gammaBins];
    // cout << "last" << xaxis_h2->GetLast() << endl;
     // rest within the loop below

    
    int i;
    float xmin[nIntervals], xmax[nIntervals];
    float xinterval;
    double integral, integralErr, scalingfactor;
    TH1D *h4[nIntervals];
    TH1D *h4_err[nIntervals];
    string histname_base = "Fiss. coinc.";
    string string_result;
    char* histname_result;
    std::stringstream sstm;

    for (i=0; i<nIntervals; i++) {
        // define the necessary variables
        xinterval = Erange / nIntervals;
        xmin[i] = Emin + i * xinterval;
        xmax[i] = xmin[i] + xinterval;
        meanEnergies[i] = 0.5*(xmax[i]+xmin[i]);
        binmin = yaxis_h2->FindBin(xmin[i]);
        binmax = yaxis_h2->FindBin(xmax[i]);


        // create the projections with correct title
        sstm.str(std::string());
        sstm << histname_base << " Emin " << xmin[i] << "keV Emax " 
             << xmax[i] <<"keV";
        string_result = sstm.str();
        histname_result = string_result.c_str();
        h4[i] = matrix->ProjectionX(histname_result,binmin,binmax);
        std::cout << "Working with bin " << i << " of " << nIntervals << " with bin limits " << binmin << " and " << binmax << std::endl;


        ///////////////////////////////////////////////
        // Rebin with bins of 100 keV
        double binSize = 90.; // 8 keV per bin            
        int numberOfBins = 14000;
        double lowestBin = -2008;
        // double lowestBin = 0;
        Double_t binLowerLimits[numberOfBins+1];
        for(int binIndex=0; binIndex<=numberOfBins; binIndex++) {
            binLowerLimits[binIndex] = lowestBin + binSize*binIndex; // First will be 0, then 100, then 200 etc
            if(binLowerLimits[binIndex] > numberOfBins-binSize) {
                numberOfBins = binIndex;
                break;
            }
        }


        h4[i] = dynamic_cast<TH1D*>(h4[i]->Rebin(numberOfBins, "", binLowerLimits));
        for(int binIndex=0; binIndex<h4[i]->GetNbinsX(); binIndex++) {
            // After rebinning, root might make an error that isn't correct due to averaging.
            // h4[i]->SetBinError(binIndex, 0.0);
        }
        ////////////////////////////////////////////////



        ////////
        // The Fission
        //start to normalize the projections
        binmin = xaxis_h1->FindBin(xmin[i]);
        binmax = xaxis_h1->FindBin(xmax[i]);
        integral = 0.;
        integralErr = 0.;
        h1->Sumw2();
        integral = h1->IntegralAndError(binmin,binmax,integralErr);

        // cout << integral << endl;
        nFissions[i] = integral;
        nFissionsErr[i] = integralErr;
        cout << "nFissions[i]: " << nFissions[i] << " +- " << nFissionsErr[i] << endl;
        ////////////

        // Further on
        // creating the m_alfna-fiss projection and then loop fill the errors array
        //TH1D *  ProjectionX (const char *name="_px", Int_t firstybin=0, Int_t lastybin=-1, Option_t *option="") const 
        TH1D *m_alfna_fiss_px = m_alfna_fiss->ProjectionX("m_alfna_fiss_px", binmin, binmax);

        m_alfna_fiss_px = dynamic_cast<TH1D*>(m_alfna_fiss_px->Rebin(numberOfBins, "", binLowerLimits));
        for(int binIndex=0; binIndex<m_alfna_fiss_px->GetNbinsX(); binIndex++) {
            // After rebinning, root might make an error that isn't correct due to averaging.
            // h4[i]->SetBinError(binIndex, 0.0);
        }
        m_alfna_fiss_px->Sumw2();
        double scalingfactor = 1/efficiency;
        m_alfna_fiss_px->Scale(1/efficiency);

        double binValue = 0.;
        double binValue_old = 0;
        double binErr = 0.;

        for(int binIndex=0; binIndex < numberOfBins; binIndex++) {
            binGamma = binIndex;
        // for(int binGamma = 0; binGamma < gammaBins; binGamma++) {
            // cout << "binGamma" << binGamma << endl;
        std_in_Interval[i][binGamma] = (m_alfna_fiss_px->GetBinError(binGamma));
        binValue = h4[i]->GetBinContent(binGamma);
        binValue_old = binValue;
        binValue = binValue /efficiency /nFissions[i];
        h4[i]->SetBinContent(binGamma,binValue);
        // and set the error 
        // h4[i]->SetBinError(binGamma,std_in_Interval[i][binGamma]);
        binErr = std_in_Interval[i][binGamma];
        if(binValue!=0){  // allow no devision  by 0
        binErr=  binValue * sqrt( pow(binErr/binValue_old,2.) +  pow(efficiencyErr/efficiency,2.) +pow(nFissionsErr[i]/nFissions[i],2.)   );
        }
        else { binErr = 0;}
        ///////////
        // Important ASSUMPTION: Set Error=Value if no counts in original spectrum (not unfolded spectrum)
        if(std_in_Interval[i][binGamma]==0)
        {
         binErr=  binValue;   
        }
        //////////
        h4[i]->SetBinError(binGamma,binErr);
        // if (i==3){
        
        // cout << "Bin: " << binGamma << ", " << binLowerLimits[binGamma] << "\t" << "Value: " << m_alfna_fiss_px->GetBinContent(binGamma)  << "\t" << "StdErr: " << std_in_Interval[i][binGamma] << endl;
        // cout << "Bin: " << binGamma << ", " << binLowerLimits[binGamma] << "\t" << "Value: " << binValue  << "\t" << "BinErr: " << binErr << endl;

        // cout << "Bin: " << binGamma << ", " << m_alfna_fiss_px->GetBinCenter(binGamma) << "\t" << "Value_be: " << m_alfna_fiss_px->GetBinContent(binGamma)  << "\t" << "StdErr: " << std_in_Interval[i][binGamma] << endl;
        // cout << "Bin: " << binGamma << ", " << m_alfna_fiss_px->GetBinCenter(binGamma) << "\t" << "Value_un: " << binValue  << "\t" << "BinErr: " << binErr << endl;
        
        // }
        }




        //find number of entries (used for the multiplicity calc.)
        binmin = h4[i]->FindBin(Egammamin);
        binmax = h4[i]->FindBin(Egammamax);
        integral = 0;
        integralErr = 0;
        // IntegralanError computes the integral and puts the error in the varibale "error
        // IntegralAndError(Int_t binx1, Int_t binx2, Double_t & error, Option_t *option)
        integral = h4[i]->IntegralAndError(binmin,binmax,integralErr);
        cout << "Integral (total count) for interval " << i << ": " << integral << "+-" << integralErr << endl;
        // cout << integral << endl;
        // nGammaCoinc[i] = 1/efficiency * integral;
        // nGammaCoincErr[i] = nGammaCoinc[i] * sqrt( pow(efficiencyErr/efficiency,2) +  pow(integralErr/integral,2) );
        // cout << "nGammaCoinc[i]: " << nGammaCoinc[i] << "+-" << nGammaCoincErr[i] << endl;

        // WE WILL DO THIS AFTER REBINNING INSTEAD
        
        // old normalization; Now normalization follows further down
        //                    - it will the Multiplicity which is equivalen to the total area!
        // scalingfactor = 0;
        // scalingfactor = 1/integral * 1/efficiency;
        // // scalingfactor = 1 ;
        // // h4[i]->Sumw2();
        // h4[i]->Scale(scalingfactor);
                

        // AverageEnergyErr[i] = h4[i]->GetMeanError(1);


        // calculate the multiplicity (before the gamma correction/extrapolation below the peak)
        // multiplicity[i] = nGammaCoinc[i] / nFissions[i];
        // multiplicityErr[i] = multiplicity[i] * sqrt( pow(nGammaCoincErr[i]/nGammaCoinc[i],2) +  pow(nFissionsErr[i]/nFissions[i],2) );
        multiplicity[i] = integral;
        multiplicityErr[i] = integralErr;
        cout << "multiplicity[i]: " << multiplicity[i] << "+-" << multiplicityErr[i] << endl;


        // float ta=nGammaCoinc[i];
        // float tb=nFissions[i];
        // cout << "nGammaCoinc: " << ta << endl;
        // cout << "nFissions: " << tb << endl;

        // float multiplicityErr[i] = multiplicity[i]*sqrt((1/ta)+(1/tb));
        // cout << xmin[i]/1000 << "-" << xmax[i]/1000 << " MeV" << "\t" 
        //      << multiplicity[i] << " "<<multiplicityErr[i] << endl;
        
        if(multiplicity[i] < 0){
            multiplicity[i] = 0;
            cout << "warning " << "multiplicity for" << xmin[i] 
                << "is < 0; now set to 0!" << endl;
        }

        // h_multi_1->Fill(xmin[i]+10,            //in order to remain in the bin 
        //                 multiplicity[i]);      //for difficult numbers like 5.33

        // double currentBinSize = 8;
        // double currentEnergy = -2008;
        // double nextEnergyRebinning = 3500;

        // for(int binIndex=0; binIndex<=numberOfBins; binIndex++) {
        //     if(currentEnergy>nextEnergyRebinning) {
        //         currentBinSize = 16; // currentBinSize = currentBinSize*2;
        //         nextEnergyRebinning = 1000000;
        //     }
        //     binLowerLimits[binIndex] = currentEnergy;
        //     // cout << binLowerLimits[binIndex] << ", ";
        //     currentEnergy += currentBinSize;
        //     if(currentEnergy > 14000) {
        //         numberOfBins = binIndex;
        //     }
        // }

        // ///////////////////////////////////////////////
        // // Rebin with bins of 100 keV
        // double binSize = 8; // 8 keV per bin            
        // int numberOfBins = 14000;
        // double lowestBin = -2008;
        // Double_t binLowerLimits[numberOfBins+1];
        // for(int binIndex=0; binIndex<=numberOfBins; binIndex++) {
        //     binLowerLimits[binIndex] = lowestBin + binSize*binIndex; // First will be 0, then 100, then 200 etc
        //     if(binLowerLimits[binIndex] > 13800) {
        //         numberOfBins = binIndex;
        //         break;
        //     }
        // }


        // h4[i] = dynamic_cast<TH1D*>(h4[i]->Rebin(numberOfBins, "", binLowerLimits));
        // for(int binIndex=0; binIndex<h4[i]->GetNbinsX(); binIndex++) {
        //     // After rebinning, root might make an error that isn't correct due to averaging.
        //     // h4[i]->SetBinError(binIndex, 0.0);
        // }
        // ///////////////////////

        binmin = h4[i]->FindBin(0);
        binmax = numberOfBins;
        // Actually, this is the Multiplicity!
        double integralBeforeExtrapolationFixErr     = multiplicityErr[i];
        double integralBeforeExtrapolationFix        = multiplicity[i];
        // double integralBeforeExtrapolationFix = h4[i]->Integral(binmin,binmax); 









        /***********************************
        * BEGIN Extrapolation from max value BEGIN
        * (Fjern hvis faktisk spekter skal vises)
        ***********************************/

        // Find maximum value and which bin that is 
        // double maxValue = h4[i]->GetBinContent(0); // Assume first value is largest and compare to all other values
        // int maxValueBin = 0;

        // for(int binIndex=1; binIndex<=numberOfBins; binIndex++) { // Loop through all bins
        //     double value = h4[i]->GetBinContent(binIndex); 
        //     // cout << binIndex << " value: " << value << ", error: " << error << endl;
        //     if(maxValue < value) {
        //         // New max value, update both maxValue and maxValueBin
        //         maxValue = value;
        //         maxValueBin = binIndex;
        //     }
        // }


        // double integralAfterExtrapolationFixErr =0.;
        // double integralAfterExtrapolationFix = h4[i]->IntegralAndError(binmin,binmax,integralAfterExtrapolationFixErr);
        // // integralAfterExtrapolationFixErr = integralAfterExtrapolationFixErr;
        // // integralAfterExtrapolationFixErr = integralAfterExtrapolationFix 
        // //                                    * sqrt( pow(integralAfterExtrapolationFixErr/integralAfterExtrapolationFix,2)
        // //                                            +  pow(efficiencyErr/efficiency,2) 
        // //                                            +  pow(nFissionsErr[i]/nFissions[i],2));
        // cout << "integralAfterExtrapolationFix: " << integralAfterExtrapolationFix << "+-" << integralAfterExtrapolationFixErr << endl;
      
        TH1D *h4_err[i] = dynamic_cast<TH1D*>( h4[i]->Clone() ); 
        double relativeErrorExtrapolatedBins = 0.3; // estimated error on the extrapolated Bins
        // double binValue = 0;
        // double binErr = 0;
        int EgammaminBin = h4[i]->FindBin(Egammamin);
        int n_binChanged = EgammaminBin - binmin;
        // cout << "binmin" << binmin << endl;
        cout << "EgammaminBin - binmin: " << EgammaminBin-binmin << endl;
        // Set all bins before EGammaMin (threshold) to that value
        for(int binIndex=0; binIndex < EgammaminBin; binIndex++) {
            if(h4[i]->GetBinCenter(binIndex) >= 0) {
               // n_binChanged = EgammaminBin - binIndex;
               // cout << n_binChanged << n_binChanged << endl;
               binValue = h4[i]->GetBinContent(EgammaminBin); // minimum bin == value that should be copies to all of the others
               h4_err[i]->SetBinContent(binIndex, binValue);
               h4[i]->SetBinContent(binIndex, binValue);
               binErr = 0;
               h4[i]->SetBinError(binIndex, binErr);
               binErr = sqrt(n_binChanged) * binValue *relativeErrorExtrapolatedBins;   // n_binChanged is
               h4_err[i]->SetBinError(binIndex, binErr);                          // necessary to get the error cal. right for summations/integrals
                                                                                  // on the error band
               // h4[i]->SetBinContent(binIndex, maxValue);
            }
        }


        // TODO : ERROR FOR TOTAL ENERGY

        double totalEnergy = 0;
        double totalEnergyErr = 0;
        double totalEnergyErrUp = 0;
        for(int binIndex=0; binIndex < numberOfBins; binIndex++) {
            if(h4[i]->GetBinCenter(binIndex) >= 0) {                            // TODO: This or rebin such that 0 is smallest bin?!
                double value = h4[i]->GetBinContent(binIndex);
                double binErr = h4[i]->GetBinError(binIndex);
                double binEnergy = h4[i]->GetBinCenter(binIndex); 
                // double binSize = binLowerLimits[binIndex+1]-binLowerLimits[binIndex];
                totalEnergy += binEnergy * value;
                totalEnergyErr += pow( binEnergy * binErr, 2.);

                binErr = h4_err[i]->GetBinError(binIndex);
                totalEnergyErrUp += pow( binEnergy * binErr, 2.);
                // cout << "binEnergy: " << binEnergy << "\t" << "Value: " << value << "\t binErr: " << binErr << endl;
            }
        }
        TotalEnergy[i] = totalEnergy;
        TotalEnergyErr[i] = sqrt(totalEnergyErr);
        totalEnergyErrUp = sqrt(totalEnergyErrUp);
        TotalEnergyUp[i] = totalEnergy + totalEnergyErrUp;
        TotalEnergyLow[i] = totalEnergy - totalEnergyErrUp;
        cout << "TotalEnergyUp/Low:" << TotalEnergyUp[i] << "\t" << TotalEnergyLow[i] << endl;
        cout << "Total energy NEW: " << TotalEnergy[i] << "+-" << TotalEnergyErr[i] << endl;

        /***********************************
        * END Extrapolation from max value END
        ***********************************/
        double integralAfterExtrapolationFixErr =0.;
        double integralAfterExtrapolationFix = h4[i]->IntegralAndError(binmin,binmax,integralAfterExtrapolationFixErr);
        // integralAfterExtrapolationFixErr = integralAfterExtrapolationFixErr;
        // integralAfterExtrapolationFixErr = integralAfterExtrapolationFix 
        //                                    * sqrt( pow(integralAfterExtrapolationFixErr/integralAfterExtrapolationFix,2)
        //                                            +  pow(efficiencyErr/efficiency,2) 
        //                                            +  pow(nFissionsErr[i]/nFissions[i],2));

        double integralErrUp =0.;
        h4_err[i]->IntegralAndError(binmin,binmax,integralErrUp);

        cout << "integralAfterExtrapolationFix: " << integralAfterExtrapolationFix << "+-" << integralAfterExtrapolationFixErr << endl;

        double multiplicityCorrectionFactor = integralAfterExtrapolationFix / integralBeforeExtrapolationFix;
        // ASSUMPTION: integralBeforeErr << integralAfterErr
        // however, probably this correction factor itself is not used elsewhere
        double multiplicityCorrectionFactorErr = multiplicityCorrectionFactor * sqrt( pow(integralAfterExtrapolationFixErr/integralAfterExtrapolationFix,2) ) ;
        
        cout << "Multiplicity correction factor: " << multiplicityCorrectionFactor << "+-" << multiplicityCorrectionFactorErr << endl;
        multiplicity[i] = integralAfterExtrapolationFix;
        multiplicityErr[i] = integralAfterExtrapolationFixErr;

        multiplicityUp[i]  = multiplicity[i] + integralErrUp;
        multiplicityLow[i] = multiplicity[i] - integralErrUp;

        cout << "Multiplicity after: " << multiplicity[i] << "+-" << multiplicityErr[i] << endl;


        //normalize everything
        binmin = 0;
        binmax = numberOfBins;
        // integral = h4[i]->Integral(binmin,binmax);
        // nFissions[i] = integral;
        // scalingfactor = 1/integral; //old scaling 

        // need this scaling here (!) as we want to transfer the plot from
        // cnt/(binsize*fiss) --*S--> cnt/(MeV*fiss)
        // so S = MeV/binsize, or 1000 * keV/binsize, with binsize in keV
        // + not needed before - because there we sum over the bins, so it's good to give values in cnt/bin
        scalingfactor = (1000/binSize);
        h4[i]->Scale(scalingfactor);
        h4_err[i]->Scale(scalingfactor);
        // cout << "Sunniva AVERAGE = " << h4[i]->GetMean(1) << endl;
        
        // double meanValue = 0;
        // for(int binIndex=0; binIndex<=numberOfBins; binIndex++) {
        //     cout << "Working on bin index = " << binIndex << " with center = " << h4[i]->GetBinCenter(binIndex) << " and value " << h4[i]->GetBinContent(binIndex) << endl;
        //     meanValue += h4[i]->GetBinContent(binIndex) * h4[i]->GetBinCenter(binIndex);
        // }
        // cout << "SunnivaHAX AVERAGE = " << meanValue / (numberOfBins+1) << endl;
        AverageEnergy[i]    = h4[i]->GetMean(1);
        AverageEnergyErr[i] = h4[i]->GetMeanError(1);
        double avEnErrUp    = h4_err[i]->GetMeanError(1);
        AverageEnergyUp[i]  = AverageEnergy[i] + avEnErrUp;
        AverageEnergyLow[i] = AverageEnergy[i] - avEnErrUp;
        cout << "AverageEnergy: " << AverageEnergy[i] << "+-" << AverageEnergyErr[i] << endl;

        // AverageEnergy[i]    = h4_err[i]->GetMean(1);
        // AverageEnergyErr[i] = h4_err[i]->GetMeanError(1);
        // cout << "AverageEnergy: " << AverageEnergy[i] << "+-" << AverageEnergyErr[i] << endl;

        //calculate the total energy
        // TotalEnergy[i] = AverageEnergy[i]*multiplicity[i];    (This is calculated before the scaling)
        // cout << "Total energy OLD = " << AverageEnergy[i]*multiplicity[i] << endl << endl;

        // // At the end of the loop: All calculations should be done;
        // // then we copy the h4_err on h4 thus that we can use the same drawing options as before
        // h4[i] = h4_err[i];
        cout << "\n" << endl;
    }


    // // Multiplicities using the raw matrixes
    // float n_hex, n_alfna_py, n_multi, nsize, nlength, nvalue;
    // for (i=0; i<h1->GetSize(); i++) {
    //     n_hex = h1->GetBinContent(i);
    //     n_alfna_py = h3->GetBinContent(i);
    //     if (n_hex==0){n_multi=0;}
    //     else {n_multi = n_alfna_py / n_hex;}

    //     nsize = h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin(); 
    //     nlength = nsize / h1->GetSize();
    //     nvalue = i*nlength - abs(h1->GetXaxis()->GetXmin()) +2; 
    //     // +2 to get into the bin for nubers
    //     //  like 5.33 (that might be difficult) 

    //     hdiv->Fill(nvalue,n_multi);
    // }

    // Multiplicities using the unfolded matrixes
    float n_hex, n_alfna_py, n_multi, nsize, nlength, nvalue;
    float varscale;
    for (i=0; i<h1->GetSize(); i++) {
        n_hex = h1->GetBinContent(i);
        n_alfna_py = h3_unf->GetBinContent(i);
        if (n_hex==0){n_multi=0;}
        else {n_multi = n_alfna_py / n_hex;}

        nsize = h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin(); 
        nlength = nsize / h1->GetSize();
        nvalue = i*nlength - abs(h1->GetXaxis()->GetXmin()) +2; 
        // +2 to get into the bin for nubers
        //  like 5.33 (that might be difficult) 

        hdiv_unf->Fill(nvalue,n_multi);
    }
    varscale = 1/efficiency;
    hdiv_unf->Scale(varscale);


    //  ///// Histogram over the muliplicites
        



    // ================== Start the Game ==================

    c1_1.cd();
    c1_1.SetLogy();
    c1_1.SetLeftMargin(0.14);
    c1_1.SetRightMargin(0.01);
    c1_1.SetBottomMargin(0.15);

    h1.GetXaxis().CenterTitle();
    h1.GetXaxis().SetTitle("Excitation energy E_{x} [keV]");
    h1.GetYaxis().CenterTitle();
    h1.GetYaxis().SetTitleOffset(1.1);
    h1.GetXaxis().SetTitleOffset(1.4);
    h1.GetYaxis().SetTitle("Counts/keV");
    h1.GetYaxis().SetTitleSize(0.05);
    h1.GetXaxis().SetTitleSize(0.05);
    h1.GetYaxis().SetLabelSize(0.06);
    h1.GetXaxis().SetLabelSize(0.06);
    h1.GetXaxis().SetRangeUser(2000,13e3);
    h1.GetYaxis().SetRangeUser(10,5000);
    h1.SetLineColor(kBlack);
    h1.SetLineWidth(2);
    h1.Scale(1/x_cal[1]);
    // cout << "calibration" << x_cal[1] << endl;
    gStyle->SetOptStat(0);
    h1.Draw();
    
    h3_unf->SetLineColor(kRed);
    h3_unf.Scale(1/x_cal[1]);
    h3_unf->Draw("same");

    TLegend *leg = new TLegend(0.67,.733051, 0.997126, 0.868644);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg->AddEntry(h1," ^{239}Pu(d,pf) ","L");
    leg->AddEntry(h3_unf," ^{239}Pu(d,pf#gamma) ","L");
    leg->SetTextSize(0.06);
    leg->SetTextFont(42);
    leg->Draw();


    // TLatex t;
    // t.SetTextSize(0.06);
    // t.SetTextFont(42);
    // t.DrawLatex(    1000,8000,"(a)");
    // t.SetTextSize(0.05);
    // t.DrawLatex(    6800,3700,"S_{n}");
    
    TLine *line1 = new TLine(Emin,0,Emin,100);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->Draw();

    TLine *line4 = new TLine(Emax,0,Emax,100);
    line4->SetLineStyle(2);
    line4->SetLineWidth(2);
    line4->Draw();

    //SN and Fiss barrier
    double Sn      = 6534;
    double FissBar = 605*10;
    TArrow *line2 = new TArrow(FissBar,100,FissBar,230,0.02,">");
    line2->SetLineStyle(1);
    line2->SetLineWidth(2);
    line2->Draw();


    TArrow *line3 = new TArrow(Sn,100,Sn,230,0.02,">");
    line3->SetLineStyle(1);
    line3->SetLineWidth(2);
    line3->Draw();

    TLatex t;
    t.SetTextSize(0.03);
    t.DrawLatex(Sn-200,90,"S_{n}");
    t.DrawLatex(FissBar-300,90,"B_{F}");


    c1_1->Update();
    c1_1->Print("multiplicitiesA1.pdf");


    TCanvas *c1_2 = new TCanvas("c1_2","Fission Gamma Spectra_All");
    // c1->Divide(1,1,0,0);
    
    c1_2.cd();
    c1_2.SetLogy();
    c1_2.SetLeftMargin(0.13);
    c1_2.SetRightMargin(0.01);
    c1_2.SetBottomMargin(0.15);
    h4[0].GetXaxis().CenterTitle();
    h4[0].GetXaxis().SetTitle("E_{#gamma} [keV]");
    h4[0].GetXaxis().SetTitleSize(0.06);
    h4[0].GetXaxis().SetTitleOffset(1.1);
    h4[0].GetYaxis().CenterTitle();
    
//    h4[0].GetXaxis().SetLabelOffset(1.4);

    h4[0].GetYaxis().SetTitle("Photons/Fission*MeV");
    h4[0].GetYaxis().SetTitleSize(0.07);
    h4[0].GetYaxis().SetTitleOffset(1.2);
    h4[0].GetYaxis().SetTitleSize(0.053);
    h4[0].GetXaxis().SetLabelSize(0.053);
    h4[0].GetYaxis().SetLabelSize(0.053);
    h4[0].GetXaxis().SetRangeUser(0,15e3);
    h4[0].GetXaxis().SetRangeUser(0,13800);
    h4[0].GetYaxis().SetRangeUser(0.00017,10);
    h4[0].SetLineColor(2);
    // h4[0].SetLineWidth(2);
     gStyle->SetOptStat(0);
    h4[0].Draw();


    
    for (i=2; i<nIntervals; i++) {
    h4[i]->SetLineColor(i+1);
    h4[i]->Draw("same");
    }


    Verbinski->SetMarkerStyle(33);
    // Verbinski->SetLineStyle(2);
    Verbinski->SetMarkerSize(0.8); 
    // Verbinski->SetLineWidth(2.8); 
    Verbinski->SetMarkerColor(kBlue); 
    Verbinski->SetLineColor(kBlue+2);

    // Verbinski->SetMarkerColor(1);
    // Verbinski->SetMarkerSize(20); 
    // Verbinski->Draw("P");

    TLine *line1 = new TLine(Egammamin,0,Egammamin,1e4);
    line1->SetLineStyle(2);
    line1->SetLineWidth(3);
    // line1->Draw();

    TLatex t;
    t.SetTextSize(0.05);
    // t.DrawLatex(-130,-0.0014,"Threshold");

    TArrow *arrow1 = new TArrow(110,-0.00089,Egammamin,-0.0002,0.02,">");
    // arrow1->Draw();
    
//     h6->SetLineColor(kRed);
// //    h6->Draw("same");
    
//     recalhist->SetLineColor(kRed);
//     recalhist->Draw("same");
    
    TLegend *leg = new TLegend(0.567529,0.563559  ,0.968391,0.891949);
    leg.SetBorderSize(0);
    leg.SetFillStyle(1);
    

    for (i=0; i<nIntervals; i++) {

        // xinterval = Erange / nIntervals;
        // xmin[i] = Emin + i * xinterval;
        // xmax[i] = xmin[i] + xinterval;

        sstm.str(std::string());
        sstm << " E_{x} = " << std::setprecision(3) << xmin[i]/1000 << " - " << xmax[i]/1000 <<" MeV";
        string_result = sstm.str();
        histname_result = string_result.c_str();
        leg->AddEntry(h4[i],histname_result,"L");
    }
    cout << "\n" << endl;

    // leg->AddEntry(Verbinski,"Verbinski","ep");

    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->Draw();

    
//     t.SetTextSize(0.053);
//     t.SetTextFont(42);
//     t.DrawLatex(    1000,2000000,"(b)");
//     t.SetTextSize(0.046);
//     t.DrawLatex(    9200,14000,"S_{n}");
    
//     t.SetTextAngle(90);
// //    t.DrawLatex(    -200,2000000,"counts");
    
    
//     TLine *line1 = new TLine(9398.1,400,9398.1,10000);
//     line1->SetLineStyle(2);
//     line1->Draw();

    c1_2->Update();
    
    c1_2->Print("multiplicitiesA2.pdf");


/////// Histogram over the muliplicites //////////////////

    //THEORY definitions
    ////////////////////////////////////////////////////////////////////////
    ifstream theory_ex_energy("Theory_Christelle/ex_energy_THEORY.txt");
    ifstream theory_multiplicity("Theory_Christelle/multiplicity_THEORY.txt");
    ifstream theory_energy_av("Theory_Christelle/energy_av_THEORY.txt");
    ifstream theory_energy_tot("Theory_Christelle/energy_tot_THEORY.txt");

    int nLs = 3;
    int nArrays = nLs +1; // +1 for energy

    float THEORY_ex_energy[100];
    float THEORY_multiplicity[nArrays][100], THEORY_energy_av[nArrays][100], THEORY_energy_tot[nArrays][100];

    //reading files
    // int numberOfValuesInFile = 8;
    int i=0;
    int n_headerlines = 1;

    string line;

    // while(i < numberOfValuesInFile) {
    // 	if (i>n_headerlines){j=i-n_headerlines;}
    // 	else j=0;
    //     theory_ex_energy >> THEORY_ex_energy[j];
    //   i++;
    // }

    i=0;

    while(getline(theory_multiplicity, line)) 
    {
        //the following line trims white space from the beginning of the string
        //line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace)))); 

        if(line[0] == '#') continue;

        stringstream(line) //>> dummy
        					>> THEORY_multiplicity[0][i] 
                         	>> THEORY_multiplicity[1][i] 
                         	>> THEORY_multiplicity[2][i]
                         	>> THEORY_multiplicity[3][i];
        // cout << THEORY_multiplicity[0][i] << " " << THEORY_multiplicity[1][i]  << endl;
        THEORY_multiplicity[0][i] *= 1000; // Multiply by 1000 to get keV from MeV
      i++;
    }
    int numberOfValuesMultiplicity = i;

    i=0;
    while(getline(theory_energy_av,line)){

    	if(line[0] == '#') continue;

        stringstream(line) //>> dummy
        				 >> THEORY_energy_av[0][i] 
                         >> THEORY_energy_av[1][i] 
                         >> THEORY_energy_av[2][i]
                         >> THEORY_energy_av[3][i];
        THEORY_energy_av[0][i] *= 1000; // Multiply by 1000 to get keV from MeV
        THEORY_energy_av[1][i] *= 1000; // Multiply by 1000 to get keV from MeV
        THEORY_energy_av[2][i] *= 1000; // Multiply by 1000 to get keV from MeV
        THEORY_energy_av[3][i] *= 1000; // Multiply by 1000 to get keV from MeV
      i++;
    }
    int numberOfValuesEnergyAv = i;
    // for(i=0;i<numberOfValuesInFile;i++){
    // 	cout<< THEORY_multiplicity[0][i] << " " << THEORY_multiplicity[1][i] << endl;
    // }
    
    i=0;
    while(getline(theory_energy_tot,line)){

    	if(line[0] == '#') continue;
        stringstream(line) //>> dummy
        				 >> THEORY_energy_tot[0][i] 
                         >> THEORY_energy_tot[1][i] 
                         >> THEORY_energy_tot[2][i]
                         >> THEORY_energy_tot[3][i];
       THEORY_energy_tot[0][i] *= 1000; // Multiply by 1000 to get keV from MeV
       THEORY_energy_tot[1][i] *= 1000; // Multiply by 1000 to get keV from MeV
       THEORY_energy_tot[2][i] *= 1000; // Multiply by 1000 to get keV from MeV
       THEORY_energy_tot[3][i] *= 1000; // Multiply by 1000 to get keV from MeV
      i++;
    }
    int numberOfValuesEnergyTot = i;

TGraph *theory_multiplicity_graph[nLs];
TGraph *theory_energy_av_graph[nLs];
TGraph *theory_energy_tot_graph[nLs];

    for(i=0;i<nLs;i++){
    theory_multiplicity_graph[i] = new TGraph(numberOfValuesMultiplicity,THEORY_multiplicity[0],THEORY_multiplicity[i+1]);
    theory_energy_av_graph[i] = new TGraph(numberOfValuesEnergyAv,       THEORY_energy_av[0],   THEORY_energy_av[i+1]);
    theory_energy_tot_graph[i] = new TGraph(numberOfValuesEnergyTot,     THEORY_energy_tot[0],  THEORY_energy_tot[i+1]);
    }
    
    

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

    //EXPERIMENTAL Multiplicities
   Double_t w = 600;
   Double_t h = 800;
    TCanvas *c2 = new TCanvas("c2","Fission Gamma Multiplicities",w,h);
    c2.SetRightMargin(0.01);
    c2.SetBottomMargin(0.15);
    c2.SetLeftMargin(0.12);
    c2->Divide(1,3,0,0);

    

    c2_1.SetRightMargin(0.008);
	c2_2.SetRightMargin(0.008);
    c2_3.SetRightMargin(0.008);
    
    c2.cd(3);
    // c2.SetLogy();
    // c2.SetLeftMargin(0.14);
    // c2.SetRightMargin(0.5);
    // c2.SetBottomMargin(0.2);
    gStyle->SetOptStat("uo");

    


    // TGraphErrors*  h_multi_1 = new TGraphErrors(nIntervals, meanEnergies, multiplicity, 0, multiplicityErr);   //errors tar ikke hensyn til usikkerhet i relative netron contribution
    TGraphErrors*  h_multi_1 = new TGraphErrors(nIntervals, meanEnergies, multiplicity, 0, multiplicityErr);  //errors med hensyn på relative neutron contribution
    TGraph *grlow = new TGraph(nIntervals, meanEnergies, multiplicityLow);
    TGraph *grup = new TGraph(nIntervals,  meanEnergies, multiplicityUp);
    // TGraph*  h_multi_1 = new TGraph(nIntervals, meanEnergies, multiplicity);  //errors med hensyn på relative neutron contribution
    h_multi_1->SetMarkerColor(kRed);
    h_multi_1->SetMarkerStyle(21);
    
    h_multi_1.GetXaxis().CenterTitle();
    h_multi_1.GetXaxis().SetTitle("Excitation energy E_{x} [keV]");
    h_multi_1.GetXaxis().SetTitleSize(0.07);
    h_multi_1.GetYaxis().CenterTitle();
    h_multi_1.GetXaxis().SetTitleOffset(1);
    h_multi_1.GetYaxis().SetTitle("Multiplicity");
    h_multi_1.GetYaxis().SetTitleSize(0.069);
    h_multi_1.GetYaxis().SetTitleOffset(0.85);
    h_multi_1.GetYaxis().SetLabelSize(0.06);
    h_multi_1.GetXaxis().SetLabelSize(0.06);
    double ymin = 5.;
    h_multi_1.GetYaxis().SetRangeUser(5.1,10.7);
    // h_multi_1.GetXaxis().SetRangeUser(4.6e3,9.9e3);

    // h1cp->Scale(2e-4);    // the scaling here is arbitrary!
    // h1cp->Draw("same"); 
    // h3->Scale(1e-4);    // the scaling here is arbitrary!
    // h3->Draw("same");
    // hdiv->SetLineColor(kBlue);
    // hdiv->Draw("same");
    h_multi_1->Draw("AP");
    
    int colour=0;
    int colour0 = 4;
    int colour1 = 8;
    int colour2 = 28;

    for(i=0;i<nLs;i++){
    	     if(i==0) {colour = colour0;}
    	else if(i==1) {colour = colour1;}
    	else if(i==2) {colour = colour2;}
    theory_multiplicity_graph[i]->SetMarkerStyle(21);
    theory_multiplicity_graph[i]->SetLineStyle(2);
    theory_multiplicity_graph[i]->SetMarkerSize(0.8); 
    theory_multiplicity_graph[i]->SetLineWidth(2.8); 
    theory_multiplicity_graph[i]->SetMarkerColor(colour); 
    theory_multiplicity_graph[i]->SetLineColor(colour);
    theory_multiplicity_graph[i]->Draw("P");
    theory_multiplicity_graph[i]->Draw("L");
    }

   grlow->SetLineStyle(4);
   grlow->Draw("same");

   grup->SetLineStyle(4);
   grup->Draw("same");


   TGraphErrors *OtherExp = new TGraphErrors(1);
   // OtherExp->SetFillColor(23);
   OtherExp->SetMarkerColor(1);
   OtherExp->SetMarkerStyle(23);
   OtherExp->SetMarkerSize(1);
   double y_value = 7.24;
   double y_err   = 0.7;
   //   (PointNumber , X , Y)
   OtherExp->SetPoint     (0,6540,y_value);
   OtherExp->SetPointError(0, 0, y_err);
   OtherExp->Draw("P");
    // theory_multiplicity_graph->SetMarkerColor(kBlue);
    

    hdiv_unf->SetLineColor(kRed);
    // hdiv_unf->Draw("same p");

    TLegend *leg = new TLegend(0.75,0.767,0.80,0.90);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    // leg->AddEntry(h_multi_1,"Hist. of Multipl.","L");
    // // leg->AddEntry(hdiv," (d,pf#gamma)/(d,pf) = Mult.,raw","L");
    // leg->AddEntry(hdiv_unf," (d,pf#gamma)/(d,pf) = Mult.","L");
    // // leg->AddEntry(h3," ^{239}Pu(d,pf#gamma), arb. scaled ","L");
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->Draw();




////////////
    TLine *line1 = new TLine(FissBar,ymin,FissBar,1e5);
    // line1->SetLineStyle(2);
    line1->SetLineWidth(1);
    line1->Draw();

    TLine *line4 = new TLine(Sn,ymin,Sn,1e5);
    // line4->SetLineStyle(2);
    line4->SetLineWidth(1);
    line4->Draw();
////////////
    // TLatex t;
    // t.SetTextSize(0.03);
    // t.DrawLatex(Sn-200,90,"S_{n}");
    // t.DrawLatex(FissBar-300,90,"B_{F}");
///////////


    // c2->Update();
    
    // c2->Print("multiplicitiesB.pdf");

    ///////////////////////////////////////////////////////
    // AVERAGE gamma energy as function of excitation energy
//     ///////////////////////////////////////////////////////

   // TCanvas *c4 = new TCanvas("c4","Average Gamma Energies");
    // c4.SetLeftMargin(0.14);
    // c1->Divide(1,1,0,0);
c2.cd(2);

// TEST ONLY
    TGraphErrors *gr = new TGraphErrors(nIntervals, meanEnergies, AverageEnergy, 0, AverageEnergyErr);
    TGraph *grlow = new TGraph(nIntervals, meanEnergies, AverageEnergyLow);
    TGraph *grup = new TGraph(nIntervals,  meanEnergies, AverageEnergyUp);

   // TGraphErrors *gr = new TGraphErrors(8, meanEnergies, AverageEnergy, 0, AverageEnergyErr);
    // TGraph *gr = new TGraph(nIntervals, meanEnergies, AverageEnergy);    
   // TGraphErrors *gr = new TGraphErrors(8, meanEnergies, AverageEnergy, 0, 0);    
	gr.GetYaxis().SetRangeUser(705,1.47e3);
	gr->SetTitle("Sunniva er hot");
	gr->GetYaxis().CenterTitle();

	gr->GetYaxis().SetTitle("Average #gamma-energy [keV]");
	gr->GetYaxis().SetTitleSize(0.078);
	gr->GetYaxis().SetTitleOffset(0.73);
	gr.GetYaxis().SetLabelSize(0.06);
	gr.GetXaxis().SetLabelSize(0.06);
	gr->SetMarkerColor(kRed);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");

   grlow->SetLineStyle(4);
   grlow->Draw("same");

   grup->SetLineStyle(4);
   grup->Draw("same");

    for(i=0;i<nLs;i++){
    	     if(i==0) {colour = colour0;}
    	else if(i==1) {colour = colour1;}
    	else if(i==2) {colour = colour2;}
    theory_energy_av_graph[i]->SetMarkerStyle(21);
    theory_energy_av_graph[i]->SetLineStyle(2);
    theory_energy_av_graph[i]->SetMarkerSize(0.8); 
    theory_energy_av_graph[i]->SetLineWidth(2.8); 
    theory_energy_av_graph[i]->SetMarkerColor(colour); 
    theory_energy_av_graph[i]->SetLineColor(colour);
    theory_energy_av_graph[i]->Draw("P");
    theory_energy_av_graph[i]->Draw("L");
    }

   TGraphErrors *OtherExp = new TGraphErrors(1);
   // OtherExp->SetFillColor(23);
   OtherExp->SetMarkerColor(1);
   OtherExp->SetMarkerStyle(23);
   OtherExp->SetMarkerSize(1);
   double y_value = 970;
   double y_err   = 50;
   //   (PointNumber , X , Y)
   OtherExp->SetPoint     (0,6540,y_value);
   OtherExp->SetPointError(0, 0, y_err);
   OtherExp->Draw("P");

    ////////////
    TLine *line1 = new TLine(FissBar,0,FissBar,1e5);
    // line1->SetLineStyle(2);
    line1->SetLineWidth(1);
    line1->Draw();

    TLine *line4 = new TLine(Sn,0,Sn,1e5);
    // line4->SetLineStyle(2);
    line4->SetLineWidth(1);
    line4->Draw();
////////////

   // c4->Print("AvgE_gamma.pdf");
//     TCanvas *c4 = new TCanvas("c4","Average Gamma Energies",900,600);



 ///////////////////////////////////////////////////////
    // TOTAL gamma energy as function of excitation energy
//     ///////////////////////////////////////////////////////

   // TCanvas *c5_1 = new TCanvas("c5_1","Total Gamma Energies");
   //  c5_1.SetLeftMargin(0.14);
c2.cd(1);
   // TGraphErrors *gr = new TGraphErrors(8, meanEnergies, TotalEnergy, 0, AverageEnergyErr);   //errors tar ikke hensyn til usikkerhet i relative netron contribution
   TGraphErrors *gr = new TGraphErrors(nIntervals, meanEnergies, TotalEnergy, 0, TotalEnergyErr); //errors med hensyn på relative neutron contribution
   TGraph *grlow = new TGraph(nIntervals, meanEnergies, TotalEnergyLow);
   TGraph *grup = new TGraph(nIntervals,  meanEnergies, TotalEnergyUp);

   // TGraph *gr = new TGraph(nIntervals, meanEnergies, TotalEnergy); //errors med hensyn på relative neutron contribution
	gr.GetYaxis().SetRangeUser(5.01e3,10.5e3);
	gr->SetTitle("Sunniva er hot");
	gr->GetYaxis().CenterTitle();
	gr->GetYaxis().SetTitle("Total #gamma-energy [keV]");
	gr->GetYaxis().SetTitleSize(0.078);
	gr->GetYaxis().SetTitleOffset(0.73);
	gr.GetYaxis().SetLabelSize(0.06);
	gr.GetXaxis().SetLabelSize(0.06);
	gr->SetMarkerColor(kRed);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");

   grlow->SetLineStyle(4);
   grlow->Draw("same");

   grup->SetLineStyle(4);
   grup->Draw("same");

    for(i=0;i<nLs;i++){
    	     if(i==0) {colour = colour0;}
    	else if(i==1) {colour = colour1;}
    	else if(i==2) {colour = colour2;}
   theory_energy_tot_graph[i]->SetMarkerStyle(21);
   theory_energy_tot_graph[i]->SetLineStyle(2);
   theory_energy_tot_graph[i]->SetMarkerSize(0.8); 
   theory_energy_tot_graph[i]->SetLineWidth(2.8); 
   theory_energy_tot_graph[i]->SetMarkerColor(colour); 
   theory_energy_tot_graph[i]->SetLineColor(colour);
   theory_energy_tot_graph[i]->Draw("P");
   theory_energy_tot_graph[i]->Draw("L");
   }


   TGraphErrors *OtherExp = new TGraphErrors(1);
   // OtherExp->SetFillColor(23);
   OtherExp->SetMarkerColor(1);
   OtherExp->SetMarkerStyle(23);
   OtherExp->SetMarkerSize(1);
   double y_value = 6510;
   double y_err   = 300;
   //   (PointNumber , X , Y)
   OtherExp->SetPoint     (0,6540,y_value);
   OtherExp->SetPointError(0, 0, y_err);
   OtherExp->Draw("P");

    TLegend *leg = new TLegend(0.152,0.582,0.463,0.951);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    // leg.SetBorderSize(2);
    leg->SetTextColor(kBlack);
    leg->SetTextSize(0.07);
    leg->AddEntry(gr,"Present exp.","ep");
    leg->AddEntry(theory_energy_tot_graph[0],"GEF L=0","lp");
    leg->AddEntry(theory_energy_tot_graph[1],"GEF L=4","lp");
    leg->AddEntry(theory_energy_tot_graph[2],"GEF L=10","lp");
    leg->AddEntry(OtherExp,"Verbinski (1973)","ep");
    // leg->AddEntry(grlow,"Error band (incl. systematic)","l");
    leg->Draw();

    ////////////
    TLine *line1 = new TLine(FissBar,0,FissBar,1e5);
    // line1->SetLineStyle(2);
    line1->SetLineWidth(1);
    line1->Draw();

    TLine *line4 = new TLine(Sn,0,Sn,1e5);
    // line4->SetLineStyle(2);
    line4->SetLineWidth(1);
    line4->Draw();
////////////

    // c5_1->Print("TotalE_Gamma.pdf");
    c2->Print("Together.pdf");
// //   c1_2.cd();
//     c4.SetLeftMargin(0.14);
// //     c1_2.SetRightMargin(0.01);
// //     c1_2.SetBottomMargin(0.14);
//     // h4[0].GetXaxis().CenterTitle();
// //     h4[0].GetXaxis().SetTitle("E(NaI) [keV]");
// //     h4[0].GetYaxis().CenterTitle();
// //     h4[0].GetYaxis().SetTitleOffset(1.2);
// // //    h4[0].GetXaxis().SetLabelOffset(1.4);
// //     h4[0].GetYaxis().SetTitle("Counts/fission*MeV");
// //     h4[0].GetXaxis().SetTitleSize(0.053);
// //     h4[0].GetYaxis().SetTitleSize(0.053);
// //     h4[0].GetXaxis().SetLabelSize(0.053);
// //     h4[0].GetYaxis().SetLabelSize(0.053);

// //     h4[0].GetXaxis().SetRangeUser(0,15e3);
//     // h4[0].GetYaxis().SetRangeUser(1000,1700);
// //     h4[0].SetLineColor(kBlack);
// //     h4[0].SetLineWidth(2);
// //     h4[0].Draw();
//     // for (i=1; i<nIntervals; i++) {
//     // h4[i]->GetMean(1);
//     // h4[i]->Draw("same");
//     // }

//     // TGraphErrors *AverageEnergyPlot = new TGraphErrors(8, xmin, AverageEnergy, AverageEnergyErr);
//     TGraph *AverageEnergyPlot = new TGraph(8, xmin, AverageEnergy);
//     AverageEnergyPlot->SetMarkerStyle(2);
//     AverageEnergyPlot->SetMarkerSize(2);
//     AverageEnergyPlot->Draw("P");

// //     TLine *line1 = new TLine(Egammamin,0,Egammamin,1e4);
// //     line1->SetLineStyle(2);
// //     line1->SetLineWidth(3);
// //     line1->Draw();

// //     TLatex t;
// //     t.SetTextSize(0.05);
// //     t.DrawLatex(-130,-0.0014,"Threshold");

// //     TArrow *arrow1 = new TArrow(110,-0.00089,Egammamin,-0.0002,0.02,">");
// //     arrow1->Draw();
//     c4->Update();
    

    //////////////////////////////////////////////////




    //////////////////////////////////////////////////////

    // Create TCanvas
    TCanvas *c3 = new TCanvas("c3","Coincidence spectra",900,420);
    c3->Divide(2,1,0,0);

    c3->cd(1);
    c3_1->SetLogz();
    c3_1->SetLeftMargin(0.14);
    c3_1->SetRightMargin(0.05);
    c3_1->SetBottomMargin(0.14);
//    c1_1->SetTopMargin(0.09);
        
    matrix->GetXaxis()->SetTitle(" #gamma-ray energy E_{#gamma} (keV)");
    matrix->GetYaxis()->SetTitle(" Excitation energy E_{x} (keV)");
    matrix->GetZaxis()->SetTitleFont(42);
    matrix->GetXaxis()->SetTitleSize(0.055);
    matrix->GetZaxis()->CenterTitle();
    matrix->GetXaxis()->SetTitleOffset(1.2);
    matrix->GetYaxis()->SetTitleOffset(1.);
    matrix->GetXaxis()->SetTitleFont(42);
    matrix->GetYaxis()->SetTitleFont(42);
    matrix->GetXaxis()->SetLabelFont(42);
    matrix->GetYaxis()->SetLabelFont(42);
    matrix->GetXaxis()->SetTitleSize(0.05);
    matrix->GetYaxis()->SetTitleSize(0.055);
    matrix->GetXaxis()->SetLabelSize(0.05);
    matrix->GetYaxis()->SetLabelSize(0.05);
    matrix->GetZaxis()->SetLabelFont(42);
//    matrix->GetXaxis()->SetRangeUser(0.,12.3);
//    matrix->GetYaxis()->SetRangeUser(0.,12.3);
    matrix->GetXaxis()->SetRangeUser(0.,11e3);
    matrix->GetYaxis()->SetRangeUser(0.,11e3);
    matrix->GetZaxis()->SetRangeUser(1,3000);
    matrix->Draw("col");

    ////////////////////////////////////////////////////////////////////////////////////////////////


    ///////////////////////////
// READ THE DATA : GEF_Gamma_Spectrum
///////////////////////////
// double x,y;
// string string_x;
// char delim;

// opening the data_file
// char *myfile = "Verbinski_239Pu.csv".c_str();

ifstream thefile;
thefile.open("Theory_Christelle/GEF_gamma_spectrum_5.85MeV.txt");
int nheaderlines = 1;

// throw error if file does not exists
if(thefile.fail()){
  cout << "\n Could not open file!!!: "  << "\n ";
  exit(1);
}
else cout << "\n Successful opening of file" << endl;


int nLines=0;
//std::ifstream myfile("main.cpp");
std::string lines;
while (std::getline(thefile, lines))
       ++nLines;
thefile.close();

thefile.open("Theory_Christelle/GEF_gamma_spectrum_5.85MeV.txt");
string line; // ignore the header line
line.resize(500); //

for(int i=0;i<nheaderlines;i++){ // (for(i=0;i<2;i++) --> 3 header lines)
getline(thefile,line); // 
} // ignore the header line

nPoints = nLines - nheaderlines;

// shift_filenames.resize(nPoints);
// shift_channels.resize(nPoints);
double energy1[nPoints], spectrum1[nPoints];

  //int i = 0;
  double d1_1, d2_2;

double binSize = 100.; // in keV
double nEvents = 1000000.;
  while(thefile){
    thefile >> x >> y >> d1_1;            // readfrist columns into "x"
    if(i<nPoints){
            energy1[i]   = x*1000. + binSize/2.; // to get to keV
            spectrum1[i] = y/nEvents*10;
            // cout << energy1[i] << "\t" << spectrum1[i] << "\t" << endl;
            }
     i++;
    }
    
  

//   // for(int i=0;i<nPoints;i++)
//   // {
//   //   cout << energy[i] << "\t" << spectrum[i] << "\t" << lower_y[i] << "\t" << higher_y[i] << endl;
//   // }

// cout << "end of reading the file " << endl;

//     // TGraphAsymmErrors (Int_t n, const Double_t *x, const Double_t *y, const Double_t *exl=0, const Double_t *exh=0, const Double_t *eyl=0, const Double_t *eyh=0)
TGraph *GEF = new TGraph(nPoints,energy1,spectrum1);

//TGraph *Verbinski = new TGraph(nPoints,energy,spectrum);

//////////////////////////////////////////
///////////////////////////////////////////



//////////////////////////////////////////
// Verbinski and a spectrum at around S_n

    TCanvas *c5 = new TCanvas("c5","Verbinski vs Sn",900,600);

    c5->cd(1);
    
    c5->SetLogy();

    c5.SetLeftMargin(0.1);
    c5.SetRightMargin(0.01);
    c5.SetBottomMargin(0.15);


    int nSelect=3;

    h4[nSelect].GetXaxis().CenterTitle();
    h4[nSelect].GetXaxis().SetTitle("E_{#gamma} [keV]");
    h4[nSelect].GetYaxis().CenterTitle();
    h4[nSelect].GetXaxis().SetTitleSize(0.06);
    h4[nSelect].GetXaxis().SetTitleOffset(1.1);
    // h4[nSelect].GetYaxis().SetTitleOffset(1);
//    h2[0].GetXaxis().SetLabelOffset(1.4);

    h4[nSelect].GetYaxis().SetTitle("Photons/Fission*MeV");
    h4[nSelect].GetYaxis().SetTitleSize(0.06);
    h4[nSelect].GetYaxis().SetTitleOffset(.8);
    h4[nSelect].GetXaxis().SetLabelSize(0.053);
    h4[nSelect].GetYaxis().SetLabelSize(0.053);
    h4[nSelect].GetXaxis().SetRangeUser(0,13.8e3);
    h4[nSelect].GetYaxis().SetRangeUser(0.00017,10);
    h4[nSelect].SetLineColor(kRed);
    // h4[0].SetLineWidth(2);
    gStyle->SetOptStat(0);
    h4[nSelect].Draw();



    ////////////
    Verbinski->SetMarkerStyle(33);
    // Verbinski->SetLineStyle(2);
    Verbinski->SetMarkerSize(0.8); 
    // Verbinski->SetLineWidth(2.8); 
    Verbinski->SetMarkerColor(kBlack); 
    Verbinski->SetLineColor(kBlack);

    // Verbinski->SetMarkerColor(1);
    // Verbinski->SetMarkerSize(20); 
    Verbinski->Draw("P");
    /////////////////

    GEF->SetMarkerStyle(22);
    GEF->SetMarkerColor(kBlue+2);
    GEF->SetMarkerSize(1);
    // GEF->SetMarkerColor(kGreen);
    GEF->Draw("P");


    
    TLegend *leg = new TLegend(0.566964,0.719388  ,0.90,0.892);
    leg.SetBorderSize(0);
    leg.SetFillStyle(1);

    for (i=nSelect; i<nSelect+1; i++) {

        // xinterval = Erange / nIntervals;
        // xmin[i] = Emin + i * xinterval;
        // xmax[i] = xmin[i] + xinterval;

        sstm.str(std::string());
        sstm << " E_x = " << std::setprecision(3) << xmin[i]/1000 << " - " << xmax[i]/1000 <<" MeV";
        string_result = sstm.str();
        histname_result = string_result.c_str();
        leg->AddEntry(h4[i],histname_result,"L");
    }
    cout << "\n" << endl;
    // leg->AddEntry(h4[i],histname_result,"L")
    leg->AddEntry(Verbinski,"Verbinski","ep");
    leg->AddEntry(GEF,"GEF at Ex=XXX","p");

    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->Draw();

    c5->Print("Compare_to_Verbinski.pdf");
    // h4[nSelect]->Print("all");

  //     char out_file_name[500];
  // snprintf(out_file_name, sizeof(out_file_name), "%s.hist", h4[nSelect]->GetName());
  // FILE *out_file;
  // out_file = fopen(out_file_name, "w");
  // cout << "Making output file: " << out_file_name << endl;
  // fprintf(out_file, "%s\t%d\n",
  //       hist_name, h[nSelect] -> GetNbinsX());
  // for (int i = 1; i <= h[nSelect] -> GetNbinsX(); i++) {
  //   fprintf(out_file, "%g\t%g\n",
  //       hist -> GetBinCenter(i), hist -> GetBinContent(i));
  // }
  // fclose(out_file);
  // cout << "Output complete" << endl;

cout << "\n" << endl;
cout << "meanEnergies[i]" << "\t" << "TotalEnergy[i]" << "\t" << "stat. unc." << "\t" << "Unc.(incl. systematic)" << "\t" << "AverageEnergy[i]"       << "\t" << "stat. unc." << "\t" << "Unc.(incl. systematic)"             << "\t" << "multiplicity[i]" << "\t" << "stat. unc." << "\t" << "Unc.(incl. systematic)" << endl;
cout << "keV" << "\t" << "keV" << "\t" << "keV" << "\t" "keV" << "\t" "keV" << "\t" << "keV"       << "\t" << "keV"             << "\t" << "ph/fiss" << "\t" << "ph/fiss" << "\t" << "ph/fiss" << endl;
for(int i=0;i<nIntervals;i++){
    cout << meanEnergies[i] << "\t" << 
    TotalEnergy[i] << "\t"      << TotalEnergyErr[i] << "\t"          << TotalEnergy[i]-TotalEnergyLow[i] << "\t" 
    << AverageEnergy[i] << "\t" << AverageEnergyErr[i] << "\t"        << AverageEnergy[i]-AverageEnergyLow[i] << "\t" 
    << multiplicity[i] << "\t"  << multiplicityErr[i]  << "\t"        << multiplicity[i]-multiplicityLow[i] <<   endl;
}

}
