#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "TXMLEngine.h"
#include <string>
#include <cctype>

int DopeType;

double epsilon_0 = 0.0885418782; // [pF/cm]
double epsilon_Si = 11.68*epsilon_0; // [pF/cm]
double epsilon_SiO2 = 3.9*epsilon_0; // [pF/cm]
double MOS_area;
double q= 1.60219E-7; // [pF.V]
double kB = 8.6173324E-5 * q; // [pF*V^2/K]
double TKelvin;
double Nintrisic;


float Temp_measurement[100];
float Humid_measurement[100];

// ------------------------------------------------------------------------------------//
// MOS arrays and parameters declartion
//------------------------------------------------------------------------------------//
double MOS_A_constant[100];
double MOS_A_slope[100];
double MOS_D_constant[100];
double MOS_D_slope[100];
double MOS_I_constant[100];
double MOS_I_slope[100];
float MOS_cv_low[100];
float MOS_cv_high[100];
float MOS_cv1_low[100];
float MOS_cv1_high[100];
float MOS_cv2_low[100];
float MOS_cv2_high[100];
float Cacc[100];
float Cinv[100];
float Tox[100];
float VFB[100];
float VFB_PQC[100];
float VFB_PQC_Derivative[100];
float CFB[100];
float CFB_PQC_Derivative[100];
float Ndop[100];
float Noxide[100];
float Noxide_PQC[100];
float Noxide_PQC_Derivative[100];
float MOS_voltage,MOS_capacitance;
float MOS_vv[500][500];
float MOS_cap[500][500];
float inv_MOS_cap[500];
double MOS_vv_d[500][500];
double MOS_CAP_d[500][500];
double inv_MOS_CAP_d[500];
double MOS_dev_spline_simple[500];
double MOS_secdev_spline_simple[500];
double MOS_dev_spline_InvCV2[500];
vector<vector <string> > MOS_timestamp_meas;
vector<vector <double> > MOS_rsstnc_mhom;
vector<vector <double> > MOS_temp_meas;
vector<vector <double> > MOS_air_temp_meas;
vector<vector <double> > MOS_rh_prcnt_meas;
vector<std::string> MOS_Operators;
vector<std::string> MOS_Begin_Timestamp;
vector<std::string> MOS_name_labels;
vector<std::string> MOS_kind_of_parts;
vector<std::string> MOS_kind_of_HM_flute_id;
vector<std::string> MOS_struct_id;
vector<std::string> MOS_set_id;
vector<std::string> MOS_config_id;
vector<std::string> MOS_equipment;
vector<double> MOS_waiting_time;
vector<double> MOS_ac_freq;
vector<double> MOS_ac_ampl;
vector<double> MOS_temp;
vector<double> MOS_av_temp;
vector<double> MOS_compliance;
vector<double> MOS_Nmeas;
TGraph *VDP_cv[600];
TGraph *MOS_cv[100];
TGraph *MOS_invC2_vs_V[100];
TSpline3 * CuSpl_MOS_cv[100];
TSpline3 * CuSpl_firstDerivative_MOS_cv[100];
TSpline3 * CuSpl_MOS_InvCV2[100];
// ------------------------------------------------------------------------------------
// VDP arrays and parameters declartion
//--------------------------------------------------------------------------------
float VDP_current[500][500];
float VDP_voltage[500][500];
float VDP_current_d[500][500];
float VDP_voltage_d[500][500];
double VDP_A_constant[600];
double VDP_A_slope[600];
float VDP_cv_low[600];
float VDP_cv_high[600];
double VDPoly_s_Rsh_E[100];
double VDPoly_r_Rsh_E[100];
double VDPstop_s_Rsh_E[100];
double VDPstop_r_Rsh_E[100];
double VDPStrip_s_Rsh_E[100];
double VDPStrip_r_Rsh_E[100];
double VDPoly_s_Rsh_W[100];
double VDPoly_r_Rsh_W[100];
double VDPstop_s_Rsh_W[100];
double VDPstop_r_Rsh_W[100];
double VDPStrip_s_Rsh_W[100];
double VDPStrip_r_Rsh_W[100];
double VDP_Rsh_all[100];    // all Rsh values in one array
vector<vector <string> > VDP_timestamp_meas;
vector<vector <double> > VDP_temp_meas;
vector<vector <double> > VDP_air_temp_meas;
vector<vector <double> > VDP_rh_prcnt_meas;
vector<std::string> VDP_Operators;
vector<std::string> VDP_Begin_Timestamp;
vector<std::string> VDP_name_labels;
vector<std::string> VDP_kind_of_parts;
vector<std::string> VDP_kind_of_HM_flute_id;
vector<std::string> VDP_struct_id;
vector<std::string> VDP_set_id;
vector<std::string> VDP_config_id;
vector<std::string> VDP_equipment;
vector<double> VDP_waiting_time;
vector<double> VDP_temp;
vector<double> VDP_av_temp;
vector<double> VDP_Nmeas;
// ------------------------------------------------------------------------------------
// FET arrays and parameters declartion
//--------------------------------------------------------------------------------
float FET_Gate_Voltage_Threshold[40];
float FET_Voltage,FET_Current;
float FET_vv[500][500];
float FET_Curr[500][500];
vector<vector <string> > FET_timestamp_meas;
vector<vector <double> > FET_temp_meas;
vector<vector <double> > FET_air_temp_meas;
vector<vector <double> > FET_rh_prcnt_meas;
double FET_vv_d[500][500];
double FET_Curr_d[500][500];
double FET_IV_CubicSpline_first_Derivative[500];
double FET_IV_CubicSpline_second_Derivative[500];
double FET_IV_CubicSpline_second_Derivative_ext[2000];
float FET_Temp_measurement[500];
float FET_Humid_measurement[500];
vector<std::string> FET_Operators;
vector<std::string> FET_Begin_Timestamp;
vector<std::string> FET_name_labels;
vector<std::string> FET_kind_of_parts;
vector<std::string> FET_struct_id;
vector<std::string> FET_set_id;
vector<std::string> FET_config_id;
vector<std::string> FET_equipment;
vector<double> FET_waiting_time;
vector<double> FET_bias;
vector<double> FET_temp;
vector<double> FET_av_temp;
vector<double> FET_compliance;
vector<double> FET_Nmeas;

TGraph *FET_iv[40];
TGraph *FET_invC2_vs_V[40];
TSpline3 * FET_CuSpl_iv[40];
TSpline3 * FET_CuSpl_firstDerivative_iv[40];
// ------------------------------------------------------------------------------------
// Capacitors arrays and parameters declartion
//--------------------------------------------------------------------------------
float CAP_CAP_Capacitance[100];
float CAP_vv[500][500];
float CAP_Capacitance[500][500];
double CAP_vv_d[500][500];
double CAP_Capacitance_d[500][500];
float CAP_Temp_measurement[500];
float CAP_Humid_measurement[500];
TGraph *CAP_cv[100];
float CAP_Value[100];
float CAP_dox[100];
vector<vector <string> > CAP_timestamp_meas;
vector<vector <double> > CAP_rsstnc_mhom;
vector<vector <double> > CAP_temp_meas;
vector<vector <double> > CAP_air_temp_meas;
vector<vector <double> > CAP_rh_prcnt_meas;
vector<std::string> CAP_Operators;
vector<std::string> CAP_Begin_Timestamp;
vector<std::string> CAP_name_labels;
vector<std::string> CAP_kind_of_parts;
vector<std::string> CAP_kind_of_HM_flute_id;
vector<std::string> CAP_struct_id;
vector<std::string> CAP_set_id;
vector<std::string> CAP_config_id;
vector<std::string> CAP_equipment;
vector<double> CAP_waiting_time;
vector<double> CAP_ac_freq;
vector<double> CAP_ac_ampl;
vector<double> CAP_temp;
vector<double> CAP_av_temp;
vector<double> CAP_compliance;
vector<double> CAP_Nmeas;
//-------------------------------------------------------------------------------//
ofstream CVS_Output_File;

TString Data_Dir = "/home/akyriakis/MOS-Measurements_Analysis/PQC_Measurements/Data_New_Format/";

TString Analysis_DataFileName[100];
TString MOS_DataFileName[100];
TString VDP_DataFileName[100];
TString CAP_DataFileName[100];

//int HM_Id = 35954;
//string Structure_Id[40] = {"01","15","31","38","42"};
//int HM_Id = 35955; // Attention: from this Device ID all CV measurements save also the Conductance so we should read 4 values from the file: Voltage, Capacitance, Conductance
//string Structure_Id[40] = {"08","18","26","35","44"};
//int HM_Id = 36793;
//string Structure_Id[40] = {"03","18","29","38"};
//int HM_Id = 36794;
//string Structure_Id[40] = {"02","15","26","35"};
//int HM_Id = 37077;
//string Structure_Id[40] = {"06","19","28","37"};
//int HM_Id = 37077;
//string Structure_Id[40] = {"09","14","27","48"};
//string Structure_Id[40] = {"09"};
//int HM_Id = 37400;
//string Structure_Id[40] = {"03","15","27","37"}; 
//int HM_Id = 37401;
//string Structure_Id[40] = {"05","14","28","35"}; 
//int HM_Id = 37408;
//string Structure_Id[40] = {"04","16","27","34"}; 
int HM_Id = 37897;
string Structure_Id[40] = {"01","12","26","39"}; 


int xml_on = 1 ;  // set 1 for xml production

// This function converts a string to char
string & Conv_Char_To_String(char c[20])
{
  static string * s;
    s = new string;
    //std::string s;
    std::stringstream ss;
    ss << c;
    ss >> *s;                // or, use `s = ss.str()
    return *s;
}
// This function converts a float to string without extra zeros at the end
string & Conv_float_to_string_scientific_format(float num, int prec_num){

  static string *strObj;
  strObj = new string;
  // Create an output string stream
  std::stringstream streamObj;
  streamObj << std::scientific;
  //Add double to stream
  streamObj << std::setprecision(prec_num) << num;
  // Get string from output string stream
  streamObj >> *strObj;
  return *strObj;
}
// This function converts a float to string without extra zeros at the end
string & Conv_float_to_string(float num, int prec_num){

  static string *strObj;
  strObj = new string;
  // Create an output string stream
  std::stringstream streamObj;
  streamObj << std::fixed;
  //Add double to stream
  streamObj << std::setprecision(prec_num) << num;
  // Get string from output string stream
  streamObj >> *strObj;
  return *strObj;
}
   // define a function with 2 parameters
Double_t fitf(Double_t *x,Double_t *par) {
//   Double_t NoxFit = par[0] + par[1]*TMath::Log(x[0]);
//   Double_t NoxFit = par[0] + par[1]*x[0];
   Double_t NoxFit = par[0] + par[1]*log(x[0]);
   return NoxFit;
}
//--------------------------------------------------------------------------------------------------//
// MOS Read data
// --------------------------------------------------------------------------------------------------//
void MOS_Read_Info_from_File(TString MOS_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(MOS_DataFileName);

   DopeType = 1; // +1 for p-type Substrate, -1 for n-type Substrate
   MOS_area = 0.129*0.129; // [cm^2]   q = 1.60219E-7; // [pF.V]


   cout << endl << "=============>>>>>> Analyze File : " << MOS_DataFileName << endl;

   int nlines = 0;
   int N_meas = 0;
   string line;
   float MOS_voltage, MOS_capacitance, MOS_conductance;
   while (in.good()) {
     getline (in,line);
      if(nlines < 17) {
         if(nlines == 1) {
           cout << line << '\n';
           char First_Name[20], Last_Name[20];
           sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
           cout<<First_Name<<Last_Name<<endl;
           MOS_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
           cout << line << '\n';
           char Date[20], Time[20];
           sscanf(line.c_str(), "%*s %s %s", Date, Time);
           cout<<Date<<Time<<endl;
           MOS_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
           cout << line << '\n';
           char Name_Label[25];
           sscanf(line.c_str(), "%*s %s", Name_Label);
           cout<<Name_Label<<endl;
           MOS_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
           cout << line << '\n';
           char batch_type[20], loc[20];
           sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
           string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
           cout<<Kind_of_part<<endl;
           MOS_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
           cout << line << '\n';
           char Kind_of_HM_flute_id[20];
           sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
           cout<<Kind_of_HM_flute_id<<endl;
           MOS_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
           cout << line << '\n';
           char Kind_of_HM_struct_id[20];
           sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
           cout<<Kind_of_HM_struct_id<<endl;
           MOS_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
           cout << line << '\n';
           char Kind_of_HM_set_id[20];
           sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
           cout<<Kind_of_HM_set_id<<endl;
           MOS_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
           cout << line << '\n';
           char str1[20], str2[20];
           sscanf(line.c_str(), "%*s %s %s", str1, str2);
           string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
           cout<<Kind_of_HM_config_id<<endl;
           MOS_config_id.push_back(Kind_of_HM_config_id);
         }
         if(nlines == 9) {
           cout << line << '\n';
           char Procedure_type[20];
           sscanf(line.c_str(), "%*s %s", Procedure_type);
           cout<<Procedure_type<<endl;
         }
         if(nlines == 10) {
           cout << line << '\n';
           double Waiting_time_s;
           sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
           cout<<Waiting_time_s<<endl;
           MOS_waiting_time.push_back(Waiting_time_s);
         }
         if(nlines == 11) {
           cout << line << '\n';
           double ac_freq;
           sscanf(line.c_str(), "%*s %lf", &ac_freq);
           cout<<ac_freq<<endl;
           MOS_ac_freq.push_back(ac_freq);
         }
         if(nlines == 12) {
           cout << line << '\n';
           double ac_ampl;
           sscanf(line.c_str(), "%*s %lf", &ac_ampl);
           cout<<ac_ampl<<endl;
           MOS_ac_ampl.push_back(ac_ampl);
        }
        if(nlines == 13) {
          cout << line << '\n';
          double temp_degC;
          sscanf(line.c_str(), "%*s %lf", &temp_degC);
          cout<<temp_degC<<endl;
          MOS_temp.push_back(temp_degC);
          TKelvin = 273 + temp_degC; // [K]
        }
        if(nlines == 14) {
          cout << line << '\n';
          double Av_temp_degC;
          sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
          cout<<Av_temp_degC<<endl;
          MOS_av_temp.push_back(Av_temp_degC);
        }
      }
      else {

      Nintrisic= (5.29*1E+19)*(TMath::Power((TKelvin/300), 2.5))*exp(-6726/TKelvin); // [cm^-3]
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double capctnc_pfrd, rsstnc_mhom;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf %lf", date, time, &voltage, &capctnc_pfrd, &rsstnc_mhom, &temp_degC, &air_temp_degC, &rh_prcnt);
      // /cout<<date<<","<<time<<","<<voltage<<","<<(float) capctnc_pfrd<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
      //        in >>  MOS_voltage >> MOS_capacitance >> MOS_conductance;
      //  in >>  MOS_voltage >> MOS_capacitance;
    	MOS_vv_d[id][N_meas] =  -voltage;
    	MOS_CAP_d[id][N_meas] = capctnc_pfrd;

      MOS_temp_meas.push_back(std::vector<double>());
      MOS_temp_meas[id].push_back(temp_degC);
      MOS_air_temp_meas.push_back(std::vector<double>());
      MOS_air_temp_meas[id].push_back(air_temp_degC);
      MOS_rh_prcnt_meas.push_back(std::vector<double>());
      MOS_rh_prcnt_meas[id].push_back(rh_prcnt);
      MOS_timestamp_meas.push_back(std::vector<string>());
      MOS_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      MOS_rsstnc_mhom.push_back(std::vector<double>());
      MOS_rsstnc_mhom[id].push_back(rsstnc_mhom);

    	MOS_vv[id][N_meas] = float(MOS_vv_d[id][N_meas]);
      MOS_cap[id][N_meas] = float(MOS_CAP_d[id][N_meas]);

      inv_MOS_cap[N_meas] = 1/(TMath::Power(MOS_cap[id][N_meas],2));
      inv_MOS_CAP_d[N_meas] = 1/(TMath::Power(MOS_CAP_d[id][N_meas],2));
      //        cout << N_meas << ")  "  <<  MOS_vv[N_meas] << "   " << MOS_cap[N_meas] << endl;
      //cout << N_meas << ")  "  <<  MOS_vv[N_meas] << "   " << MOS_cap[N_meas] << "   " << inv_MOS_cap[N_meas] << endl;
    	N_meas++;
      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   cout  << "----->>>>>>   Number of Measurements = " << N_meas << endl;
   MOS_Nmeas.push_back(N_meas); // insert the number of measurments to an array
   //---------- Cubic spline interpolation to the data---------------------------------------------------------------//
   int  N_meas_reduced = N_meas-2;
   if(N_meas_reduced > 300 ) {
      cout << "Serious Problem : Measurements > 300 that TSpline3 can handle...  " << endl;
      cout << " Solution: Reduce number of measurements for TSpline3 rootine.... " << endl;
      cout << " Now Abort..." << endl;
      abort();
   }
   double vv_d_reduced[300];
   double CAP_d_reduced[300];
   double inv_CAP_d_reduced[300];
   for (Int_t i=0; i<N_meas_reduced; i++) {
      vv_d_reduced[i] = MOS_vv_d[id][N_meas_reduced-i]; // reverse the order for min to max to work the TSpline
      CAP_d_reduced[i] = MOS_CAP_d[id][N_meas_reduced-i];
      inv_CAP_d_reduced[i] = inv_MOS_CAP_d[N_meas_reduced-i];
      //cout << i << ")  "  <<  vv_d_reduced[i] << "   " << CAP_d_reduced[i] << endl;
   }
   // TSpline workes withnordered array from Min to Max and not vice versa
   // Also Seems that works well for less than 150 points
   CuSpl_MOS_cv[id] = new TSpline3("Cubic Spline", vv_d_reduced, CAP_d_reduced , N_meas_reduced, "cb2e2", 0, 0);
   //-------------------------------------------------------------------------------------------------------------//
  //   TCanvas *cc_spline = new TCanvas();
   MOS_cv[id] =  new TGraph(N_meas_reduced,vv_d_reduced,CAP_d_reduced);
   if((id % 2) == 0) MOS_cv[id]->SetTitle(Form("MOS%d_0%s CV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) MOS_cv[id]->SetTitle(Form("MOS%d_0%s CV",HM_Id,Structure_Id[file_id].c_str()));
   MOS_cv[id]->GetXaxis()->SetTitle("Gate MOS_voltage [V]");
   MOS_cv[id]->GetYaxis()->SetTitle("MOS MOS_capacitance [pF]");
   MOS_cv[id]->SetDrawOption("AP");
   MOS_cv[id]->SetMarkerStyle(20+id);
   int colorid = 1+id;
   if(colorid == 5) colorid = 46;
   MOS_cv[id]->SetMarkerColor(colorid);
    //   MOS_cv[id]->Draw();
    //   CuSpl_MOS_cv[id]->SetLineColor(kBlue);
    //   CuSpl_MOS_cv[id]->Draw("same");

    //   for (Int_t i=0; i<N_meas_reduced; i++) {
    //
    //      cout << "Spline : " << i << ")  "  <<  MOS_vv_d[i] << "   " << CuSpl_MOS_cv[id]->Eval(MOS_vv_d[i]) << endl;
    //
    //   }
    // Get V_Fb grosso - modo
    // --------------------------

    // Get the First Derivative of CV curve
    for(int i = 0;i<N_meas_reduced; i++) {
       MOS_dev_spline_simple[i] = CuSpl_MOS_cv[id]->Derivative(vv_d_reduced[i]);
    }
    // Make Cubic Spline fit of the First Derivative
    CuSpl_firstDerivative_MOS_cv[id] = new TSpline3("Cubic Spline", vv_d_reduced, MOS_dev_spline_simple , N_meas_reduced, "cb2e2", 0, 0);
    //    for (Int_t i=0; i<N_meas; i++) {
    //
    //      cout << i << ")  "  <<  MOS_vv_d[i] << "   " << CuSpl_firstDerivative_MOS_cv[id]->Eval(vv_d_reduced[i]) << endl;
    //
    //    }
    // get the minimum of the derivative
    int min_index_simple = 0;
    float min_MOS_dev_spline_simple = MOS_dev_spline_simple[0];
    for (int i = 1;i<N_meas; i++) {
       if((MOS_dev_spline_simple[i] <  min_MOS_dev_spline_simple) && vv_d_reduced[i] <0) {
        //  cout<< min_index_simple<<endl;
          min_MOS_dev_spline_simple = MOS_dev_spline_simple[i];
	        min_index_simple = i;
       }
    }
    double V_FB_derivative_simple = double(MOS_vv_d[id][min_index_simple]);
    double C_FB_derivative_simple =CuSpl_MOS_cv[id]->Eval(V_FB_derivative_simple);
    cout  << "---------------- >>>>> Min derivative Simple at i = " << min_index_simple
          << " with value = " << MOS_dev_spline_simple[min_index_simple]
	        << " V_flatband_Simple = " <<  V_FB_derivative_simple
	        << " and C_flatband_Simple = " <<  C_FB_derivative_simple<< endl;
    // Get V_Fb fine tuned
    // -------------------------
    int Nub_steps = 1000;
    float vv_initial = vv_d_reduced[0];
    float vv_final = vv_d_reduced[N_meas_reduced-1];
    Double_t dev_spline[2000];
    Float_t dev_spline_f[2000];
    Float_t MOS_voltage_extented[2000];
    for(int i = 0;i<Nub_steps; i++) {
      MOS_voltage_extented[i] =  vv_initial + i*(vv_final - vv_initial)/(Nub_steps-1);
      dev_spline[i] = CuSpl_MOS_cv[id]->Derivative(MOS_voltage_extented[i]);
      dev_spline_f[i] = dev_spline[i];
    }
    //    TCanvas *cc_spline_der = new TCanvas();
    TGraph *cv_deriv =  new TGraph(Nub_steps,MOS_voltage_extented,dev_spline_f);
    cv_deriv->SetTitle(Form("First Derivative of MOS CV"));
    cv_deriv->GetXaxis()->SetTitle("Gate MOS_voltage [V]");
    //    cv_deriv->Draw();
    CuSpl_firstDerivative_MOS_cv[id]->SetLineColor(kBlue);
    //    CuSpl_firstDerivative_MOS_cv[id]->Draw("same");
    Double_t secdev_spline[2000];
    Float_t secdev_spline_f[2000];
    for(int i = 0;i<Nub_steps; i++) {
      MOS_voltage_extented[i] =  vv_initial + i*(vv_final - vv_initial)/(Nub_steps-1);
      secdev_spline[i] = CuSpl_firstDerivative_MOS_cv[id]->Derivative(MOS_voltage_extented[i]);
      secdev_spline_f[i] = secdev_spline[i];
    }
    //    TCanvas *cc_spline_secder = new TCanvas();
    TGraph *cv_secderiv =  new TGraph(Nub_steps,MOS_voltage_extented,secdev_spline_f);
    cv_secderiv->SetTitle(Form("Second Derivative of MOS CV"));
    cv_secderiv->GetXaxis()->SetTitle("Gate MOS_voltage [V]");
    //    cv_secderiv->Draw();
    // get the minimum of the derivative
    int min_index = 0;
    float min_dev_spline_f = dev_spline_f[0];
    for (int i = 1;i<Nub_steps; i++) {
       if((dev_spline_f[i] < min_dev_spline_f) &&  MOS_voltage_extented[i]<0) {
         min_dev_spline_f = dev_spline_f[i];
	 min_index = i;
       }
    }
    double V_FB_derivative = double(MOS_voltage_extented[min_index]);
    double C_FB_derivative =CuSpl_MOS_cv[id]->Eval(V_FB_derivative);
    cout  << "---------------- >>>>> Min derivative accurate at i = " << min_index
          << "  with value = " << dev_spline_f[min_index]
	        << " V_flatband = " <<  V_FB_derivative
	        << " and C_flatband = " << C_FB_derivative << endl;
   // Set linear fit ranges
   // --------------------------------------------------------------------------
   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------------------------------
   MOS_cv_low[id] = vv_d_reduced[0];
   MOS_cv_high[id] = vv_d_reduced[40];
   // Depletion Region Linear Fit
   // --------------------------------------------------------------------------
   MOS_cv1_low[id] = vv_d_reduced[min_index_simple-3];
   MOS_cv1_high[id] = vv_d_reduced[min_index_simple+3];
   if(id>0) {
     MOS_cv1_high[id] = vv_d_reduced[min_index_simple+6];
   }
   // ---------------------------------------------------
   // Inversion Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   MOS_cv2_low[id] = vv_d_reduced[N_meas_reduced-20];;
   MOS_cv2_high[id] = vv_d_reduced[N_meas_reduced-1];
   // Get tree branches to find the fit points
   // -------------------------------------------
   //   Int_t nentries = (Int_t)T->GetEntries();
   //   cout << "Number of Measurements = " << nentries << endl;
   // --------------------------------------------------------------------------
   TF1 *f_cv = new TF1("f_cv","pol0", MOS_cv_low[id], MOS_cv_high[id]);
   //   f_cv->SetParameter(1,0.0);
   MOS_cv[id]->Fit("f_cv","0R+");
   TF1 *f_cv1 = new TF1("f_cv1","pol1",MOS_cv1_low[id], MOS_cv1_high[id]);
   MOS_cv[id]->Fit("f_cv1","0R+");
   TF1 *f_cv2 = new TF1("f_cv2","pol0", MOS_cv2_low[id], MOS_cv2_high[id]);
   MOS_cv[id]->Fit("f_cv2","0R+");

   // Accumulation Region
   // ---------------------
   MOS_A_constant[id] = MOS_cv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   MOS_A_slope[id]  = 0.0; // Accumulation line slope
   double C_ox = MOS_A_constant[id]; // MOS Accumulation MOS_capacitance
   double C_ox_Error = MOS_cv[id]->GetFunction("f_cv")->GetParError(0);
   cout << "----->>>   Accumulation MOS_capacitance = " << C_ox << " +/- " << C_ox_Error << " [pF]" << endl;
   double MOS_Oxide_Thickness = epsilon_SiO2 * MOS_area/C_ox;
   cout << "----->>>   Oxide Layer Thickness = " << MOS_Oxide_Thickness
        << " [cm]" << " = " << MOS_Oxide_Thickness*1E4 << " [um]" << endl;

   // Depletion Region
   // ---------------------
   MOS_D_constant[id] = MOS_cv[id]->GetFunction("f_cv1")->GetParameter(0); // Depletion Region constant term
   MOS_D_slope[id] = MOS_cv[id]->GetFunction("f_cv1")->GetParameter(1); // Depletion Region slope

   // Inversion Region
   // ---------------------
   MOS_I_constant[id] = MOS_cv[id]->GetFunction("f_cv2")->GetParameter(0); // Inversion constant term
   MOS_I_slope[id]  = 0.0; // Inversion line slope
   double C_inv = MOS_I_constant[id]; // MOS High Frequency MOS_capacitance - Inversion MOS_capacitance
   double C_inv_Error = MOS_cv[id]->GetFunction("f_cv2")->GetParError(0);
   cout << "----->>>   Inversion MOS_capacitance = " << C_inv << " +/- " << C_inv_Error << " [pF]" << endl;
   // ------------------------------------------------------------------------------------------------------------------------------------

   MOS_invC2_vs_V[id] =  new TGraph(N_meas,MOS_vv[id],inv_MOS_cap);
   MOS_invC2_vs_V[id]->SetTitle("MOS Flute 1");
   MOS_invC2_vs_V[id]->SetDrawOption("AP");
   MOS_invC2_vs_V[id]->SetMarkerStyle(20+id);
   int colorid1 = 1+id;
   if(colorid1 == 5) colorid1 = 46;
   MOS_invC2_vs_V[id]->SetMarkerColor(colorid1);

   CuSpl_MOS_InvCV2[id] = new TSpline3("Cubic Spline", vv_d_reduced, inv_CAP_d_reduced, N_meas_reduced, "cb2e2", 0, 0);

    // Find the Best Fit point for 1/CV^2
    // --------------------------------
    for(int i = 0;i<N_meas; i++) {
       MOS_dev_spline_InvCV2[i] = CuSpl_MOS_InvCV2[id]->Derivative(MOS_vv[id][i]);
    }
    // get the minimum of the derivative
    int max_index_InvCV2 = 0;
    float max_MOS_dev_spline_InvCV2 = MOS_dev_spline_InvCV2[0];
    for (int i = 1;i<N_meas; i++) {
       if(MOS_dev_spline_InvCV2[i] >  max_MOS_dev_spline_InvCV2) {
         max_MOS_dev_spline_InvCV2 = MOS_dev_spline_InvCV2[i];
	 max_index_InvCV2 = i;
       }
    }
    cout  << "---------------- >>>>> Max derivative  InvCV2at i = " << max_index_InvCV2
          << " with InvCV2 value = " << inv_MOS_cap[max_index_InvCV2] << endl;

   float invc2v_low;
   float invc2v_high;

   invc2v_low = MOS_vv_d[id][max_index_InvCV2 - 2];
   invc2v_high = MOS_vv_d[id][max_index_InvCV2 + 2];
   if(id>0) {
     invc2v_low = MOS_vv_d[id][max_index_InvCV2 - 4];
     invc2v_high = MOS_vv_d[id][max_index_InvCV2 + 4];
   }
   cout << "----->>> invc2v_low = " << invc2v_low << " invc2v_high = " << invc2v_high <<  endl;


   TF1 *f_invc2v = new TF1("f_invc2v","pol1",invc2v_low , invc2v_high);

   MOS_invC2_vs_V[id]->Fit("f_invc2v","0R+");

   double invc2v_fit_slope = MOS_invC2_vs_V[id]->GetFunction("f_invc2v")->GetParameter(1);
   double invc2v_fit_slopeError = MOS_invC2_vs_V[id]->GetFunction("f_invc2v")->GetParError(1);
   cout << "d(1/C^2)/dV = " << invc2v_fit_slope << " +/- " << invc2v_fit_slopeError << " [pF^-2/V]" << endl;
   double Nsub = 2/(q*epsilon_Si*MOS_area*MOS_area*invc2v_fit_slope);
//   cout << "----->>>   Substrate Dopping Concentration = " << Nsub << " [cm-3]" << endl;
   cout << "kB = " <<  kB << endl;
   cout << "kB * T = " <<  kB * TKelvin << endl;
   double L_Debey = sqrt(epsilon_Si*kB * TKelvin/(q*q*Nsub));
   double C_Debey = epsilon_Si*MOS_area/L_Debey;
//   cout << "----->>>   L_Debey = " << L_Debey << " [cm]" << " = " << L_Debey*1E4 << " [um]"<< endl;
   double C_FB_Term = epsilon_Si*MOS_area/L_Debey;
   double C_FB = C_ox*C_FB_Term/(C_ox + C_FB_Term );
//   cout << "----->>>   Flatband MOS_capacitance = " << C_FB << " [pF]" << endl;
   double C_ox1 = MOS_cv[id]->GetFunction("f_cv1")->GetParameter(0);
   double cv_fit_slope1 = MOS_cv[id]->GetFunction("f_cv1")->GetParameter(1);
   double V_FB = (C_FB - C_ox1)/cv_fit_slope1;
   double Bulk_Potential  =  - kB * TKelvin * log(Nsub/Nintrisic)*DopeType/q; // [V]
   double MS_Work_Function = 4.1 -(4.15 + 1.12/2) + Bulk_Potential; // [V] 4.1 -> Al work function, 4.15 -> Si work function, 1.12 -> Si Energy gap
   double Threshold_MOS_voltage  = V_FB + 2*fabs(Bulk_Potential) + MOS_area*sqrt(4*epsilon_Si*q*Nsub*fabs(Bulk_Potential))/C_ox;
   double Eff_Oxide_Charge = C_ox*(MS_Work_Function - V_FB)/MOS_area;
   double W_depletion = MOS_area*epsilon_Si*(1/C_inv - 1/C_ox); // [cm]
   double V_FlatBand_PQC = (MOS_A_constant[id] - MOS_D_constant[id])/(MOS_D_slope[id]  - MOS_A_slope[id]);
   double N_oxide = C_ox*(MS_Work_Function - V_FB)/(q*MOS_area);
   double Nox_PQC = C_ox*(MS_Work_Function - V_FlatBand_PQC)/(q*MOS_area);
   double Nox_PQC_Derivative = C_ox*(MS_Work_Function - V_FB_derivative)/(q*MOS_area);

   cout << "----->>>   Accumulation MOS_capacitance = " << C_ox << " +/- " << C_ox_Error << " [pF]" << endl;
   cout << "----->>>   High Frequency MOS_capacitance = " << C_inv << " +/- " << C_inv_Error << " [pF]" << endl;
   cout << "----->>>   Oxide Layer Thickness = " << MOS_Oxide_Thickness << " [cm]" << " = " << MOS_Oxide_Thickness*1E7 << " [nm]" << endl;
   cout << "----->>>   Substrate Dopping Concentration = " << Nsub << " [cm^{-3}]" << endl;
   cout << "----->>>   Oxide Charge Concentration = " << Nox_PQC << " [cm^{-3}]" << endl;
   cout << "----->>>   L_Debey = " << L_Debey << " [cm]" << " = " << L_Debey*1E7 << " [nm]"<< endl;
   cout << "----->>>   C_Debey = " << C_Debey <<  " [pF]" << endl;
   cout << "----->>>   Flatband MOS_capacitance = " << C_FB << " [pF]" << endl;
   cout << "----->>>   Flatband MOS_capacitance from derivative method= " << C_FB_derivative << " [pF]" << endl;
   cout << "----->>>   Flatband MOS_voltage = " << V_FB << " [V]" << endl;
   cout << "----->>>   Flatband MOS_voltage from derivative method = " << V_FB_derivative << " [V]" << endl;
   cout << "----->>>   Flatband MOS_voltage from the intersection of the two Accumulation - Depletion lines = " << V_FlatBand_PQC << " [V]" << endl;
   cout << "----->>>   Bulk Potential (Phi_B) = " << Bulk_Potential << " [V]" << endl;
   cout << "----->>>   Metal Semiconductor Work Function (W_MS) = " << MS_Work_Function << " [V]" << endl;
   cout << "----->>>   Threshold MOS_voltage  = " << Threshold_MOS_voltage << " [V]" << endl;
   cout << "----->>>   Effective Oxide Charge   = " << Eff_Oxide_Charge << " [pF.V]  = " <<  Eff_Oxide_Charge*1E-12 << " [Cb]"<< endl;
   cout << "----->>>   Depletion Depth   = " << W_depletion << " [cm]  = " <<  W_depletion*1E4 << " [um]"<< endl;

   Cacc[id] = C_ox; // in pF
   Cinv[id] = C_inv; // in pFMOS_Flute1_MOSt_Sample_Analysis_Irradiated_Aug2020_EEAE_PQC_TSpline3_v6.C
   Tox[id] = MOS_Oxide_Thickness*1E4; // in um
   VFB[id] = V_FB; // in V
   VFB_PQC[id] = V_FlatBand_PQC; // in V
   VFB_PQC_Derivative[id] = V_FB_derivative; // in V
   CFB[id] = C_FB; // in pf
   CFB_PQC_Derivative[id] = C_FB_derivative; // in pf
   Ndop[id] = Nsub; // cm^{-2}
   Noxide[id] = N_oxide; // cm^{-3}
   Noxide_PQC[id] = Nox_PQC; // cm^{-2}
   Noxide_PQC_Derivative[id] = Nox_PQC_Derivative; // cm^{-2}


}
//--------------------------------------------------------------------------------------------------//
// MOS Abalysis Final
// --------------------------------------------------------------------------------------------------//
void MOS_PQC_Analysis_Final(int nFiles)
{
 //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1
 //   gROOT->Reset();
   gStyle->SetOptStat(0);
   for(int i=0;i<nFiles;i++) {
      MOS_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/CV/VPX%d_0%s_2-S_HM_E_Left_PQC1_MOSQuarter_CV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
      MOS_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/CV/VPX%d_0%s_2-S_HM_W_Left_PQC1_MOSQuarter_CV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
      cout << "Analyze file : " << MOS_DataFileName[i] << endl;
   }
   for (int i=0;i<nFiles;i++) {
     cout << endl << endl << "MOS Structure: Analyze file : " << MOS_DataFileName[2*i] << endl;
     MOS_Read_Info_from_File(MOS_DataFileName[2*i],2*i,i);
     cout << endl << endl << "MOS Structure: Analyze file : " << MOS_DataFileName[2*i+1] << endl;
     MOS_Read_Info_from_File(MOS_DataFileName[2*i+1],2*i+1,i);
   }

   TCanvas *MOS_cc_final[2];
   TMultiGraph *MOS_mg[2];
   for (int i=0;i<2;i++) {
     MOS_cc_final[i] = new TCanvas(Form("MOS_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();
     MOS_mg[i] = new TMultiGraph(Form("MOS_mg_%d",i),"MOS CV");
     if(i == 0) MOS_mg[i]->SetTitle("MOS_E CV ");
     if(i == 1) MOS_mg[i]->SetTitle("MOS_W CV ");
     for(int id=0;id<nFiles;id++) {
        MOS_mg[i]->Add(MOS_cv[2*id+i]);
     }

     MOS_mg[i]->Draw("LPsame");
     MOS_mg[i]->GetXaxis()->SetTitle("Gate MOS_voltage [V]");
     MOS_mg[i]->GetYaxis()->SetTitle("MOS MOS_capacitance [pF]");
     MOS_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     MOS_cc_final[i]->BuildLegend(0.6,0.6,0.85,0.85);
     gPad->Modified();
     gPad->Update();
   }
/*
   for(int id=nFiles-2;id>-1;id--) {
       TF1 *func_A = new TF1("func_A","[0] + [1]*x",MOS_cv_low[id],MOS_cv1_low[id]);
       func_A->SetParameters(MOS_A_constant[id], MOS_A_slope[id]);
       func_A->Draw("same");

       TF1 *func_D = new TF1("func_D","[0] + [1]*x",MOS_cv_high[id],MOS_cv1_high[id]);
       func_D->SetParameters(MOS_D_constant[id], MOS_D_slope[id]);
       func_D->Draw("same");

       TF1 *func_I = new TF1("func_I","[0] + [1]*x",MOS_cv1_high[id],MOS_cv2_high[id]);
       func_I->SetParameters(MOS_I_constant[id], MOS_I_slope[id]);
       func_I->Draw("same");
   }
*/
//   cc_final->Print(Form("CV_Flute1_MOS_SN_%d_250mV_%d_kHz.png",Structure_Id,freq),"png");
//   cc_final->Print(Form("CV_Flute1_MOS_SN_%d_250mV_%d_kHz.pdf",Structure_Id,freq),"pdf");
//   cc_final->Print(Form("CV_Flute1_MOS_SN_%d_250mV_%d_kHz.eps",Structure_Id,freq),"eps");
/*
   auto *cc_final_inv = new TCanvas("cc_final_inv","cc_final_inv",700,50,600,600);
   gPad->SetGrid();

   auto mg_invC2 = new TMultiGraph("mg_invC2","mg_invC2");

   for(int id=nFiles-2;id>-1;id--) {

      mg_invC2->Add(MOS_invC2_vs_V[id]);

   }
   MOS_invC2_vs_V[nFiles-1]->Draw("ALP");
   //mg_invC2->GetHistogram()->GetXaxis()->SetRangeUser(-110.,10);
   mg_invC2->SetTitle("MOS CV vs Radiation Dose");
   mg_invC2->Draw("LP");
   cc_final_inv->BuildLegend(0.6,0.6,0.85,0.85);
   gPad->Modified();
   gPad->Update();

//   cc_final_inv->Print(Form("1overC2_vs_V_Flute1_MOS_SN_%d_250mV_%d_kHz.png",Structure_Id,freq),"png");
//   cc_final_inv->Print(Form("1overC2_vs_V_Flute1_MOS_SN_%d_250mV_%d_kHz.pdf",Structure_Id,freq),"pdf");
//   cc_final_inv->Print(Form("1overC2_vs_V_Flute1_MOS_SN_%d_250mV_%d_kHz.eps",Structure_Id,freq),"eps");
*/
   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;

   cout << "        MOS Id                 | C_acc[pF]     | C_inv[pF]  |  t_ox [um] | N_Bulk[cm^{-3}] | VFB_Deb [V] | VFB_Lines [V] | VFB_Deriv [V] | Nox_Deb[cm^{-2}] | Nox_Lines[cm^{-2}] | Nox_Deriv[cm^{-2}]" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	      << setw(15) <<  setprecision(4) << Cacc[2*id]
	      << setw(15) <<  setprecision(4) << Cinv[2*id]
	      << setw(12) <<  setprecision(2) << Tox[2*id]
	      << setw(17) <<  setprecision(3)<< Ndop[2*id]
	      << setw(14) <<  VFB[2*id]
	      << setw(15) <<  VFB_PQC[2*id] << setw(17)  << VFB_PQC_Derivative[2*id]
	      << setw(20) <<  Noxide[2*id]  <<  setw(20) << Noxide_PQC[2*id] <<  setw(20)  << Noxide_PQC_Derivative[2*id]
	      << endl;

	 cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(15) <<  setprecision(4) << Cacc[2*id+1]
	       << setw(15) <<  setprecision(4) << Cinv[2*id+1]
	       << setw(12) <<  setprecision(2) << Tox[2*id+1]
	       << setw(17) <<  setprecision(3) << Ndop[2*id+1]
	       << setw(14) <<  VFB[2*id+1]
	       << setw(15) <<  VFB_PQC[2*id+1] << setw(17)  << VFB_PQC_Derivative[2*id+1]
	       << setw(20) <<  Noxide[2*id+1]  <<  setw(20) << Noxide_PQC[2*id+1] <<  setw(20)  << Noxide_PQC_Derivative[2*id+1]
	       << endl;

   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------------------------------------------//
// MOS XML production (set xml_on to 1)
// --------------------------------------------------------------------------------------------------//
void MOS_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      cout << "MOS Xml production is on" << endl;
      for(int j=0; j<=1; j++){
 //       TString xml_filename[nFiles];
        char delimeter1 = '/';
        int pos1 = MOS_DataFileName[2*i+j].Last(delimeter1);
        char delimeter2 = '.';
        int pos2 = MOS_DataFileName[2*i+j].Last(delimeter2);
        xml_filename = MOS_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1))+".xml";
//        cout << MOS_DataFileName[2*i+j](0,pos1) << endl;

        cout << xml_filename << endl;
        ofstream xml_file;

        xml_file.open(Form("./XML_Info_Flute1/VPX%d/",HM_Id)+xml_filename);

        // First create engine
        TXMLEngine MOS_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = MOS_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = MOS_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = MOS_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = MOS_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = MOS_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = MOS_xml.NewChild(HEADER, 0, "RUN");
              MOS_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              MOS_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              MOS_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              MOS_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              MOS_xml.NewChild(RUN, 0, "INITIATED_BY_USER", MOS_Operators[2*i+j].c_str());
              MOS_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", MOS_Begin_Timestamp[2*i+j].c_str());
              MOS_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = MOS_xml.NewChild(ROOT, 0, "DATA_SET");
            MOS_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            MOS_xml.NewChild(DATA_SET, 0, "VERSION", "v2");
            XMLNodePointer_t PART = MOS_xml.NewChild(DATA_SET, 0, "PART");
              MOS_xml.NewChild(PART, 0, "NAME_LABEL", MOS_name_labels[2*i+j].c_str());
              MOS_xml.NewChild(PART, 0, "KIND_OF_PART", MOS_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = MOS_xml.NewChild(DATA_SET, 0, "DATA");
              MOS_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC1");
              MOS_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", MOS_struct_id[2*i+j].c_str());
              MOS_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", MOS_set_id[2*i+j].c_str());
              MOS_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", MOS_config_id[2*i+j].c_str());
              MOS_xml.NewChild(DATA, 0, "EQUIPMENT", "HP 4192A");
              MOS_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(MOS_waiting_time[2*i+j],3).c_str());
              MOS_xml.NewChild(DATA, 0, "AC_FREQ_HZ", Conv_float_to_string_scientific_format(MOS_ac_freq[2*i+j],3).c_str());
              MOS_xml.NewChild(DATA, 0, "AC_AMPL_V", Conv_float_to_string(MOS_ac_ampl[2*i+j],3).c_str());
              MOS_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(MOS_temp[2*i+j],3).c_str());
              MOS_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(MOS_av_temp[2*i+j],3).c_str());
              MOS_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "MOS1");
//              MOS_xml.NewChild(DATA, 0, "FILE_NAME", xml_filename[2*i+j]);
              MOS_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = MOS_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = MOS_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = MOS_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  MOS_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_CV");
                  MOS_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon CV Test");
              XMLNodePointer_t DATA_SET_CDS1 = MOS_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                MOS_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                MOS_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "V2");
                XMLNodePointer_t PART_CDS1 = MOS_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  MOS_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", MOS_name_labels[2*i+j].c_str());
                  MOS_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", MOS_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<MOS_Nmeas[2*i+j]-1; k++){
                XMLNodePointer_t DATA_CDS1 = MOS_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  MOS_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(MOS_vv[2*i+j][k],3).c_str());
                  MOS_xml.NewChild(DATA_CDS1, 0, "CAPCTNC_PFRD", Conv_float_to_string(MOS_cap[2*i+j][k],3).c_str());
                  MOS_xml.NewChild(DATA_CDS1, 0, "RESSTNC_MOHM", Conv_float_to_string(MOS_rsstnc_mhom[2*i+j][k],3).c_str());
                  MOS_xml.NewChild(DATA_CDS1, 0, "TIME", MOS_timestamp_meas[2*i+j][k].c_str());
                  MOS_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(MOS_temp_meas[2*i+j][k],3).c_str());
                  MOS_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(MOS_air_temp_meas[2*i+j][k],3).c_str());
                  MOS_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(MOS_rh_prcnt_meas[2*i+j][k],3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = MOS_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = MOS_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = MOS_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      MOS_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_CV_PAR");
                      MOS_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon CV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = MOS_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    MOS_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    MOS_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "CV_measurement-004");
                    XMLNodePointer_t PART_CDS2 = MOS_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      MOS_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", MOS_name_labels[2*i+j].c_str());
                      MOS_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", MOS_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = MOS_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      MOS_xml.NewChild(DATA_CDS2, 0, "VFB_V", Conv_float_to_string(VFB_PQC[2*i+j],3).c_str());
                      MOS_xml.NewChild(DATA_CDS2, 0, "CACC_PFRD", Conv_float_to_string(Cacc[2*i+j],3).c_str());
                      MOS_xml.NewChild(DATA_CDS2, 0, "TOX_NM", Conv_float_to_string(Tox[2*i+j]*1E+3,3).c_str());
                      MOS_xml.NewChild(DATA_CDS2, 0, "NOX", Conv_float_to_string_scientific_format(Noxide[2*i+j],3).c_str());

        XMLDocPointer_t MOS_xmldoc = MOS_xml.NewDoc();
        MOS_xml.DocSetRootElement(MOS_xmldoc, ROOT);
        // Save document to file
        MOS_xml.SaveDoc(MOS_xmldoc, Form("./XML_Info_Flute1/VPX%d/",HM_Id)+xml_filename);
        // Release memory before exit
        MOS_xml.FreeDoc(MOS_xmldoc);
      }
    }
  }
}
// ----
// =====================================================================================================================================
//--------------------------------------------------------------------------------------------------//
// VDP Read Data
// --------------------------------------------------------------------------------------------------//
double VDP_Read_Info_from_File(TString VDP_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(VDP_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << VDP_DataFileName << " with Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;

   cout << in.good() << endl;
   while (in.good()) {
     getline (in,line);
      if(nlines < 16) {
         if(nlines == 1) {
         cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         VDP_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         VDP_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         VDP_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         VDP_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         VDP_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         VDP_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         VDP_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20];
         sscanf(line.c_str(), "%*s %s", str1);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1);
         cout<<Kind_of_HM_config_id<<endl;
         VDP_config_id.push_back(Kind_of_HM_config_id);
         }
         if(nlines == 9) {
         cout << line << '\n';
         char Procedure_type[20];
         sscanf(line.c_str(), "%*s %s", Procedure_type);
         cout<<Procedure_type<<endl;
         }
         if(nlines == 10) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string  Equipment =  Conv_Char_To_String(str1) +" "+Conv_Char_To_String(str2);
         cout<<Equipment<<endl;
         VDP_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         VDP_waiting_time.push_back(Waiting_time_s);
         }
        if(nlines == 12) {
        cout << line << '\n';
        double temp_degC;
        sscanf(line.c_str(), "%*s %lf", &temp_degC);
        cout<<temp_degC<<endl;
        VDP_temp.push_back(temp_degC);
       }
        if(nlines == 13) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        VDP_av_temp.push_back(Av_temp_degC);
       }
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double current_namp;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf ", date, time, &current_namp, &voltage, &temp_degC, &air_temp_degC, &rh_prcnt);
      //cout<<date<<","<<time<<","<<voltage<<","<<(float) current_namp<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
//        in >>  VDP_voltage >> VDP_capacitance >> VDP_conductance;
      //  in >>  VDP_voltage >> VDP_capacitance;
    	VDP_voltage_d[id][N_meas] = voltage*1E-9;  // to 1E-9 exei mpei gia na diorthwsei to lathos pou yparxei sta txt. tha prepei na aferethei an diorthwthei
    	VDP_current_d[id][N_meas] = current_namp;
      VDP_temp_meas.push_back(std::vector<double>());
      VDP_temp_meas[id].push_back(temp_degC);
      VDP_air_temp_meas.push_back(std::vector<double>());
      VDP_air_temp_meas[id].push_back(air_temp_degC);
      VDP_rh_prcnt_meas.push_back(std::vector<double>());
      VDP_rh_prcnt_meas[id].push_back(rh_prcnt);
      VDP_timestamp_meas.push_back(std::vector<string>());
      VDP_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
    	VDP_voltage[id][N_meas] = float(VDP_voltage_d[id][N_meas]);
      VDP_current[id][N_meas] = float(VDP_current_d[id][N_meas]);


    //        cout << N_meas << ")  "  <<  VDP_vv[N_meas] << "   " << VDP_cap[N_meas] << endl;
    //  cout << N_meas << ")  "  <<  VDP_voltage[id][N_meas] << "   " << VDP_current[id][N_meas]  << endl;
    	N_meas++;

      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();

   cout  << "----->>>>>>   Number of Measurements = " << N_meas << endl;
   VDP_Nmeas.push_back(N_meas);

   TGaxis::SetMaxDigits(3);

//   TCanvas *cc_spline = new TCanvas();
   VDP_cv[id] =  new TGraph(N_meas,VDP_current_d[id], VDP_voltage_d[id]);
   if((id % 12) == 0) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPPoly_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 1) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPPoly_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 2) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPPoly_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 3) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPPoly_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 4) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStop_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 5) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStop_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 6) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStop_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 7) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStop_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 8) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStrip_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 9) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStrip_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 10) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStrip_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 11) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStrip_r",HM_Id,Structure_Id[file_id].c_str()));
   VDP_cv[id]->GetXaxis()->SetTitle("VDP Current [A]");
   VDP_cv[id]->GetYaxis()->SetTitle("VDP Voltage [V]");
   VDP_cv[id]->SetDrawOption("AP");
   VDP_cv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   VDP_cv[id]->SetMarkerColor(colorid);
   //   VDP_cv[id]->Draw();
   // Set linear fit ranges
   // -------------------------
   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   VDP_cv_low[id] = VDP_current_d[id][3];
   VDP_cv_high[id] = VDP_current_d[id][38];
// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", VDP_cv_low[id], VDP_cv_high[id]);
   //   VDP_cv[id]->Fit("f_cv","0R+");
   VDP_cv[id]->Fit("f_cv","0R+");
   // Accumulation Region
   // ---------------------
   VDP_A_constant[id] = VDP_cv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   VDP_A_slope[id]  = VDP_cv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double Resistance = (TMath::Pi()/TMath::Log(2))*VDP_A_slope[id]; // VDP Resistance
   double Resistance_error = (TMath::Pi()/TMath::Log(2))*(VDP_cv[id]->GetFunction("f_cv")->GetParError(1));

   cout << "----->>>   Resistance Voltage = " << Resistance << " +/- " <<  Resistance_error<< " [Ohm/sq]" << endl;
   return Resistance;

}
//--------------------------------------------------------------------------------------------------//
// VDP Analysis Final
// --------------------------------------------------------------------------------------------------//
void VDP_PQC_Flute1_Analysis_Final(int nFiles)
{    //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1
   //   gROOT->Reset();
   gStyle->SetOptStat(0);

   //TString VDP_DataFileName[600];
   for(int i=0;i<nFiles;i++) {
    VDP_DataFileName[12*i+0] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPPoly_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPPoly_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPPoly_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPPoly_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+4] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStop_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+5] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStop_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+6] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStop_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+7] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStop_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+8] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStrip_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+9] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStrip_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+10] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStrip_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+11] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStrip_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {

     cout << endl << endl << "VDP_Poly Structure : Analyze file : " << VDP_DataFileName[12*i] << endl;
     VDPoly_s_Rsh_E[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+0],(12*i),i);
     cout << endl << endl << "VDP_Poly Structure : Analyze file : " << VDP_DataFileName[12*i+1] << endl;
     VDPoly_s_Rsh_W[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+1],(12*i+1),i);

     cout << endl << endl << "VDP_Poly Structure : Analyze file : " << VDP_DataFileName[12*i+2] << endl;
     VDPoly_r_Rsh_E[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+2],(12*i+2),i);
     cout << endl << endl << "VDP_Poly Structure : Analyze file : " << VDP_DataFileName[12*i+3] << endl;
     VDPoly_r_Rsh_W[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+3],(12*i+3),i);

     cout << endl << endl << "VDP_Stop Structure : Analyze file : " << VDP_DataFileName[12*i+4] << endl;
     VDPstop_s_Rsh_E[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+4],(12*i+4),i);
     cout << endl << endl << "VDP_Stop Structure : Analyze file : " << VDP_DataFileName[12*i+5] << endl;
     VDPstop_s_Rsh_W[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+5],(12*i+5),i);

     cout << endl << endl << "VDP_Stop Structure : Analyze file : " << VDP_DataFileName[12*i+6] << endl;
     VDPstop_r_Rsh_E[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+6],(12*i+6),i);
     cout << endl << endl << "VDP_Stop Structure : Analyze file : " << VDP_DataFileName[12*i+7] << endl;
     VDPstop_r_Rsh_W[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+7],(12*i+7),i);

     cout << endl << endl << "VDP_Strip Structure : Analyze file : " << VDP_DataFileName[12*i+8] << endl;
     VDPStrip_s_Rsh_E[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+8],(12*i+8),i);
     cout << "VDP_Strip Structure : Analyze file : " << VDP_DataFileName[12*i+9] << endl;
     VDPStrip_s_Rsh_W[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+9],(12*i+9),i);

     cout << endl << endl << "VDP_Strip Structure : Analyze file : " << VDP_DataFileName[12*i+10] << endl;
     VDPStrip_r_Rsh_E[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+10],(12*i+10),i);
     cout << endl << endl << "VDP_Strip Structure : Analyze file : " << VDP_DataFileName[12*i+11] << endl;
     VDPStrip_r_Rsh_W[i] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+11],(12*i+11),i);
   }
   for (int i=0;i<nFiles;i++) {
     for (int j=0;j<12;j++){
       VDP_Rsh_all[12*i+j] = VDP_Read_Info_from_File(VDP_DataFileName[12*i+j],(12*i+j),i);
     }
   }
   TCanvas *VDP_cc_final[12];
   TMultiGraph *VDP_mg[12];
   for (int i=0;i<12;i++) {
     VDP_cc_final[i] = new TCanvas(Form("cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     VDP_mg[i] = new TMultiGraph(Form("mg_%d",i),"VDP Resistance");
     if(i == 0) VDP_mg[i]->SetTitle("VDPPoly_s_E Resistance ");
     if(i == 1) VDP_mg[i]->SetTitle("VDPPoly_s_W Resistance ");
     if(i == 2) VDP_mg[i]->SetTitle("VDPPoly_r_E Resistance ");
     if(i == 3) VDP_mg[i]->SetTitle("VDPPoly_r_W Resistance ");
     if(i == 4) VDP_mg[i]->SetTitle("VDPStop_s_E Resistance ");
     if(i == 5) VDP_mg[i]->SetTitle("VDPStop_s_W Resistance ");
     if(i == 6) VDP_mg[i]->SetTitle("VDPStop_r_E Resistance ");
     if(i == 7) VDP_mg[i]->SetTitle("VDPStop_r_W Resistance ");
     if(i == 8) VDP_mg[i]->SetTitle("VDPStrip_s_E Resistance ");
     if(i == 9) VDP_mg[i]->SetTitle("VDPStrip_s_W Resistance ");
     if(i == 10) VDP_mg[i]->SetTitle("VDPStrip_r_E Resistance ");
     if(i == 11) VDP_mg[i]->SetTitle("VDPStrip_r_W Resistance ");

    // VDP_cv[12*i+]->Draw("ALP");
     for(int id=0;id<nFiles;id++) {
        VDP_mg[i]->Add(VDP_cv[12*id+i]);
     }
     VDP_mg[i]->Draw("LPsame");
     VDP_mg[i]->GetXaxis()->SetTitle("Current [A]");
     VDP_mg[i]->GetYaxis()->SetTitle("Voltage [V]");
     VDP_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();
     VDP_cc_final[i]->BuildLegend(0.15,0.5,0.5,0.85);
     gPad->Modified();
     gPad->Update();
   }

/*   old stuff
   auto cc_final_2 = new TCanvas("cc_final_2","",100,100,1000,1000);
   gPad->SetGrid();

   auto mg2 = new TMultiGraph("mg2","VDPPoly_r_W Resistance");
   mg2->SetTitle("VDPPoly_r_W Resistance ");

   VDP_cv[1]->Draw("ALP");
   for(int id=0;id<nFiles;id++) {

      mg2->Add(VDP_cv[12*id+1]);

   }
   mg2->Draw("LPsame");
   mg2->GetXaxis()->SetTitle("Current [A]");
   mg2->GetYaxis()->SetTitle("Voltage [V]");
   mg2->GetYaxis()->SetTitleOffset(1.35);
   gPad->Modified();
   gPad->Update();

   cc_final_2->BuildLegend(0.15,0.5,0.5,0.85);


   gPad->Modified();
   gPad->Update();


   for(int id=nFiles-2;id>-1;id--) {
       TF1 *func_A = new TF1("func_A","[0] + [1]*x",cv_low[id],cv1_low[id]);
       func_A->SetParameters(A_constant[id], A_slope[id]);
       func_A->Draw("same");

       TF1 *func_D = new TF1("func_D","[0] + [1]*x",cv_high[id],cv1_high[id]);
       func_D->SetParameters(D_constant[id], D_slope[id]);
       func_D->Draw("same");

       TF1 *func_I = new TF1("func_I","[0] + [1]*x",cv1_high[id],cv2_high[id]);
       func_I->SetParameters(I_constant[id], I_slope[id]);
       func_I->Draw("same");
   }

//   cc_final->Print(Form("CV_Flute1_MOS_SN_%d_250mV_%d_kHz.png",Structure_Id,freq),"png");
//   cc_final->Print(Form("CV_Flute1_MOS_SN_%d_250mV_%d_kHz.pdf",Structure_Id,freq),"pdf");
//   cc_final->Print(Form("CV_Flute1_MOS_SN_%d_250mV_%d_kHz.eps",Structure_Id,freq),"eps");
*/
   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;

   cout << "        VDP Id                 |  VDPoly_s_Rsh [kOhm/sq]   | VDPoly_r_Rsh [kOhm/sq]| VDPstop_s_Rsh [kOhm/sq]| VDPstop_r_Rsh [kOhm/sq] | VDPStrip_s_Rsh [Ohm/sq] | VDPStrip_r_Rsh [Ohm/sq]" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) << setw(20) <<  setprecision(3) << VDPoly_s_Rsh_E[id]/1000.  <<  setw(30)
	       << VDPoly_r_Rsh_E[id]/1000. <<  setw(25) <<  setprecision(4) << VDPstop_s_Rsh_E[id]/1000.  <<  setw(25)
	       << VDPstop_r_Rsh_E[id]/1000.  <<  setw(25)  << VDPStrip_s_Rsh_E[id] <<  setw(25)  << VDPStrip_r_Rsh_E[id]   << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << setw(20) <<  setprecision(3) << VDPoly_s_Rsh_W[id]/1000.  <<  setw(30)
	       << VDPoly_r_Rsh_W[id]/1000. <<  setw(25) <<  setprecision(4) << VDPstop_s_Rsh_W[id]/1000.  <<  setw(25)
	       << VDPstop_r_Rsh_W[id]/1000.  <<  setw(25)  << VDPStrip_s_Rsh_W[id] <<  setw(25)  << VDPStrip_r_Rsh_W[id]   << endl;
   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------------------------------------------//
// VDP xml production
// --------------------------------------------------------------------------------------------------//
void VDP_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  cout << "VDP Xml production is on" << endl;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<12; j++){
          //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = VDP_DataFileName[12*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = VDP_DataFileName[12*i+j].Last(delimeter2);
         xml_filename = VDP_DataFileName[12*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[4*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute1/VPX%d/",HM_Id) +xml_filename);

        // First create engine
        TXMLEngine VDP_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = VDP_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = VDP_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = VDP_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = VDP_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = VDP_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = VDP_xml.NewChild(HEADER, 0, "RUN");
              VDP_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              VDP_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              VDP_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              VDP_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              VDP_xml.NewChild(RUN, 0, "INITIATED_BY_USER", VDP_Operators[12*i+j].c_str());
              VDP_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", VDP_Begin_Timestamp[12*i+j].c_str());
              VDP_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = VDP_xml.NewChild(ROOT, 0, "DATA_SET");
            VDP_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            VDP_xml.NewChild(DATA_SET, 0, "VERSION", "IV_measurement");
            XMLNodePointer_t PART = VDP_xml.NewChild(DATA_SET, 0, "PART");
              VDP_xml.NewChild(PART, 0, "NAME_LABEL", VDP_name_labels[12*i+j].c_str());
              VDP_xml.NewChild(PART, 0, "KIND_OF_PART", VDP_kind_of_parts[12*i+j].c_str());
            XMLNodePointer_t DATA = VDP_xml.NewChild(DATA_SET, 0, "DATA");
              VDP_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC1");
              VDP_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", VDP_struct_id[12*i+j].c_str());
              VDP_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", VDP_set_id[12*i+j].c_str());
              VDP_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", VDP_config_id[12*i+j].c_str());
              VDP_xml.NewChild(DATA, 0, "EQUIPMENT", VDP_equipment[12*i+j].c_str());
              VDP_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(VDP_waiting_time[12*i+j], 3).c_str());
              VDP_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(VDP_temp[12*i+j], 3).c_str());
              VDP_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(VDP_av_temp[12*i+j], 3).c_str());
              VDP_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", VDP_struct_id[12*i+j].c_str());
              VDP_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = VDP_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = VDP_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = VDP_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  VDP_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  VDP_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = VDP_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                VDP_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                VDP_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "IV_measurement");
                XMLNodePointer_t PART_CDS1 = VDP_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  VDP_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", VDP_name_labels[12*i+j].c_str());
                  VDP_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", VDP_kind_of_parts[12*i+j].c_str());
                for (int k=0; k<VDP_Nmeas[12*i+j]-1; k++){
                    XMLNodePointer_t DATA_CDS1 = VDP_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                      VDP_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string_scientific_format(VDP_voltage[12*i+j][k], 3).c_str());
                      VDP_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(VDP_current[12*i+j][k]*1E+9, 3).c_str());
                      VDP_xml.NewChild(DATA_CDS1, 0, "TIME", VDP_timestamp_meas[12*i+j][k].c_str());
                      VDP_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(VDP_temp_meas[12*i+j][k], 3).c_str());
                      VDP_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(VDP_air_temp_meas[12*i+j][k], 3).c_str());
                      VDP_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(VDP_rh_prcnt_meas[12*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = VDP_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = VDP_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = VDP_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      VDP_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      VDP_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = VDP_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    VDP_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    VDP_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement");
                    XMLNodePointer_t PART_CDS2 = VDP_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      VDP_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", VDP_name_labels[12*i+j].c_str());
                      VDP_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", VDP_kind_of_parts[12*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = VDP_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      VDP_xml.NewChild(DATA_CDS2, 0, "RSH_OHMSQR", Conv_float_to_string(VDP_Rsh_all[12*i+j], 3).c_str());
                      VDP_xml.NewChild(DATA_CDS2, 0, "R_OHM", Conv_float_to_string(VDP_Rsh_all[12*i+j]*(TMath::Log(2)/TMath::Pi()), 3).c_str());

        XMLDocPointer_t VDP_xmldoc = VDP_xml.NewDoc();
        VDP_xml.DocSetRootElement(VDP_xmldoc, ROOT);
        // Save document to file
        VDP_xml.SaveDoc(VDP_xmldoc, Form("./XML_Info_Flute1/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        VDP_xml.FreeDoc(VDP_xmldoc);
      }
    }
  }
}
//====================================================================================================//
//--------------------------------------------------------------------------------------------------//
// FET Read data
// --------------------------------------------------------------------------------------------------//
void FET_Read_Info_from_File(TString Analysis_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(Analysis_DataFileName);

   int nlines = 0;
   int N_meas = 0;

   string line;

   float FET_Voltage, FET_Current;
   cout << in.good()<< '\n';
   unsigned int lines_read = 0U;
   //while ((lines_read < 10) && (std::getline(text_file_stream, text_from_file)))
   //{
  //++lines_read;
  //}
   while (in.good()) {
     getline (in,line);
      if(nlines < 18) {
         if(nlines == 1) {
         //cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         FET_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         FET_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         FET_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         FET_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         FET_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         FET_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         FET_config_id.push_back(Kind_of_HM_config_id);
         }
         if(nlines == 9) {
         cout << line << '\n';
         char Procedure_type[20];
         sscanf(line.c_str(), "%*s %s", Procedure_type);
         cout<<Procedure_type<<endl;
         }
         if(nlines == 10) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string  Equipment =  Conv_Char_To_String(str1) +" "+Conv_Char_To_String(str2);
         cout<<Equipment<<endl;
         FET_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         FET_waiting_time.push_back(Waiting_time_s);
         }
         if(nlines == 12) {
         cout << line << '\n';
         double Bias_V;
         sscanf(line.c_str(), "%*s %lf", &Bias_V);
         cout<<Bias_V<<endl;
         FET_bias.push_back(Bias_V);
         }
         if(nlines == 13) {
         cout << line << '\n';
         double Temp_set_degC;
         sscanf(line.c_str(), "%*s %lf", &Temp_set_degC);
         cout<<Temp_set_degC<<endl;
         FET_temp.push_back(Temp_set_degC);
        }
        if(nlines == 14) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        FET_av_temp.push_back(Av_temp_degC);
       }
        if(nlines == 15) {
        cout << line << '\n';
        double ComplianceA;
        sscanf(line.c_str(), "%*s %lf", &ComplianceA);
        cout<<ComplianceA<<endl;
        FET_compliance.push_back(ComplianceA);
      }
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double current_namp;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf", date, time, &voltage, &current_namp, &temp_degC, &air_temp_degC, &rh_prcnt);
      //sscanf(line.c_str(), "%s", date, time);
      //in >> date >> time >> voltage >> current_namp >> temp_degC >> air_temp_degC >> rh_prcnt;
      //cout<<date<<","<<time<<","<<voltage<<","<<(float) current_namp<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;

    	FET_vv_d[id][N_meas] = voltage;
    	FET_Curr_d[id][N_meas] = current_namp;
      FET_temp_meas.push_back(std::vector<double>());
      FET_temp_meas[id].push_back(temp_degC);
      FET_air_temp_meas.push_back(std::vector<double>());
      FET_air_temp_meas[id].push_back(air_temp_degC);
      FET_rh_prcnt_meas.push_back(std::vector<double>());
      FET_rh_prcnt_meas[id].push_back(rh_prcnt);
      FET_timestamp_meas.push_back(std::vector<string>());
      FET_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
    	FET_vv[id][N_meas] = (float) FET_vv_d[id][N_meas];
      FET_Curr[id][N_meas] = (float) FET_Curr_d[id][N_meas];
      //cout<< FET_vv[N_meas] << "," << FET_Curr[N_meas] << endl;
      N_meas++;
      }
      if (!in.good()) break;
      nlines++;
   }
   //cout << N_meas << ")  "  <<  FET_vv[N_meas] << "   " << FET_Curr[N_meas]  << endl;
   printf(" Found %d lines\n",nlines);
   in.close();

   cout  << "----->>>>>>   Number of Measurements = " << N_meas << endl;
   FET_Nmeas.push_back(N_meas);

   int  N_meas_reduced = N_meas-1;

   if(N_meas_reduced > 300 ) {
      cout << "Serious Problem : Measurements > 300 that TSpline3 can handle...  " << endl;
      cout << " Solution: Reduce number of measurements for TSpline3 rootine.... " << endl;
      cout << " Now Abort..." << endl;
      abort();
   }
   double FET_vv_d_reduced[300];
   double FET_Current_d_reduced[300];

   for (Int_t i=0; i<N_meas_reduced; i++) {
      FET_vv_d_reduced[i] = FET_vv_d[id][i];
      FET_Current_d_reduced[i] = FET_Curr_d[id][i];
   }

   FET_CuSpl_iv[id] = new TSpline3("Cubic Spline", FET_vv_d_reduced, FET_Current_d_reduced, N_meas_reduced, "cb2e2", 0, 0);  // Seems that works well for less than 150 points
   //TCanvas *cc_spline = new TCanvas();
   FET_iv[id] =  new TGraph(N_meas_reduced,FET_vv_d_reduced,FET_Current_d_reduced);

   if((id % 2) == 0) FET_iv[id]->SetTitle(Form("FET%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) FET_iv[id]->SetTitle(Form("FET%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));

   FET_iv[id]->GetXaxis()->SetTitle("Gate FET_Voltage [V]");
   FET_iv[id]->GetYaxis()->SetTitle("Drain FET_Current [#muA]");
   FET_iv[id]->SetDrawOption("AP");
   FET_iv[id]->SetMarkerStyle(20+id);
   int colorid = 1+id;
   if(colorid == 5) colorid = 46;
   FET_iv[id]->SetMarkerColor(colorid);
   //FET_iv[id]->Draw();
   FET_CuSpl_iv[id]->SetLineColor(kBlue); // Plot the CSpline curve as well
   //FET_CuSpl_iv[id]->Draw("same");

    // Get the First Derivative of IV curve
    // --------------------------------------
    for(int i = 0;i<N_meas_reduced; i++) {
       FET_IV_CubicSpline_first_Derivative[i] = FET_CuSpl_iv[id]->Derivative(FET_vv_d_reduced[i]);
    }

    // Make Cubic Spline fit of the First Derivative
    FET_CuSpl_firstDerivative_iv[id] = new TSpline3("Cubic Spline", FET_vv_d_reduced, FET_IV_CubicSpline_first_Derivative , N_meas_reduced, "cb2e2", 0, 0);

    // get the maximum of the first derivative
    int max_index_simple = 0;
    float max_FET_IV_CubicSpline_first_Derivative = FET_IV_CubicSpline_first_Derivative[0];
    for (int i = 1;i<N_meas_reduced; i++) {
       if(FET_IV_CubicSpline_first_Derivative[i] >  max_FET_IV_CubicSpline_first_Derivative && FET_vv_d_reduced[i]>0 && FET_vv_d_reduced[i] < 5) {
         max_FET_IV_CubicSpline_first_Derivative = FET_IV_CubicSpline_first_Derivative[i];
	       max_index_simple = i;
       }
    }
    double FET_Gate_FET_Voltage_max_first_Derivative = double(FET_vv_d[id][max_index_simple]);
    cout  << "---------------- >>>>> Max First derivative of FET IV curve at FET_Voltage index = " << max_index_simple
	  << " with FET Gate FET_Voltage   = " <<  FET_Gate_FET_Voltage_max_first_Derivative << " V " << endl;

    //TGraph *cv_deriv =  new TGraph(N_meas_reduced,FET_vv_d_reduced,FET_IV_CubicSpline_first_Derivative);
    //cv_deriv->SetTitle(Form("First Derivative of FET IV"));
   // cv_deriv->GetXaxis()->SetTitle("Gate FET_Voltage [V]");
    //cv_deriv->Draw("same");
    FET_CuSpl_firstDerivative_iv[id]->SetLineColor(kRed);
    //FET_CuSpl_firstDerivative_iv[id]->Draw("same");

    // Get the Second Derivative of IV curve
    // --------------------------------------
    for(int i = 0;i<N_meas_reduced; i++) {
       FET_IV_CubicSpline_second_Derivative[i] = FET_CuSpl_firstDerivative_iv[id]->Derivative(FET_vv_d_reduced[i]);
    }

    //TGraph *cv_secderiv =  new TGraph(N_meas_reduced,FET_vv_d_reduced,FET_IV_CubicSpline_second_Derivative);
   // cv_secderiv->SetTitle(Form("Second Derivative of FET CV"));
    //cv_secderiv->GetXaxis()->SetTitle("Gate FET_Voltage [V]");
    //cv_secderiv->SetLineColor(kBlue);
    //cv_secderiv->Draw("same");

    int Nub_steps = 1000;
    float FET_Voltage_extented[2000];
    float secdev_spline[2000];
    float secdev_spline_f[2000];
    float vv_initial = FET_vv_d_reduced[0];
    float vv_final = FET_vv_d_reduced[N_meas_reduced-1];
    for(int i = 0;i<Nub_steps; i++) {
      FET_Voltage_extented[i] =  vv_initial + i*(vv_final - vv_initial)/(Nub_steps-1);
      FET_IV_CubicSpline_second_Derivative_ext[i] = FET_CuSpl_firstDerivative_iv[id]->Derivative(FET_Voltage_extented[i]);
    }

    // get the maximum of the second derivative
    int max_index_secdev = 0;
    float max_FET_IV_CubicSpline_second_Derivative = FET_IV_CubicSpline_second_Derivative_ext[0];
    for (int i = 1;i<Nub_steps; i++) {
       if(FET_IV_CubicSpline_second_Derivative_ext[i] >  max_FET_IV_CubicSpline_second_Derivative && FET_Voltage_extented[i] < FET_Gate_FET_Voltage_max_first_Derivative) {
         max_FET_IV_CubicSpline_second_Derivative = FET_IV_CubicSpline_second_Derivative_ext[i];
	       max_index_secdev = i;
       }
    }
    FET_Gate_Voltage_Threshold[id] = double(FET_Voltage_extented[max_index_secdev]);
    cout  << "---------------- >>>>> Max Second derivative of FET IV curve at FET_Voltage index = " <<  max_index_secdev
	  << " and FET Gate FET_Voltage Threshold = " <<  FET_Gate_Voltage_Threshold[id] << " V " << endl;

}
//--------------------------------------------------------------------------------------------------//
// FET Analysis Final
// --------------------------------------------------------------------------------------------------//
void FET_PQC_Analysis_Final(int nFiles)
{

 //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1

 //   gROOT->Reset();
   gStyle->SetOptStat(0);


   for(int i=0;i<nFiles;i++) {
    Analysis_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Left_PQC1_FETPSs_trans.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    Analysis_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Left_PQC1_FETPSs_trans.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {
     cout << endl << endl << "FET Structure : Analyze file : " << Analysis_DataFileName[2*i] << endl;
     FET_Read_Info_from_File(Analysis_DataFileName[2*i],2*i,i);
     cout << endl << endl << "FET Structure : Analyze file : " << Analysis_DataFileName[2*i+1] << endl;
     FET_Read_Info_from_File(Analysis_DataFileName[2*i+1],2*i+1,i);
  }

   TCanvas *FET_cc_final[2];
   TMultiGraph *FET_mg[2];
   for (int i=0;i<2;i++) {
     FET_cc_final[i] = new TCanvas(Form("FET_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     FET_mg[i] = new TMultiGraph(Form("mg_%d",i),"FET IV");
     if(i == 0) FET_mg[i]->SetTitle("FET_E Gate_FET_Voltage_Threshold ");
     if(i == 1) FET_mg[i]->SetTitle("FET_W Gate_FET_Voltage_Threshold ");

     for(int id=0;id<nFiles;id++) {
        FET_mg[i]->Add(FET_iv[2*id+i]);
     }
     FET_mg[i]->Draw("LPsame");
     FET_mg[i]->GetXaxis()->SetTitle("FET_Current [A]");
     FET_mg[i]->GetYaxis()->SetTitle("FET_Voltage [V]");
     FET_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     FET_cc_final[i]->BuildLegend(0.15,0.5,0.35,0.85);


     gPad->Modified();
     gPad->Update();
   }

   cout <<endl;
   cout << "===========================+==================================================================" << endl;
   cout << "---------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  ---------------------------------" << endl;
   cout << "==============================================================================================" << endl;

   cout << "        FET Id                 |  VFB_Deb [V] |" <<endl;
   cout << "----------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
      cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) << setw(15)  <<  setprecision(3) << FET_Gate_Voltage_Threshold[2*id]   << endl;
      cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << setw(15)  <<  setprecision(3) << FET_Gate_Voltage_Threshold[2*id+1]   << endl;

   }
   cout << "----------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------------------------------------------//
// FET xml production
// --------------------------------------------------------------------------------------------------//
void FET_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    cout << "FET Xml production is on" << endl;
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<=1; j++){
  //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = Analysis_DataFileName[2*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = Analysis_DataFileName[2*i+j].Last(delimeter2);
         xml_filename = Analysis_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
 //        cout << MOS_DataFileName[2*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute1/VPX%d/",HM_Id) +xml_filename);


         // First create engine
        TXMLEngine FET_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = FET_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = FET_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = FET_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = FET_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = FET_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = FET_xml.NewChild(HEADER, 0, "RUN");
              FET_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              FET_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              FET_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              FET_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              FET_xml.NewChild(RUN, 0, "INITIATED_BY_USER", FET_Operators[2*i+j].c_str());
              FET_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", FET_Begin_Timestamp[2*i+j].c_str());
              FET_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = FET_xml.NewChild(ROOT, 0, "DATA_SET");
            FET_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            FET_xml.NewChild(DATA_SET, 0, "VERSION", "v2");
            XMLNodePointer_t PART = FET_xml.NewChild(DATA_SET, 0, "PART");
              FET_xml.NewChild(PART, 0, "NAME_LABEL", FET_name_labels[2*i+j].c_str());
              FET_xml.NewChild(PART, 0, "KIND_OF_PART", FET_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = FET_xml.NewChild(DATA_SET, 0, "DATA");
              FET_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC1");
              FET_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", FET_struct_id[2*i+j].c_str());
              FET_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", FET_set_id[2*i+j].c_str());
              FET_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", FET_config_id[2*i+j].c_str());
              FET_xml.NewChild(DATA, 0, "EQUIPMENT", FET_equipment[2*i+j].c_str());
              FET_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(FET_waiting_time[2*i+j],3).c_str());
              FET_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(FET_temp[2*i+j],3).c_str());
              FET_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(FET_av_temp[2*i+j],3).c_str());
              FET_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "FET1");
              FET_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = FET_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = FET_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = FET_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  FET_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  FET_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = FET_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                FET_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                FET_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "V2");
                XMLNodePointer_t PART_CDS1 = FET_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  FET_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", FET_name_labels[2*i+j].c_str());
                  FET_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", FET_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<FET_Nmeas[2*i+j]-1; k++){
                XMLNodePointer_t DATA_CDS1 = FET_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  FET_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(FET_vv[2*i+j][k],3).c_str());
                  //FET_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", to_string(FET_Curr[k]).c_str());
                  FET_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(FET_Curr[2*i+j][k],3).c_str());
                  FET_xml.NewChild(DATA_CDS1, 0, "TIME", FET_timestamp_meas[2*i+j][k].c_str());
                  FET_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(FET_temp_meas[2*i+j][k], 3).c_str());
                  FET_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(FET_air_temp_meas[2*i+j][k], 3).c_str());
                  FET_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(FET_rh_prcnt_meas[2*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = FET_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = FET_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = FET_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      FET_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_TC_PAR");
                      FET_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon TC Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = FET_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    FET_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    FET_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "TC_measurement-004");
                    XMLNodePointer_t PART_CDS2 = FET_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      FET_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", FET_name_labels[2*i+j].c_str());
                      FET_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", FET_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = FET_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      FET_xml.NewChild(DATA_CDS2, 0, "VTH_V", Conv_float_to_string(FET_Gate_Voltage_Threshold[2*i+j], 5).c_str());
        XMLDocPointer_t FET_xmldoc = FET_xml.NewDoc();
        FET_xml.DocSetRootElement(FET_xmldoc, ROOT);
        // Save document to file
        FET_xml.SaveDoc(FET_xmldoc, Form("./XML_Info_Flute1/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        FET_xml.FreeDoc(FET_xmldoc);
      }
    }
  }
}
// ==================================================================================================//
//--------------------------------------------------------------------------------------------------//
// CAP Read data
// --------------------------------------------------------------------------------------------------//
void CAP_Read_Info_from_File(TString CAP_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(CAP_DataFileName);

   int nlines = 0;
   int N_meas = 0;

   string line;

   float CAP_voltage, CAP_capacitance, CAP_conductance;
   while (in.good()) {
     getline (in,line);
      if(nlines < 17) {
         if(nlines == 1) {
         cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         TString Operator = TString(First_Name)+TString(Last_Name);
         cout<<"The Operator is:" << Operator<<endl;
         CAP_Operators.push_back(Conv_Char_To_String(First_Name)
                                +Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         //cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         //cout<<Date<<Time<<endl;
         CAP_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         //cout << Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time) << endl;
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         CAP_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         CAP_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         CAP_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         CAP_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         CAP_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         CAP_config_id.push_back(Kind_of_HM_config_id);
         }
         if(nlines == 9) {
         cout << line << '\n';
         char Procedure_type[20];
         sscanf(line.c_str(), "%*s %s", Procedure_type);
         cout<<Procedure_type<<endl;
         }
         if(nlines == 10) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         CAP_waiting_time.push_back(Waiting_time_s);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double ac_freq;
         sscanf(line.c_str(), "%*s %lf", &ac_freq);
         cout<<ac_freq<<endl;
         CAP_ac_freq.push_back(ac_freq);
         }
         if(nlines == 12) {
         cout << line << '\n';
         double ac_ampl;
         sscanf(line.c_str(), "%*s %lf", &ac_ampl);
         cout<<ac_ampl<<endl;
         CAP_ac_ampl.push_back(ac_ampl);
        }
        if(nlines == 13) {
        cout << line << '\n';
        double temp_degC;
        sscanf(line.c_str(), "%*s %lf", &temp_degC);
        cout<<temp_degC<<endl;
        CAP_temp.push_back(temp_degC);
       }
        if(nlines == 14) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        CAP_av_temp.push_back(Av_temp_degC);
       }
    //    nlines++;
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double capctnc_pfrd, rsstnc_mhom;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf %lf", date, time, &voltage, &capctnc_pfrd, &rsstnc_mhom, &temp_degC, &air_temp_degC, &rh_prcnt);

    	CAP_vv_d[id][N_meas] = voltage;
    	CAP_Capacitance_d[id][N_meas] = capctnc_pfrd;

      CAP_temp_meas.push_back(std::vector<double>());
      CAP_temp_meas[id].push_back(temp_degC);
      CAP_air_temp_meas.push_back(std::vector<double>());
      CAP_air_temp_meas[id].push_back(air_temp_degC);
      CAP_rh_prcnt_meas.push_back(std::vector<double>());
      CAP_rh_prcnt_meas[id].push_back(rh_prcnt);
      CAP_timestamp_meas.push_back(std::vector<string>());
      CAP_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      CAP_rsstnc_mhom.push_back(std::vector<double>());
      CAP_rsstnc_mhom[id].push_back(rsstnc_mhom);
    //        cout << N_meas << ")  "  <<  CAP_vv[N_meas] << "   " << CAP_Capacitance[N_meas]  << endl;
      N_meas++;
    //  nlines++;
      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   CAP_Nmeas.push_back(N_meas);
   cout  << "----->>>>>>   Number of Measurements = " << N_meas << endl;
   //TCanvas *cc_spline = new TCanvas();
   CAP_cv[id] =  new TGraph(N_meas,CAP_vv_d[id],CAP_Capacitance_d[id]);

   if((id % 2) == 0) CAP_cv[id]->SetTitle(Form("CAP%d_0%s CV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) CAP_cv[id]->SetTitle(Form("CAP%d_0%s CV",HM_Id,Structure_Id[file_id].c_str()));

   CAP_cv[id]->GetXaxis()->SetTitle("CAP_voltage [V]");
   CAP_cv[id]->GetYaxis()->SetTitle("CAP_capacitance [pF]");
   CAP_cv[id]->SetDrawOption("AP");
   CAP_cv[id]->SetMarkerStyle(20+id);
   int colorid = 1+id;
   if(colorid == 5) colorid = 46;
   CAP_cv[id]->SetMarkerColor(colorid);
   //CAP_cv[id]->Draw();

   TF1 *f_cv = new TF1("f_cv","pol0", -5, 5);
   CAP_cv[id]->Fit("f_cv","0R+");

   double A_cap = 0.0130*0.0130; //cm^2
   CAP_Value[id] = CAP_cv[id]->GetFunction("f_cv")->GetParameter(0); //
   double dox = ((epsilon_Si*A_cap)/CAP_Value[id])*1E+7;
   CAP_dox[id] =  dox ;
}
//--------------------------------------------------------------------------------------------------//
// CAP Analysis Final
// --------------------------------------------------------------------------------------------------//
void CAP_PQC_Analysis_Final(int nFiles)
{
   //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1
   //   gROOT->Reset();
   gStyle->SetOptStat(0);

  // TString CAP_DataFileName[100];
   for(int i=0;i<nFiles;i++) {
    CAP_DataFileName[4*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/CV/VPX%d_0%s_2-S_HM_E_Left_PQC1_CapW_CV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    CAP_DataFileName[4*i+1] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/CV/VPX%d_0%s_2-S_HM_E_Left_PQC1_CapE_CV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    CAP_DataFileName[4*i+2] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/CV/VPX%d_0%s_2-S_HM_W_Left_PQC1_CapW_CV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    CAP_DataFileName[4*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/CV/VPX%d_0%s_2-S_HM_W_Left_PQC1_CapE_CV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }
   for (int i=0;i<nFiles;i++) {
     cout << endl << endl << "CAP structure : Analyze file : " << CAP_DataFileName[4*i] << endl;
     CAP_Read_Info_from_File(CAP_DataFileName[4*i],4*i,i);
     cout << endl << endl << "CAP structure : Analyze file : " << CAP_DataFileName[4*i+1] << endl;
     CAP_Read_Info_from_File(CAP_DataFileName[4*i+1],4*i+1,i);
     cout << endl << endl << "CAP structure : Analyze file : " << CAP_DataFileName[4*i+2] << endl;
     CAP_Read_Info_from_File(CAP_DataFileName[4*i+2],4*i+2,i);
     cout << endl << endl << "CAP structure : Analyze file : " << CAP_DataFileName[4*i+3] << endl;
     CAP_Read_Info_from_File(CAP_DataFileName[4*i+3],4*i+3,i);
  }
   TCanvas *CAP_cc_final[4];
   TMultiGraph *CAP_mg[4];
   for (int i=0;i<4;i++) {
     CAP_cc_final[i] = new TCanvas(Form("CAP_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();
     CAP_mg[i] = new TMultiGraph(Form("CAP_mg_%d",i),"CAP IV");
     if(i == 0) CAP_mg[i]->SetTitle("East HM CAP_W Gate_CAP_voltage_Threshold ");
     if(i == 1) CAP_mg[i]->SetTitle("East HM CAP_E Gate_CAP_voltage_Threshold ");
     if(i == 2) CAP_mg[i]->SetTitle("West HM CAP_W Gate_CAP_voltage_Threshold ");
     if(i == 3) CAP_mg[i]->SetTitle("West HM CAP_E Gate_CAP_voltage_Threshold ");
     for(int id=0;id<nFiles;id++) {
        CAP_mg[i]->Add(CAP_cv[4*id+i]);
     }
     CAP_mg[i]->Draw("LPsame");
     CAP_mg[i]->GetXaxis()->SetTitle("CAP_capacitance [pF]");
     CAP_mg[i]->GetYaxis()->SetTitle("CAP_voltage [V]");
     CAP_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();
     CAP_cc_final[i]->BuildLegend(0.15,0.5,0.35,0.85);
     gPad->Modified();
     gPad->Update();
   }
   cout <<endl;
   cout << "===========================+==================================================================" << endl;
   cout << "---------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  ---------------------------------" << endl;
   cout << "==============================================================================================" << endl;
   cout << "        CAP Id                 |  CAP_capacitance W [pF] |   CAP_capacitance E [pF] |" <<endl;
   cout << "----------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
      cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
            << setw(15)  <<  setprecision(3) << CAP_Value[4*id]  << setw(15)  <<  setprecision(3) << CAP_Value[4*id+1] << endl;
      cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
            << setw(15)  <<  setprecision(3) << CAP_Value[4*id+2] << setw(15)  <<  setprecision(3) << CAP_Value[4*id+3]  << endl;

   }
   cout << "----------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------------------------------------------//
// CAP xml production
// --------------------------------------------------------------------------------------------------//
void CAP_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    cout << "CAP Xml production is on" << endl;
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<4; j++){
        //       TString xml_filename[nFiles];
        char delimeter1 = '/';
        int pos1 = CAP_DataFileName[4*i+j].Last(delimeter1);
        char delimeter2 = '.';
        int pos2 = CAP_DataFileName[4*i+j].Last(delimeter2);
        xml_filename = CAP_DataFileName[4*i+j](pos1+1, ((pos2-pos1)-1))+".xml";
        //        cout << CAP_DataFileName[4*i+j](0,pos1) << endl;

        cout << xml_filename << endl;
        ofstream xml_file;

        xml_file.open(Form("./XML_Info_Flute1/VPX%d/",HM_Id)+xml_filename);


        //cout << MOS_Operators[0] << endl;
        // First create engine
        TXMLEngine CAP_xml;
        // Create main node of document tree
        XMLNodePointer_t ROOT = CAP_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = CAP_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = CAP_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = CAP_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = CAP_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = CAP_xml.NewChild(HEADER, 0, "RUN");
              CAP_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              CAP_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              CAP_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              CAP_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              CAP_xml.NewChild(RUN, 0, "INITIATED_BY_USER", CAP_Operators[4*i+j].c_str());
              CAP_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", CAP_Begin_Timestamp[4*i+j].c_str());
              CAP_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = CAP_xml.NewChild(ROOT, 0, "DATA_SET");
            CAP_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            CAP_xml.NewChild(DATA_SET, 0, "VERSION", "CV-measurement");
            XMLNodePointer_t PART = CAP_xml.NewChild(DATA_SET, 0, "PART");
              CAP_xml.NewChild(PART, 0, "NAME_LABEL", CAP_name_labels[4*i+j].c_str());
              CAP_xml.NewChild(PART, 0, "KIND_OF_PART", CAP_kind_of_parts[4*i+j].c_str());
            XMLNodePointer_t DATA = CAP_xml.NewChild(DATA_SET, 0, "DATA");
              CAP_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC1");
              CAP_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", CAP_struct_id[4*i+j].c_str());
              CAP_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", CAP_set_id[4*i+j].c_str());
              CAP_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", CAP_config_id[4*i+j].c_str());
              CAP_xml.NewChild(DATA, 0, "EQUIPMENT", "HP 4192A");
              CAP_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(CAP_waiting_time[4*i+j],4).c_str());
              CAP_xml.NewChild(DATA, 0, "AC_FREQ_HZ", Conv_float_to_string(CAP_ac_freq[4*i+j],5).c_str());
              CAP_xml.NewChild(DATA, 0, "AC_AMPL_V", Conv_float_to_string(CAP_ac_ampl[4*i+j],3).c_str());
              CAP_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(CAP_temp[4*i+j],3).c_str());
              CAP_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(CAP_av_temp[4*i+j],3).c_str());
              CAP_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "Cap1");
//              CAP_xml.NewChild(DATA, 0, "FILE_NAME", xml_filename[4*i+j]);
              CAP_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = CAP_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = CAP_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = CAP_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  CAP_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_CV");
                  CAP_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon CV Test");
              XMLNodePointer_t DATA_SET_CDS1 = CAP_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                CAP_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                CAP_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "CV-measurement");
                XMLNodePointer_t PART_CDS1 = CAP_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  CAP_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", CAP_name_labels[4*i+j].c_str());
                  CAP_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", CAP_kind_of_parts[4*i+j].c_str());
                for (int k=0; k<CAP_Nmeas[4*i+j]-1; k++){
                XMLNodePointer_t DATA_CDS1 = CAP_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  CAP_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(CAP_vv_d[4*i+j][k],3).c_str());
                  CAP_xml.NewChild(DATA_CDS1, 0, "CAPCTNC_PFRD", Conv_float_to_string(CAP_Capacitance_d[4*i+j][k],3).c_str());
                  CAP_xml.NewChild(DATA_CDS1, 0, "RESSTNC_MOHM", Conv_float_to_string(CAP_rsstnc_mhom[4*i+j][k],3).c_str());
                  CAP_xml.NewChild(DATA_CDS1, 0, "TIME", CAP_timestamp_meas[4*i+j][k].c_str());
                  CAP_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(CAP_temp_meas[4*i+j][k],3).c_str());
                  CAP_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(CAP_air_temp_meas[4*i+j][k],3).c_str());
                  CAP_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(CAP_rh_prcnt_meas[4*i+j][k],3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = CAP_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = CAP_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = CAP_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      CAP_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_CV_PAR");
                      CAP_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon CV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = CAP_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    CAP_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    CAP_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "CV-measurement");
                    XMLNodePointer_t PART_CDS2 = CAP_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      CAP_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", CAP_name_labels[4*i+j].c_str());
                      CAP_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", CAP_kind_of_parts[4*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = CAP_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      CAP_xml.NewChild(DATA_CDS2, 0, "CAC_PFRD", Conv_float_to_string(CAP_Value[4*i+j],5).c_str());
                      CAP_xml.NewChild(DATA_CDS2, 0, "DOX_NM", Conv_float_to_string(CAP_dox[4*i+j],5).c_str());

        XMLDocPointer_t CAP_xmldoc = CAP_xml.NewDoc();
        CAP_xml.DocSetRootElement(CAP_xmldoc, ROOT);

        // Save document to file
        CAP_xml.SaveDoc(CAP_xmldoc, Form("./XML_Info_Flute1/VPX%d/",HM_Id)+xml_filename);

        // Release memory before exit
        CAP_xml.FreeDoc(CAP_xmldoc);
      }
    }
  }
}
// ========================================================================================================================================

void Flute1_2_S_Quick_Characterization_with_xml(int nFiles)
{
   CAP_PQC_Analysis_Final(nFiles);
   CAP_PQC_xml_production(xml_on, nFiles);

   VDP_PQC_Flute1_Analysis_Final(nFiles);
   VDP_PQC_xml_production(xml_on, nFiles);

   MOS_PQC_Analysis_Final(nFiles);
   MOS_PQC_xml_production(xml_on, nFiles);

   FET_PQC_Analysis_Final(nFiles);
   FET_PQC_xml_production(xml_on, nFiles);

   CVS_Output_File.open(Form("CSV_Info/VPX%d/Flute1_VPX%d_0xx_2-S_HM_E_W_Left_Quick_Info.csv",HM_Id,HM_Id));


   cout << endl << endl;
   cout << "==============================================================================================================================================" << endl;
   cout << "------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  --------------------------------------------------" << endl;
   cout << "==============================================================================================================================================" << endl;
   cout << "     Flute1  Id            |  CAP_E  | CAP_W  |   VDPPoly   |   VDPStop   |  VDPStrip   | FET_Vth |  C_acc  |  t_ox  |  VFB_Deb  | Nox_Deb   | " << endl;
   cout << "                           |  [pF]   |  [pF]  |  [kOhm/sq]  |  [kOhm/sq]  |  [Ohm/sq]   |   [V]   |   [pF]  |  [um]  |    [V]    | [cm^{-2}] | " << endl;
   cout << "                           |         |        |  stad | rot | stad |  rot | stad |  rot |         |         |        |           |           | " << endl;
   cout << "     Spec Limit            | >2.028  | >2.028 |  2.2- |2.2- |18-20 |18-20 |32-35 |32-35 |         |  ~82    |  0.7   |    <5     |           | " << endl;
   cout << "     Spec Limit            |         |        |   2.4 | 2.4 |      |      |      |      |         |         |        |           |           | " << endl;
   cout << "---------------------------------------------------------------------------------------------------------------------------------------------- " << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(8) <<  setprecision(3) << CAP_Value[4*id]  << setw(10)  <<  setprecision(3) << CAP_Value[4*id+1]
	       << setw(8) <<  setprecision(3) << VDPoly_s_Rsh_E[id]/1000.  <<  setw(7)
	       << VDPoly_r_Rsh_E[id]/1000. <<  setw(7) <<  setprecision(4) << VDPstop_s_Rsh_E[id]/1000.  <<  setw(7)
	       << VDPstop_r_Rsh_E[id]/1000.  <<  setw(7)  << VDPStrip_s_Rsh_E[id] <<  setw(7)  << VDPStrip_r_Rsh_E[id]
	       << setw(9)  <<  setprecision(3) << FET_Gate_Voltage_Threshold[2*id]
	       << setw(10) << setprecision(4) << Cacc[2*id] << setw(9) <<  setprecision(2) << Tox[2*id]
	       << setw(10) << VFB[2*id] << setw(14)  << Noxide[2*id]
	       << endl;


         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
               << setw(8)  <<  setprecision(3) << CAP_Value[4*id+2] << setw(10)  <<  setprecision(3) << CAP_Value[4*id+3]
	       << setw(8) <<  setprecision(3) << VDPoly_s_Rsh_W[id]/1000.  <<  setw(7)
	       << VDPoly_r_Rsh_W[id]/1000. <<  setw(7) <<  setprecision(4) << VDPstop_s_Rsh_W[id]/1000.  <<  setw(7)
	       << VDPstop_r_Rsh_W[id]/1000.  <<  setw(7)  << VDPStrip_s_Rsh_W[id] <<  setw(7)  << VDPStrip_r_Rsh_W[id]
	       << setw(9)  <<  setprecision(3) << FET_Gate_Voltage_Threshold[2*id+1]
	       << setw(10) <<  setprecision(4) << Cacc[2*id+1] <<  setw(9) <<  setprecision(2) << Tox[2*id+1]
	       << setw(10) << VFB[2*id+1] << setw(14)  << Noxide[2*id+1]
	       << endl;
   }

   cout << "----------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   CVS_Output_File<< "     Flute1  Id            ;Mode ;CAP_E  ; CAP_W ; VDPPoly;   ;   VDPStop;  ;  VdPStrip;  ; FET_Vth ; C_acc  ;  t_ox  ;  VFB_Deb  ; Nox_Deb" << endl;
   CVS_Output_File<< "                           ; ; [pF] ;  [pF] ;[kOhm/sq];  ;  [kOhm/sq]; ;  [Ohm/sq];  ;   [V] ;   [pF]  ;  [um]  ;    [V]    ; [cm^{-2}] " <<endl;
   CVS_Output_File<< "                           ; ;      ;       ; stad ; rot ; stad ;  rot ; stad ;  rot ;        ; ; ; ; " << endl;
   CVS_Output_File<< "     Spec Limit            ; ;>2.028 ;>2.028 ; 2.2-2.4 ;2.2-2.4 ;18- 20 ;18 - 20 ;32 - 35 ;32 - 35 ; ;  ~82 ;  0.7 ;  <5 ;  " << endl;

   for(int id=0;id<nFiles;id++) {
         CVS_Output_File  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())<< ";" << "Quick" << ";"
	 	          << CAP_Value[4*id]  << ";"  << CAP_Value[4*id+1] << ";"
                          << VDPoly_s_Rsh_E[id]/1000.   <<  ";" << VDPoly_r_Rsh_E[id]/1000.  << ";"
	                  << VDPstop_s_Rsh_E[id]/1000.  <<  ";" << VDPstop_r_Rsh_E[id]/1000. << ";"
	                  << VDPStrip_s_Rsh_E[id] <<  ";"  << VDPStrip_r_Rsh_E[id] << ";" << FET_Gate_Voltage_Threshold[2*id] << ";"
	                  << Cacc[2*id] << ";"  << Tox[2*id] <<";" <<  VFB[2*id] << ";"  << Noxide[2*id]
                	  << endl;


         CVS_Output_File  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << ";" << "Quick" << ";"
	 	 	  << CAP_Value[4*id+2]  << ";"  << CAP_Value[4*id+3] << ";"
	                  << VDPoly_s_Rsh_W[id]/1000.   <<  ";" << VDPoly_r_Rsh_W[id]/1000.  << ";"
	                  << VDPstop_s_Rsh_W[id]/1000.  <<  ";" << VDPstop_r_Rsh_W[id]/1000. << ";"
	                  << VDPStrip_s_Rsh_W[id] <<  ";"  << VDPStrip_r_Rsh_W[id] << ";" << FET_Gate_Voltage_Threshold[2*id+1]<< ";"
	                  << Cacc[2*id+1] << ";"  << Tox[2*id+1] <<";" <<  VFB[2*id+1] << ";"  << Noxide[2*id+1]
			  << endl;
   }

    CVS_Output_File.close();

}
