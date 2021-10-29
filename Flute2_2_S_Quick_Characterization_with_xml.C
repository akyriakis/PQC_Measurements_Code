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

double epsilon_0;
double epsilon_Si;
double epsilon_SiO2;
double GCD_area;
double q;
double kB;
double TKelvin;
double Nintrisic;

// ---------------------------------------------------

float Diel_IV_vv[100][1000];
float Diel_IV_Curr[100][1000];
double Diel_vv_d[1000][1000];
double Diel_Curr_d[1000][1000];
float Diel_IV_Temp_measurement[1000];
float Diel_IV_Humid_measurement[1000];
TGraph *Diel_IV[500];
float Diel_BreakDown_Voltage[500];
vector< vector <string > > Diel_IV_timestamp_meas;
vector< vector <double > > Diel_IV_temp_meas;
vector< vector <double > > Diel_IV_air_temp_meas;
vector< vector <double > > Diel_IV_rh_prcnt_meas;
vector<std::string> Diel_IV_Operators;
vector<std::string> Diel_IV_Begin_Timestamp;
vector<std::string> Diel_IV_name_labels;
vector<std::string> Diel_IV_kind_of_parts;
vector<std::string> Diel_IV_struct_id;
vector<std::string> Diel_IV_set_id;
vector<std::string> Diel_IV_config_id;
vector<std::string> Diel_IV_equipment;
vector<double> Diel_IV_waiting_time;
vector<double> Diel_IV_bias;
vector<double> Diel_IV_temp;
vector<double> Diel_IV_av_temp;
vector<double> Diel_IV_compliance;
vector<double> Diel_IV_Nmeas;
// -------------------------------------------

float VDPStop_current_d[500];
float VDPStop_voltage_d[500];
float VDPStop_Temp_measurement[600];
float VDPStop_Humid_measurement[600];
double VDPStop_A_constant[600];
double VDPStop_A_slope[600];
float VDPStop_cv_low[600];
float VDPStop_cv_high[600];
double VDPStop_s_Rsh_E[100];
double VDPStop_r_Rsh_E[100];
double VDPStop_s_Rsh_W[100];
double VDPStop_r_Rsh_W[100];
TGraph *VDPStop_iv[600];

float LineWidthStop_current_d[500][500];
float LineWidthStop_voltage_d[500][500];
float LineWidthStop_current[500][500];
float LineWidthStop_voltage[500][500];
float LineWidthStop_Temp_measurement[600];
float LineWidthStop_Humid_measurement[600];
double LineWidthStop_A_constant[600];
double LineWidthStop_A_slope[600];
float LineWidthStop_cv_low[600];
float LineWidthStop_cv_high[600];
double LineWidthStop_Rsh_E[100];
double LineWidthStop_Rsh_W[100];
double LineWidthStop_all[100];
double LineWidthStop_res_all[100];
double LineWidthStop_s_Rsh_Rgeom_Corrected_E[100];
double LineWidthStop_r_Rsh_Rgeom_Corrected_E[100];
double LineWidthStop_s_Rsh_Rgeom_Corrected_W[100];
double LineWidthStop_r_Rsh_Rgeom_Corrected_W[100];
TGraph *LineWidthStop_iv[600];
vector<vector <string > > LineWidthStop_timestamp_meas;
vector<vector <double > > LineWidthStop_temp_meas;
vector<vector <double > > LineWidthStop_air_temp_meas;
vector<vector <double > > LineWidthStop_rh_prcnt_meas;
vector<std::string> LineWidthStop_Operators;
vector<std::string> LineWidthStop_Begin_Timestamp;
vector<std::string> LineWidthStop_name_labels;
vector<std::string> LineWidthStop_kind_of_parts;
vector<std::string> LineWidthStop_kind_of_HM_flute_id;
vector<std::string> LineWidthStop_struct_id;
vector<std::string> LineWidthStop_set_id;
vector<std::string> LineWidthStop_config_id;
vector<std::string> LineWidthStop_equipment;
vector<double> LineWidthStop_waiting_time;
vector<double> LineWidthStop_temp;
vector<double> LineWidthStop_av_temp;
vector<double> LineWidthStop_Nmeas;

// --------------------------------------------------------

float VDPStrip_current_d[500];
float VDPStrip_voltage_d[500];

float VDPStrip_Temp_measurement[600];
float VDPStrip_Humid_measurement[600];

double VDPStrip_A_constant[600];
double VDPStrip_A_slope[600];

float VDPStrip_cv_low[600];
float VDPStrip_cv_high[600];

double VDPStrip_s_Rsh_E[100];
double VDPStrip_r_Rsh_E[100];

double VDPStrip_s_Rsh_W[100];
double VDPStrip_r_Rsh_W[100];

TGraph *VDPStrip_iv[600];


float LineWidthStrip_current_d[500][500];
float LineWidthStrip_voltage_d[500][500];
float LineWidthStrip_current[500][500];
float LineWidthStrip_voltage[500][500];
float LineWidthStrip_Temp_measurement[600];
float LineWidthStrip_Humid_measurement[600];
double LineWidthStrip_A_constant[600];
double LineWidthStrip_A_slope[600];
float LineWidthStrip_cv_low[600];
float LineWidthStrip_cv_high[600];
double LineWidthStrip_Rsh_E[100];
double LineWidthStrip_Rsh_W[100];
double LineWidthStrip_all[100];
double LineWidthStrip_res_all[100];
double LineWidthStrip_s_Rsh_Rgeom_Corrected_E[100];
double LineWidthStrip_r_Rsh_Rgeom_Corrected_E[100];
double LineWidthStrip_s_Rsh_Rgeom_Corrected_W[100];
double LineWidthStrip_r_Rsh_Rgeom_Corrected_W[100];
TGraph *LineWidthStrip_iv[600];
vector<vector <string > > LineWidthStrip_timestamp_meas;
vector<vector <double > > LineWidthStrip_temp_meas;
vector<vector <double > > LineWidthStrip_air_temp_meas;
vector<vector <double > > LineWidthStrip_rh_prcnt_meas;
vector<std::string> LineWidthStrip_Operators;
vector<std::string> LineWidthStrip_Begin_Timestamp;
vector<std::string> LineWidthStrip_name_labels;
vector<std::string> LineWidthStrip_kind_of_parts;
vector<std::string> LineWidthStrip_kind_of_HM_flute_id;
vector<std::string> LineWidthStrip_struct_id;
vector<std::string> LineWidthStrip_set_id;
vector<std::string> LineWidthStrip_config_id;
vector<std::string> LineWidthStrip_equipment;
vector<double> LineWidthStrip_waiting_time;
vector<double> LineWidthStrip_temp;
vector<double> LineWidthStrip_av_temp;
vector<double> LineWidthStrip_Nmeas;

// --------------------------------------------------------

float RPoly_current_d[500][500];
float RPoly_voltage_d[500][500];
float RPoly_current[500][500];
float RPoly_voltage[500][500];
float RPoly_Temp_measurement[600];
float RPoly_Humid_measurement[600];
double RPoly_A_constant[600];
double RPoly_A_slope[600];
float RPoly_cv_low[600];
float RPoly_cv_high[600];
double RPoly_E[100];
double RPoly_r_Rsh_E[100];
double RPoly_W[100];
double RPoly_r_Rsh_W[100];
double RPoly_all[100];
double RPoly_Resistivity_s_E[100];
double RPoly_Resistivity_r_E[100];
double RPoly_Resistivity_s_W[100];
double RPoly_Resistivity_r_W[100];
TGraph *RPoly_iv[600];
vector< vector <string > > RPoly_timestamp_meas;
vector< vector <double > > RPoly_temp_meas;
vector< vector <double > > RPoly_air_temp_meas;
vector< vector <double > > RPoly_rh_prcnt_meas;
vector<std::string> RPoly_Operators;
vector<std::string> RPoly_Begin_Timestamp;
vector<std::string> RPoly_name_labels;
vector<std::string> RPoly_kind_of_parts;
vector<std::string> RPoly_struct_id;
vector<TString> RPoly_set_id;
vector<std::string> RPoly_config_id;
vector<std::string> RPoly_equipment;
vector<double> RPoly_waiting_time;
vector<double> RPoly_bias;
vector<double> RPoly_temp;
vector<double> RPoly_av_temp;
vector<double> RPoly_compliance;
vector<double> RPoly_Nmeas;

// --------------------------------------------------------

float GCD_Temp_measurement[100];
float GCD_Humid_measurement[100];
float GCD_Vbias_measurement[100];
double GCD_A_constant[100];
double GCD_A_slope[100];
double GCD_D_constant[100];
double GCD_D_slope[100];
double GCD_I_constant[100];
double GCD_I_slope[100];
float GCD_iv_low[100];
float GCD_iv_high[100];
float GCD_iv1_low[100];
float GCD_iv1_high[100];
float GCD_iv2_low[100];
float GCD_iv2_high[100];
float GCD_I_accumulation[100];
float GCD_I_depletion[100];
float GCD_I_inversion[100];
float GCD_VFB_Derivative[100];
float GCD_I_surface[100];
float GCD_S_interface_recombination_Velocity[100];
float GCD_D_it[100];
float GCD_N_it[100];
double GCD_vv[500][500];
double GCD_current[500][500];
double GCD_vv_d[500][500];
double GCD_current_d[500][500];
double GCD_vv_d_reduced[500];
double GCD_current_d_reduced[500];
double GCD_dev_spline_simple[500];
TCanvas *GCD_cc_spline[100];
TCanvas *GCD_cc_spline_der[100];
TGraph *GCD_iv[100];
TSpline3 * GCD_CuSpl_iv[100];
vector< vector <string > > GCD_timestamp_meas;
vector< vector <double > > GCD_temp_meas;
vector< vector <double > > GCD_air_temp_meas;
vector< vector <double > > GCD_rh_prcnt_meas;
vector<std::string> GCD_Operators;
vector<std::string> GCD_Begin_Timestamp;
vector<std::string> GCD_name_labels;
vector<std::string> GCD_kind_of_parts;
vector<std::string> GCD_struct_id;
vector<TString> GCD_set_id;
vector<std::string> GCD_config_id;
vector<std::string> GCD_equipment;
vector<double> GCD_waiting_time;
vector<double> GCD_bias;
vector<double> GCD_temp;
vector<double> GCD_av_temp;
vector<double> GCD_compliance;
vector<double> GCD_Nmeas;



ofstream CVS_Output_File;

TString Data_Dir = "/home/akyriakis/MOS-Measurements_Analysis/PQC_Measurements/Data_New_Format/";

TString Diel_IV_DataFileName[100];
TString LineWidthStop_DataFileName[600];
TString LineWidthStrip_DataFileName[600];
TString VDP_DataFileName[100];
TString RPoly_DataFileName[600];
TString GCD_DataFileName[100];

int xml_on = 1;

// ==================================================================================

//int HM_Id = 34355;
//string Structure_Id[100] = {"02","03","04","05","06","10","12","14","16","22","26","30","31","41","42","43","45","46","47","48"};
//int HM_Id = 34356;
//string Structure_Id[100] = {"01","02","03","06","14","15","20","22","23","24","30","31","34","46","47"};
//int HM_Id = 35715;
//string Structure_Id[40] = {"03","16","28","41"};
//int HM_Id = 35716;
//string Structure_Id[40] = {"01","10","22","41","49"};
//int HM_Id = 35720;
//string Structure_Id[40] = {"01","11","18","29","44"};
//int HM_Id = 35721;
//string Structure_Id[40] = {"06","16","25","39","47"};
//int HM_Id = 36243;
//string Structure_Id[40] = {"03","11","20","31","40"};
//int HM_Id = 36244;
//string Structure_Id[40] = {"07","14","23","35","44"};
//int HM_Id = 36244;
//string Structure_Id[40] = {"07","14","23","35","44"};
//int HM_Id = 37400;
//string Structure_Id[40] = {"03","15","27","37"}; 
//int HM_Id = 37401;
//string Structure_Id[40] = {"05","14","28","35"}; 
//int HM_Id = 37408;
//string Structure_Id[40] = {"04","16","27","34"}; 
int HM_Id = 37897;
string Structure_Id[40] = {"01","12","26","39"}; 


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


// ------------------------------------------------------------------------
void Diel_IV_Read_Info_from_File(TString Diel_IV_DataFileName, int id, int file_id)
{

   float Diel_half_area = 0.25*0.25; // in cm^2
   float Diel_width = 0.03; // in cm = 300um
   float Diel_Volume = Diel_half_area*Diel_width; //in cm^3

   float Diel_Current_600V = 0;
   ifstream in;
   in.open(Diel_IV_DataFileName);

   int nlines = 0;
   int N_meas = 0;

   string line;

   char date[20],time[20];
   float Diel_Voltage, Diel_Current, Diel_Temperature, Diel_AirTemperature, Diel_Hymidity;;
   while (in.good()) {
     getline (in,line);
      if(nlines < 18) {
         if(nlines == 1) {
         //cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         Diel_IV_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         Diel_IV_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         Diel_IV_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         Diel_IV_kind_of_parts.push_back(Kind_of_part);
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
         Diel_IV_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         Diel_IV_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         Diel_IV_config_id.push_back(Kind_of_HM_config_id);
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
         Diel_IV_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         Diel_IV_waiting_time.push_back(Waiting_time_s);
         }
         if(nlines == 12) {
         cout << line << '\n';
         double Bias_V;
         sscanf(line.c_str(), "%*s %lf", &Bias_V);
         cout<<Bias_V<<endl;
         Diel_IV_bias.push_back(Bias_V);
         }
         if(nlines == 13) {
         cout << line << '\n';
         double Temp_set_degC;
         sscanf(line.c_str(), "%*s %lf", &Temp_set_degC);
         cout<<Temp_set_degC<<endl;
         Diel_IV_temp.push_back(Temp_set_degC);
        }
        if(nlines == 14) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        Diel_IV_av_temp.push_back(Av_temp_degC);
       }
        if(nlines == 15) {
        cout << line << '\n';
        double ComplianceA;
        sscanf(line.c_str(), "%*s %lf", &ComplianceA);
        cout<<ComplianceA<<endl;
        Diel_IV_compliance.push_back(ComplianceA);
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

    if(Diel_Voltage>0) {
      Diel_vv_d[id][N_meas] = voltage;
      Diel_Curr_d[id][N_meas] = current_namp;
      Diel_IV_vv[id][N_meas] = float(Diel_vv_d[id][N_meas]);
      Diel_IV_Curr[id][N_meas] = float(Diel_Curr_d[id][N_meas]);
      Diel_IV_temp_meas.push_back(std::vector<double>());
      Diel_IV_temp_meas[id].push_back(temp_degC);
      Diel_IV_air_temp_meas.push_back(std::vector<double>());
      Diel_IV_air_temp_meas[id].push_back(air_temp_degC);
      Diel_IV_rh_prcnt_meas.push_back(std::vector<double>());
      Diel_IV_rh_prcnt_meas[id].push_back(rh_prcnt);
      Diel_IV_timestamp_meas.push_back(std::vector<string>());
      Diel_IV_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
        N_meas++;
       }
      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   Diel_IV_Nmeas.push_back(N_meas);
   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;


   //TCanvas *cc = new TCanvas("cc","Diode Current vs Diel_Voltage",50,50,1000,1000);
   Diel_IV[id] =  new TGraph(N_meas, Diel_vv_d[id], Diel_Curr_d[id]);

   if((id % 6) == 0) Diel_IV[id]->SetTitle(Form("Diel_NE_E%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 6) == 1) Diel_IV[id]->SetTitle(Form("Diel_NW_E%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 6) == 2) Diel_IV[id]->SetTitle(Form("Diel_SW_E%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 6) == 3) Diel_IV[id]->SetTitle(Form("Diel_NE_W%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 6) == 4) Diel_IV[id]->SetTitle(Form("Diel_NW_W%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 6) == 5) Diel_IV[id]->SetTitle(Form("Diel_SW_W%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));

   Diel_IV[id]->GetXaxis()->SetTitle("Diel_Voltage [V]");
   Diel_IV[id]->GetYaxis()->SetTitle("Diel_Current [nA]");
   Diel_IV[id]->SetDrawOption("AP");
   Diel_IV[id]->SetMarkerStyle(20+id);
   int colorid = 1+id;
   if(colorid == 5) colorid = 46;
   Diel_IV[id]->SetMarkerColor(colorid);
   //Diel_IV[id]->Draw();

   Diel_BreakDown_Voltage[id] = Diel_vv_d[id][N_meas-1];
   cout << "----->>>  Diel_BreakDown_Voltage = " <<  Diel_BreakDown_Voltage[id] << " V"<<  endl;

}
void Diel_IV_PQC2_Analysis_Final(int nFiles)
{

 //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1

 //   gROOT->Reset();
   gStyle->SetOptStat(0);

   //TString Diel_IV_DataFileName[500];
   for(int i=0;i<nFiles;i++) {
    Diel_IV_DataFileName[6*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Left_PQC2_DielNE_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    Diel_IV_DataFileName[6*i+1] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Left_PQC2_DielNW_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    Diel_IV_DataFileName[6*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Left_PQC2_DielSW_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());


    Diel_IV_DataFileName[6*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Left_PQC2_DielNE_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    Diel_IV_DataFileName[6*i+4] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Left_PQC2_DielNW_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    Diel_IV_DataFileName[6*i+5] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Left_PQC2_DielSW_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {
     cout << endl;
     cout << "Analyze file : " << Diel_IV_DataFileName[6*i] << endl;
     Diel_IV_Read_Info_from_File(Diel_IV_DataFileName[6*i],6*i,i);

     cout << "Analyze file : " << Diel_IV_DataFileName[6*i+1] << endl;
     Diel_IV_Read_Info_from_File(Diel_IV_DataFileName[6*i+1],6*i+1,i);

     cout << "Analyze file : " << Diel_IV_DataFileName[6*i+2] << endl;
     Diel_IV_Read_Info_from_File(Diel_IV_DataFileName[6*i+2],6*i+2,i);

     cout << "Analyze file : " << Diel_IV_DataFileName[6*i+3] << endl;
     Diel_IV_Read_Info_from_File(Diel_IV_DataFileName[6*i+3],6*i+3,i);


     cout << "Analyze file : " << Diel_IV_DataFileName[6*i+4] << endl;
     Diel_IV_Read_Info_from_File(Diel_IV_DataFileName[6*i+4],6*i+4,i);

     cout << "Analyze file : " << Diel_IV_DataFileName[6*i+5] << endl;
     Diel_IV_Read_Info_from_File(Diel_IV_DataFileName[6*i+5],6*i+5,i);

  }

   TCanvas *Diel_cc_final[6];
   TMultiGraph *mg[6];
   for (int i=0;i<6;i++) {
     Diel_cc_final[i] = new TCanvas(Form("cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     mg[i] = new TMultiGraph(Form("mg_%d",i),"Diode  Current vs V");
     if(i == 0) mg[i]->SetTitle("Diel_NE E Current vs V");
     if(i == 1) mg[i]->SetTitle("Diel_NW E Current vs V");
     if(i == 2) mg[i]->SetTitle("Diel_SW E Current vs V");
     if(i == 3) mg[i]->SetTitle("Diel_NE W Current vs V");
     if(i == 4) mg[i]->SetTitle("Diel_NW W Current vs V");
     if(i == 5) mg[i]->SetTitle("Diel_SW W Current vs V");

     for(int id=0;id<nFiles;id++) {

        mg[i]->Add(Diel_IV[6*id+i]);

     }
     mg[i]->Draw("LPsame");
     mg[i]->GetXaxis()->SetTitle("Diel_Voltage [V]");
     mg[i]->GetYaxis()->SetTitle("I [nA]");
     mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     Diel_cc_final[i]->BuildLegend(0.15,0.45,0.55,0.85);


     gPad->Modified();
     gPad->Update();
   }

   cout <<endl;
   cout << "===========================+=================================================================================" << endl;
   cout << "---------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  ----------------------------------------------" << endl;
   cout << "=============================================================================================================" << endl;

   cout << "       PQC2  Diel Id        |  DielNE_BreakVoltage [V] |  DielNW_BreakVoltage [V] |  DielSW_BreakVoltage [V] " <<endl;
   cout << "-------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
      cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
            << setw(17)  <<  setprecision(3) <<  Diel_BreakDown_Voltage[6*id]
            << setw(27)  <<  setprecision(3) <<  Diel_BreakDown_Voltage[6*id+1]
            << setw(27)  <<  setprecision(3) <<  Diel_BreakDown_Voltage[6*id+2]
	    << endl;


      cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
            << setw(17)  <<  setprecision(3) <<  Diel_BreakDown_Voltage[6*id+3]
	    << setw(27)  <<  setprecision(3) <<  Diel_BreakDown_Voltage[6*id+4]
            << setw(27)  <<  setprecision(3) <<  Diel_BreakDown_Voltage[6*id+5]
	    << endl;

     }
   cout << "----------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------------------------------------------//
// Diel_IV xml production
// --------------------------------------------------------------------------------------------------//
void Diel_IV_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    cout << "Diel IV Xml production is on" << endl;
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<6; j++){
  //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = Diel_IV_DataFileName[6*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = Diel_IV_DataFileName[6*i+j].Last(delimeter2);
         xml_filename = Diel_IV_DataFileName[6*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
 //        cout << MOS_DataFileName[6*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);

         // First create engine
        TXMLEngine Diel_IV_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = Diel_IV_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = Diel_IV_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = Diel_IV_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = Diel_IV_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = Diel_IV_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = Diel_IV_xml.NewChild(HEADER, 0, "RUN");
              Diel_IV_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              Diel_IV_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              Diel_IV_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              Diel_IV_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              Diel_IV_xml.NewChild(RUN, 0, "INITIATED_BY_USER", Diel_IV_Operators[6*i+j].c_str());
              Diel_IV_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", Diel_IV_Begin_Timestamp[6*i+j].c_str());
              Diel_IV_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = Diel_IV_xml.NewChild(ROOT, 0, "DATA_SET");
            Diel_IV_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            Diel_IV_xml.NewChild(DATA_SET, 0, "VERSION", "IV_measurement");
            XMLNodePointer_t PART = Diel_IV_xml.NewChild(DATA_SET, 0, "PART");
              Diel_IV_xml.NewChild(PART, 0, "NAME_LABEL", Diel_IV_name_labels[6*i+j].c_str());
              Diel_IV_xml.NewChild(PART, 0, "KIND_OF_PART", Diel_IV_kind_of_parts[6*i+j].c_str());
            XMLNodePointer_t DATA = Diel_IV_xml.NewChild(DATA_SET, 0, "DATA");
              Diel_IV_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC2");
              Diel_IV_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", Diel_IV_struct_id[6*i+j].c_str());
              Diel_IV_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", Diel_IV_set_id[6*i+j].c_str());
              Diel_IV_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", Diel_IV_config_id[6*i+j].c_str());
              Diel_IV_xml.NewChild(DATA, 0, "EQUIPMENT", Diel_IV_equipment[6*i+j].c_str());
              Diel_IV_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(Diel_IV_waiting_time[6*i+j],3).c_str());
              Diel_IV_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(Diel_IV_temp[6*i+j],3).c_str());
              Diel_IV_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(Diel_IV_av_temp[6*i+j],3).c_str());
              Diel_IV_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "Diel_IV1");
              Diel_IV_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = Diel_IV_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = Diel_IV_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = Diel_IV_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  Diel_IV_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  Diel_IV_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = Diel_IV_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                Diel_IV_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                Diel_IV_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "IV_measurement");
                XMLNodePointer_t PART_CDS1 = Diel_IV_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  Diel_IV_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", Diel_IV_name_labels[6*i+j].c_str());
                  Diel_IV_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", Diel_IV_kind_of_parts[6*i+j].c_str());
                for (int k=0; k<Diel_IV_Nmeas[6*i+j]-1; k++){
                XMLNodePointer_t DATA_CDS1 = Diel_IV_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  Diel_IV_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(Diel_IV_vv[6*i+j][k],2).c_str());
                  //Diel_IV_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", to_string(Diel_IV_Curr[k]).c_str());
                  Diel_IV_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string(Diel_IV_Curr[6*i+j][k],4).c_str());
                  Diel_IV_xml.NewChild(DATA_CDS1, 0, "TIME", Diel_IV_timestamp_meas[6*i+j][k].c_str());
                  Diel_IV_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string_scientific_format(Diel_IV_temp_meas[6*i+j][k], 3).c_str());
                  Diel_IV_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(Diel_IV_air_temp_meas[6*i+j][k], 3).c_str());
                  Diel_IV_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(Diel_IV_rh_prcnt_meas[6*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = Diel_IV_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = Diel_IV_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = Diel_IV_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      Diel_IV_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      Diel_IV_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = Diel_IV_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    Diel_IV_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    Diel_IV_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement");
                    XMLNodePointer_t PART_CDS2 = Diel_IV_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      Diel_IV_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", Diel_IV_name_labels[6*i+j].c_str());
                      Diel_IV_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", Diel_IV_kind_of_parts[6*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = Diel_IV_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      Diel_IV_xml.NewChild(DATA_CDS2, 0, "VBD_V", Conv_float_to_string(Diel_BreakDown_Voltage[6*i+j], 3).c_str());
        XMLDocPointer_t Diel_IV_xmldoc = Diel_IV_xml.NewDoc();
        Diel_IV_xml.DocSetRootElement(Diel_IV_xmldoc, ROOT);
        // Save document to file
        Diel_IV_xml.SaveDoc(Diel_IV_xmldoc, Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        Diel_IV_xml.FreeDoc(Diel_IV_xmldoc);
      }
    }
  }
}
// ==================================================================================================//
// ------------------------------------------------------------------------
double VDPStop_Read_Info_from_File(TString VDPStop_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(VDPStop_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << VDPStop_DataFileName << " wit Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;

    char date[20],time[20];
   float Voltage, Current, VdP_Temperature, VdP_AirTemperature, VdP_Hymidity;
   while (in.good()) {

     if(nlines < 16) {
         getline (in,line);
         if(nlines == 13) {
	   cout << "Line Number  = " << nlines << " with content : " <<  line << '\n';
	   string Temperature = line.substr(14, 7);
	   VDPStop_Temp_measurement[id] = std::stof(Temperature);
	   cout << "Temperature = " << VDPStop_Temp_measurement[id] << '\n';
        }
	} else {

          in  >> date >> time >> Current >> Voltage >> VdP_Temperature >> VdP_AirTemperature >> VdP_Hymidity;

	VDPStop_current_d[N_meas] = Current;
	VDPStop_voltage_d[N_meas] = Voltage*1E-9;

        //cout << N_meas << ")  "  <<  current_d[N_meas] << "   " << voltage_d[N_meas] << endl;
	N_meas++;
      }

      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();

   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;

   TGaxis::SetMaxDigits(3);

//   TCanvas *cc_spline = new TCanvas();
   VDPStop_iv[id] =  new TGraph(N_meas,VDPStop_current_d, VDPStop_voltage_d);
   if((id % 4) == 0) VDPStop_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStop_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 1) VDPStop_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStop_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 2) VDPStop_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStop_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 3) VDPStop_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStop_r",HM_Id,Structure_Id[file_id].c_str()));
   VDPStop_iv[id]->GetXaxis()->SetTitle("VdP Current [A]");
   VDPStop_iv[id]->GetYaxis()->SetTitle("VdP Voltage [V]");
   VDPStop_iv[id]->SetDrawOption("AP");
   VDPStop_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   VDPStop_iv[id]->SetMarkerColor(colorid);
//   VDPStop_iv[id]->Draw();


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   VDPStop_cv_low[id] = VDPStop_current_d[3];
   VDPStop_cv_high[id] = VDPStop_current_d[38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", VDPStop_cv_low[id], VDPStop_cv_high[id]);
//   VDPStop_iv[id]->Fit("f_cv","0R+");
   VDPStop_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   VDPStop_A_constant[id] = VDPStop_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   VDPStop_A_slope[id]  = VDPStop_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double Resistance = (TMath::Pi()/TMath::Log(2))*VDPStop_A_slope[id]; // VdP Resistance

   cout << "----->>> VDPStop PQC1 Structure Resistance  = " << Resistance << " [Ohm/sq]" << endl;

   return Resistance;

}
void VDPStop_PQC_Flute1_Analysis_Final(int nFiles)
{

   gStyle->SetOptStat(0);

   TString VDPStop_DataFileName[600];
   for(int i=0;i<nFiles;i++) {

    VDPStop_DataFileName[4*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStop_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDPStop_DataFileName[4*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStop_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDPStop_DataFileName[4*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStop_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDPStop_DataFileName[4*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStop_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {


     cout << "VDPStop Flute1 Structure: Analyze file : " << VDPStop_DataFileName[4*i] << endl;
     VDPStop_s_Rsh_E[i] = VDPStop_Read_Info_from_File(VDPStop_DataFileName[4*i],(4*i),i);


     cout << "VDPStop Flute1 Structure: Analyze file : " << VDPStop_DataFileName[4*i+1] << endl;
     VDPStop_s_Rsh_W[i] = VDPStop_Read_Info_from_File(VDPStop_DataFileName[4*i+1],(4*i+1),i);

     cout << "VDPStop Flute1 Structure: Analyze file : " << VDPStop_DataFileName[4*i+2] << endl;
     VDPStop_r_Rsh_E[i] = VDPStop_Read_Info_from_File(VDPStop_DataFileName[4*i+2],(4*i+2),i);

     cout << "VDPStop Flute1 Structure: Analyze file : " << VDPStop_DataFileName[4*i+3] << endl;
     VDPStop_r_Rsh_W[i] = VDPStop_Read_Info_from_File(VDPStop_DataFileName[4*i+3],(4*i+3),i);
   }


   cout <<endl;
   cout << "======================================================================================================================" << endl;
   cout << "---------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -----------------" << endl;
   cout << "==================================================================================================" << endl;

   cout << "      VDPStop   Id                  | VDPStop_s_Rsh [Ohm/sq] | VDPStop_r_Rsh [Ohm/sq]" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) <<  setw(25)  <<  setprecision(4)<< VDPStop_s_Rsh_E[id] <<  setw(25)  << VDPStop_r_Rsh_E[id]   << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) <<  setw(25) <<  setprecision(4) << VDPStop_s_Rsh_W[id] <<  setw(25)  << VDPStop_r_Rsh_W[id]   << endl;
   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;


}

// ------------------------------------------------------------------------
double LineWidthStop_Read_Info_from_File(TString LineWidthStop_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(LineWidthStop_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << LineWidthStop_DataFileName << " wit Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;
   while (in.good()) {
     getline (in,line);
      if(nlines < 16) {
         if(nlines == 1) {
         cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         LineWidthStop_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         LineWidthStop_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         LineWidthStop_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         LineWidthStop_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         LineWidthStop_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         LineWidthStop_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         LineWidthStop_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         LineWidthStop_config_id.push_back(Kind_of_HM_config_id);
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
         LineWidthStop_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         LineWidthStop_waiting_time.push_back(Waiting_time_s);
         }
        if(nlines == 12) {
        cout << line << '\n';
        double temp_degC;
        sscanf(line.c_str(), "%*s %lf", &temp_degC);
        cout<<temp_degC<<endl;
        LineWidthStop_temp.push_back(temp_degC);
       }
        if(nlines == 13) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        LineWidthStop_av_temp.push_back(Av_temp_degC);
       }
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double current_namp;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf ", date, time, &current_namp, &voltage, &temp_degC, &air_temp_degC, &rh_prcnt);
      //cout<<date<<","<<time<<","<<voltage<<","<<(float) current_namp<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
   //        in >>  LineWidthStop_voltage >> LineWidthStop_capacitance >> LineWidthStop_conductance;
      //  in >>  LineWidthStop_voltage >> LineWidthStop_capacitance;
      LineWidthStop_voltage_d[id][N_meas] = voltage*1E-9;  // to 1E-9 exei mpei gia na diorthwsei to lathos pou yparxei sta txt. tha prepei na aferethei an diorthwthei
      LineWidthStop_current_d[id][N_meas] = current_namp;
      LineWidthStop_temp_meas.push_back(std::vector<double>());
      LineWidthStop_temp_meas[id].push_back(temp_degC);
      LineWidthStop_air_temp_meas.push_back(std::vector<double>());
      LineWidthStop_air_temp_meas[id].push_back(air_temp_degC);
      LineWidthStop_rh_prcnt_meas.push_back(std::vector<double>());
      LineWidthStop_rh_prcnt_meas[id].push_back(rh_prcnt);
      LineWidthStop_timestamp_meas.push_back(std::vector<string>());
      LineWidthStop_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      LineWidthStop_voltage[id][N_meas] = float(LineWidthStop_voltage_d[id][N_meas]);
      LineWidthStop_current[id][N_meas] = float(LineWidthStop_current_d[id][N_meas]);
    //        cout << N_meas << ")  "  <<  LineWidthStop_vv[N_meas] << "   " << LineWidthStop_cap[N_meas] << endl;
      cout << N_meas << ")  "  <<  LineWidthStop_voltage[N_meas] << "   " << LineWidthStop_current[N_meas]  << endl;
      N_meas++;

      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   LineWidthStop_Nmeas.push_back(N_meas);
   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;

   TGaxis::SetMaxDigits(3);

//   TCanvas *cc_spline = new TCanvas();
   LineWidthStop_iv[id] =  new TGraph(N_meas,LineWidthStop_current_d[id], LineWidthStop_voltage_d[id]);
   if((id % 2) == 0) LineWidthStop_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC2_LineWidthStop",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) LineWidthStop_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC2_LineWidthStop",HM_Id,Structure_Id[file_id].c_str()));
   LineWidthStop_iv[id]->GetXaxis()->SetTitle("LineWidthStop Current [A]");
   LineWidthStop_iv[id]->GetYaxis()->SetTitle("LineWidthStop Voltage [V]");
   LineWidthStop_iv[id]->SetDrawOption("AP");
   LineWidthStop_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   LineWidthStop_iv[id]->SetMarkerColor(colorid);
//   LineWidthStop_iv[id]->Draw();


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   LineWidthStop_cv_low[id] = LineWidthStop_current_d[id][3];
   LineWidthStop_cv_high[id] = LineWidthStop_current_d[id][38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", LineWidthStop_cv_low[id], LineWidthStop_cv_high[id]);
//   LineWidthStop_iv[id]->Fit("f_cv","0R+");
   LineWidthStop_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   LineWidthStop_A_constant[id] = LineWidthStop_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   LineWidthStop_A_slope[id]  = LineWidthStop_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double R_sh = LineWidthStop_A_slope[id]; // LineWidthStop R_sh

   cout << "----->>>   R_sh Voltage = " << R_sh << " [Ohm/sq]" << endl;

   return R_sh;

}
void LineWidthStop_PQC_Flute2_Analysis_Final(int nFiles)
{
   gStyle->SetOptStat(0);

   VDPStop_PQC_Flute1_Analysis_Final(nFiles); // get the R_Geometry

   for(int i=0;i<nFiles;i++) {
    LineWidthStop_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC2_LinewidthStop.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << LineWidthStop_DataFileName[2*i] << endl;

    LineWidthStop_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC2_LinewidthStop.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << LineWidthStop_DataFileName[2*i+1] << endl;

   }

   float d = 128.5; // in um

   for (int i=0;i<nFiles;i++) {

     LineWidthStop_Rsh_E[i] = LineWidthStop_Read_Info_from_File(LineWidthStop_DataFileName[2*i],(2*i),i);
     LineWidthStop_s_Rsh_Rgeom_Corrected_E[i] = (1/LineWidthStop_Rsh_E[i])*VDPStop_s_Rsh_E[i]*d;
     LineWidthStop_r_Rsh_Rgeom_Corrected_E[i] = (1/LineWidthStop_Rsh_E[i])*VDPStop_r_Rsh_E[i]*d;

     LineWidthStop_Rsh_W[i] = LineWidthStop_Read_Info_from_File(LineWidthStop_DataFileName[2*i+1],(2*i+1),i);
     LineWidthStop_s_Rsh_Rgeom_Corrected_W[i] = (1/LineWidthStop_Rsh_W[i])*VDPStop_s_Rsh_W[i]*d;
     LineWidthStop_r_Rsh_Rgeom_Corrected_W[i] = (1/LineWidthStop_Rsh_W[i])*VDPStop_r_Rsh_W[i]*d;

     LineWidthStop_all[2*i+0] = LineWidthStop_s_Rsh_Rgeom_Corrected_E[i]; // put all std values in one array; for the xml production only std needed
     LineWidthStop_all[2*i+1] = LineWidthStop_s_Rsh_Rgeom_Corrected_W[i];
     LineWidthStop_res_all[2*i+0] = LineWidthStop_Rsh_E[i]; // put all std values in one array; for the xml production only std needed
     LineWidthStop_res_all[2*i+1] = LineWidthStop_Rsh_W[i];
   }


   TCanvas *LineWidthStop_cc_final[4];
   TMultiGraph *LineWidthStop_mg[4];
   for (int i=0;i<2;i++) {
     LineWidthStop_cc_final[i] = new TCanvas(Form("LineWidthStop_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     LineWidthStop_mg[i] = new TMultiGraph(Form("LineWidthStop_mg_%d",i),"LineWidthStop");
     if(i == 0) LineWidthStop_mg[i]->SetTitle("LineWidthStop_E");
     if(i == 1) LineWidthStop_mg[i]->SetTitle("LineWidthStop_W");

    // LineWidthStop_iv[4*i+]->Draw("ALP");
     for(int id=0;id<nFiles;id++) {

        LineWidthStop_mg[i]->Add(LineWidthStop_iv[2*id+i]);

     }
     LineWidthStop_mg[i]->Draw("LPsame");
     LineWidthStop_mg[i]->GetXaxis()->SetTitle("Current [A]");
     LineWidthStop_mg[i]->GetYaxis()->SetTitle("Voltage [V]");
     LineWidthStop_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     LineWidthStop_cc_final[i]->BuildLegend(0.15,0.5,0.5,0.85);


     gPad->Modified();
     gPad->Update();
   }


   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;

   cout << "      LineWidthStop   Id    |  LineWidthStop_s_Rsh [um]   | LineWidthStop_r_Rsh [um]|" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(18) <<  setprecision(4) << LineWidthStop_s_Rsh_Rgeom_Corrected_E[id]  <<  setw(30)  << LineWidthStop_r_Rsh_Rgeom_Corrected_E[id]
	       << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(18) <<  setprecision(4) << LineWidthStop_s_Rsh_Rgeom_Corrected_W[id]  <<  setw(30)  << LineWidthStop_r_Rsh_Rgeom_Corrected_W[id]
	       << endl;
   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

}

void LineWidthStop_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  cout << "LineWidthStop Xml production is on" << endl;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<=1; j++){
          //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = LineWidthStop_DataFileName[2*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = LineWidthStop_DataFileName[2*i+j].Last(delimeter2);
         xml_filename = LineWidthStop_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[4*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);

        // First create engine
        TXMLEngine LineWidthStop_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = LineWidthStop_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = LineWidthStop_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = LineWidthStop_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = LineWidthStop_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = LineWidthStop_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = LineWidthStop_xml.NewChild(HEADER, 0, "RUN");
              LineWidthStop_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              LineWidthStop_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              LineWidthStop_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              LineWidthStop_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              LineWidthStop_xml.NewChild(RUN, 0, "INITIATED_BY_USER", LineWidthStop_Operators[2*i+j].c_str());
              LineWidthStop_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", LineWidthStop_Begin_Timestamp[2*i+j].c_str());
              LineWidthStop_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = LineWidthStop_xml.NewChild(ROOT, 0, "DATA_SET");
            LineWidthStop_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            LineWidthStop_xml.NewChild(DATA_SET, 0, "VERSION", "IV_measurement");
            XMLNodePointer_t PART = LineWidthStop_xml.NewChild(DATA_SET, 0, "PART");
              LineWidthStop_xml.NewChild(PART, 0, "NAME_LABEL", LineWidthStop_name_labels[2*i+j].c_str());
              LineWidthStop_xml.NewChild(PART, 0, "KIND_OF_PART", LineWidthStop_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = LineWidthStop_xml.NewChild(DATA_SET, 0, "DATA");
              LineWidthStop_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC2");
              LineWidthStop_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", LineWidthStop_struct_id[2*i+j].c_str());
              LineWidthStop_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", LineWidthStop_set_id[2*i+j].c_str());
              LineWidthStop_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", LineWidthStop_config_id[2*i+j].c_str());
              LineWidthStop_xml.NewChild(DATA, 0, "EQUIPMENT", LineWidthStop_equipment[2*i+j].c_str());
              LineWidthStop_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(LineWidthStop_waiting_time[2*i+j], 3).c_str());
              LineWidthStop_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(LineWidthStop_temp[2*i+j], 3).c_str());
              LineWidthStop_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(LineWidthStop_av_temp[2*i+j], 3).c_str());
              LineWidthStop_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", LineWidthStop_struct_id[2*i+j].c_str());
              LineWidthStop_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = LineWidthStop_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = LineWidthStop_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = LineWidthStop_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  LineWidthStop_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  LineWidthStop_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = LineWidthStop_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                LineWidthStop_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                LineWidthStop_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "IV_measurement");
                XMLNodePointer_t PART_CDS1 = LineWidthStop_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  LineWidthStop_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", LineWidthStop_name_labels[2*i+j].c_str());
                  LineWidthStop_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", LineWidthStop_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<LineWidthStop_Nmeas[2*i+j]-1; k++){
                    XMLNodePointer_t DATA_CDS1 = LineWidthStop_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                      LineWidthStop_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string_scientific_format(LineWidthStop_voltage[2*i+j][k], 3).c_str());
                      LineWidthStop_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(LineWidthStop_current[2*i+j][k]*1E+9, 3).c_str());
                      LineWidthStop_xml.NewChild(DATA_CDS1, 0, "TIME", LineWidthStop_timestamp_meas[2*i+j][k].c_str());
                      LineWidthStop_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(LineWidthStop_temp_meas[2*i+j][k], 3).c_str());
                      LineWidthStop_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(LineWidthStop_air_temp_meas[2*i+j][k], 3).c_str());
                      LineWidthStop_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(LineWidthStop_rh_prcnt_meas[2*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = LineWidthStop_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = LineWidthStop_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = LineWidthStop_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      LineWidthStop_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      LineWidthStop_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = LineWidthStop_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    LineWidthStop_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    LineWidthStop_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement");
                    XMLNodePointer_t PART_CDS2 = LineWidthStop_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      LineWidthStop_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", LineWidthStop_name_labels[2*i+j].c_str());
                      LineWidthStop_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", LineWidthStop_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = LineWidthStop_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      LineWidthStop_xml.NewChild(DATA_CDS2, 0, "T_UM", Conv_float_to_string(LineWidthStop_all[2*i+j], 3).c_str());
                      LineWidthStop_xml.NewChild(DATA_CDS2, 0, "R_OHM", Conv_float_to_string_scientific_format(LineWidthStop_res_all[2*i+j], 3).c_str());

        XMLDocPointer_t LineWidthStop_xmldoc = LineWidthStop_xml.NewDoc();
        LineWidthStop_xml.DocSetRootElement(LineWidthStop_xmldoc, ROOT);
        // Save document to file
        LineWidthStop_xml.SaveDoc(LineWidthStop_xmldoc, Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        LineWidthStop_xml.FreeDoc(LineWidthStop_xmldoc);
      }
    }
  }
}
// ==================================================================================================//
// ------------------------------------------------------------------------
double VDPStrip_Read_Info_from_File(TString VDPStrip_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(VDPStrip_DataFileName);

   cout << endl << "=============>>>>>> Analyze File : " << VDPStrip_DataFileName << " wit Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;

   char date[20],time[20];
   float Voltage, Current, VdP_Temperature, VdP_AirTemperature, VdP_Hymidity;
   while (in.good()) {

     if(nlines < 16) {
         getline (in,line);
         if(nlines == 13) {
	   cout << "Line Number  = " << nlines << " with content : " <<  line << '\n';
	   string Temperature = line.substr(14, 7);
	   VDPStrip_Temp_measurement[id] = std::stof(Temperature);
	   cout << "Temperature = " << VDPStrip_Temp_measurement[id] << '\n';
        }
	} else {

          in  >> date >> time >> Current >> Voltage >> VdP_Temperature >> VdP_AirTemperature >> VdP_Hymidity;

	VDPStrip_current_d[N_meas] = Current;
	VDPStrip_voltage_d[N_meas] = Voltage*1E-9;

        //cout << N_meas << ")  "  <<  current_d[N_meas] << "   " << voltage_d[N_meas] << endl;
	N_meas++;
      }

      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();

   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;

   TGaxis::SetMaxDigits(3);

//   TCanvas *cc_spline = new TCanvas();
   VDPStrip_iv[id] =  new TGraph(N_meas,VDPStrip_current_d, VDPStrip_voltage_d);
   if((id % 4) == 0) VDPStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStrip_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 1) VDPStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStrip_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 2) VDPStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStrip_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 3) VDPStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStrip_r",HM_Id,Structure_Id[file_id].c_str()));
   VDPStrip_iv[id]->GetXaxis()->SetTitle("VdP Current [A]");
   VDPStrip_iv[id]->GetYaxis()->SetTitle("VdP Voltage [V]");
   VDPStrip_iv[id]->SetDrawOption("AP");
   VDPStrip_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   VDPStrip_iv[id]->SetMarkerColor(colorid);
//   VDPStrip_iv[id]->Draw();


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   VDPStrip_cv_low[id] = VDPStrip_current_d[3];
   VDPStrip_cv_high[id] = VDPStrip_current_d[38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", VDPStrip_cv_low[id], VDPStrip_cv_high[id]);
//   VDPStrip_iv[id]->Fit("f_cv","0R+");
   VDPStrip_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   VDPStrip_A_constant[id] = VDPStrip_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   VDPStrip_A_slope[id]  = VDPStrip_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double Resistance = (TMath::Pi()/TMath::Log(2))*VDPStrip_A_slope[id]; // VdP Resistance

   cout << "----->>> VDPStrip PQC1 Structure Resistance  = " << Resistance << " [Ohm/sq]" << endl;

   return Resistance;

}
void VDPStrip_PQC_Flute1_Analysis_Final(int nFiles)
{

   gStyle->SetOptStat(0);


   TString VDPStrip_DataFileName[600];
   for(int i=0;i<nFiles;i++) {

    VDPStrip_DataFileName[4*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStrip_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDPStrip_DataFileName[4*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStrip_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDPStrip_DataFileName[4*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPStrip_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDPStrip_DataFileName[4*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPStrip_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {
     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPStrip_DataFileName[4*i] << endl;
     VDPStrip_s_Rsh_E[i] = VDPStrip_Read_Info_from_File(VDPStrip_DataFileName[4*i],(4*i),i);
     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPStrip_DataFileName[4*i+1] << endl;
     VDPStrip_s_Rsh_W[i] = VDPStrip_Read_Info_from_File(VDPStrip_DataFileName[4*i+1],(4*i+1),i);

     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPStrip_DataFileName[4*i+2] << endl;
     VDPStrip_r_Rsh_E[i] = VDPStrip_Read_Info_from_File(VDPStrip_DataFileName[4*i+2],(4*i+2),i);

     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPStrip_DataFileName[4*i+3] << endl;
     VDPStrip_r_Rsh_W[i] = VDPStrip_Read_Info_from_File(VDPStrip_DataFileName[4*i+3],(4*i+3),i);
   }


   cout <<endl;
   cout << "======================================================================================================================" << endl;
   cout << "---------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -----------------" << endl;
   cout << "==================================================================================================" << endl;

   cout << "      VDPStrip   Id                  | VDPStrip_s_Rsh [Ohm/sq] | VDPStrip_r_Rsh [Ohm/sq]" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) <<  setw(25)  <<  setprecision(4)<< VDPStrip_s_Rsh_E[id] <<  setw(25)  << VDPStrip_r_Rsh_E[id]   << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) <<  setw(25) <<  setprecision(4) << VDPStrip_s_Rsh_W[id] <<  setw(25)  << VDPStrip_r_Rsh_W[id]   << endl;
   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;


}
double LineWidthStrip_Read_Info_from_File(TString LineWidthStrip_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(LineWidthStrip_DataFileName);

   cout << endl << "=============>>>>>> Analyze File : " << LineWidthStrip_DataFileName << " with Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;
   while (in.good()) {
     getline (in,line);
      if(nlines < 16) {
         if(nlines == 1) {
         cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         LineWidthStrip_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         LineWidthStrip_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         LineWidthStrip_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         LineWidthStrip_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         LineWidthStrip_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         LineWidthStrip_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         LineWidthStrip_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         LineWidthStrip_config_id.push_back(Kind_of_HM_config_id);
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
         LineWidthStrip_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         LineWidthStrip_waiting_time.push_back(Waiting_time_s);
         }
        if(nlines == 12) {
        cout << line << '\n';
        double temp_degC;
        sscanf(line.c_str(), "%*s %lf", &temp_degC);
        cout<<temp_degC<<endl;
        LineWidthStrip_temp.push_back(temp_degC);
       }
        if(nlines == 13) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        LineWidthStrip_av_temp.push_back(Av_temp_degC);
       }
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double current_namp;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf ", date, time, &current_namp, &voltage, &temp_degC, &air_temp_degC, &rh_prcnt);
      //cout<<date<<","<<time<<","<<voltage<<","<<(float) current_namp<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
   //        in >>  LineWidthStrip_voltage >> LineWidthStrip_capacitance >> LineWidthStrip_conductance;
      //  in >>  LineWidthStrip_voltage >> LineWidthStrip_capacitance;
      LineWidthStrip_voltage_d[id][N_meas] = voltage*1E-9;  // to 1E-9 exei mpei gia na diorthwsei to lathos pou yparxei sta txt. tha prepei na aferethei an diorthwthei
      LineWidthStrip_current_d[id][N_meas] = current_namp;
      LineWidthStrip_temp_meas.push_back(std::vector<double>());
      LineWidthStrip_temp_meas[id].push_back(temp_degC);
      LineWidthStrip_air_temp_meas.push_back(std::vector<double>());
      LineWidthStrip_air_temp_meas[id].push_back(air_temp_degC);
      LineWidthStrip_rh_prcnt_meas.push_back(std::vector<double>());
      LineWidthStrip_rh_prcnt_meas[id].push_back(rh_prcnt);
      LineWidthStrip_timestamp_meas.push_back(std::vector<string>());
      LineWidthStrip_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      LineWidthStrip_voltage[id][N_meas] = float(LineWidthStrip_voltage_d[id][N_meas]);
      LineWidthStrip_current[id][N_meas] = float(LineWidthStrip_current_d[id][N_meas]);
    //        cout << N_meas << ")  "  <<  LineWidthStrip_vv[N_meas] << "   " << LineWidthStrip_cap[N_meas] << endl;
    //  cout << N_meas << ")  "  <<  LineWidthStrip_voltage[N_meas] << "   " << LineWidthStrip_current[N_meas]  << endl;
      N_meas++;

      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   LineWidthStrip_Nmeas.push_back(N_meas);
   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;

   TGaxis::SetMaxDigits(3);

   //   TCanvas *cc_spline = new TCanvas();
   LineWidthStrip_iv[id] =  new TGraph(N_meas,LineWidthStrip_current_d[id], LineWidthStrip_voltage_d[id]);
   if((id % 2) == 0) LineWidthStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC2_LineWidthStrip",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) LineWidthStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC2_LineWidthStrip",HM_Id,Structure_Id[file_id].c_str()));
   LineWidthStrip_iv[id]->GetXaxis()->SetTitle("LineWidthStrip Current [A]");
   LineWidthStrip_iv[id]->GetYaxis()->SetTitle("LineWidthStrip Voltage [V]");
   LineWidthStrip_iv[id]->SetDrawOption("AP");
   LineWidthStrip_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   LineWidthStrip_iv[id]->SetMarkerColor(colorid);
   //   LineWidthStrip_iv[id]->Draw();
   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   LineWidthStrip_cv_low[id] = LineWidthStrip_current_d[id][3];
   LineWidthStrip_cv_high[id] = LineWidthStrip_current_d[id][38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", LineWidthStrip_cv_low[id], LineWidthStrip_cv_high[id]);
//   LineWidthStrip_iv[id]->Fit("f_cv","0R+");
   LineWidthStrip_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   LineWidthStrip_A_constant[id] = LineWidthStrip_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   LineWidthStrip_A_slope[id]  = LineWidthStrip_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double R_sh = LineWidthStrip_A_slope[id]; // LineWidthStrip R_sh

   cout << "----->>>   R_sh Voltage = " << R_sh << " [Ohm/sq]" << endl;

   return R_sh;

}
void LineWidthStrip_PQC_Flute2_Analysis_Final(int nFiles)
{
   gStyle->SetOptStat(0);

   VDPStrip_PQC_Flute1_Analysis_Final(nFiles); // get the R_Geometry

   for(int i=0;i<nFiles;i++) {
    LineWidthStrip_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC2_LinewidthStrip.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << LineWidthStrip_DataFileName[2*i] << endl;

    LineWidthStrip_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC2_LinewidthStrip.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << LineWidthStrip_DataFileName[2*i+1] << endl;
       cout<<"caacacqwcqwc"<< endl;
   }


   float d = 128.5; // in um
   for (int i=0;i<nFiles;i++) {

     LineWidthStrip_Rsh_E[i] = LineWidthStrip_Read_Info_from_File(LineWidthStrip_DataFileName[2*i],(2*i),i);
     cout << "VDPStrip_s_Rsh_E[i]" << VDPStrip_s_Rsh_E[i] << endl;
     LineWidthStrip_s_Rsh_Rgeom_Corrected_E[i] = (1/LineWidthStrip_Rsh_E[i])*VDPStrip_s_Rsh_E[i]*d;
     LineWidthStrip_r_Rsh_Rgeom_Corrected_E[i] = (1/LineWidthStrip_Rsh_E[i])*VDPStrip_r_Rsh_E[i]*d;

     LineWidthStrip_Rsh_W[i] = LineWidthStrip_Read_Info_from_File(LineWidthStrip_DataFileName[2*i+1],(2*i+1),i);
     LineWidthStrip_s_Rsh_Rgeom_Corrected_W[i] = (1/LineWidthStrip_Rsh_W[i])*VDPStrip_s_Rsh_W[i]*d;
     LineWidthStrip_r_Rsh_Rgeom_Corrected_W[i] = (1/LineWidthStrip_Rsh_W[i])*VDPStrip_r_Rsh_W[i]*d;

     LineWidthStrip_all[2*i+0] = LineWidthStrip_s_Rsh_Rgeom_Corrected_E[i]; // put all std values in one array; for the xml production only std needed
     LineWidthStrip_all[2*i+1] = LineWidthStrip_s_Rsh_Rgeom_Corrected_W[i];
     LineWidthStrip_res_all[2*i+0] = LineWidthStrip_Rsh_E[i]; // put all std values in one array; for the xml production only std needed
     LineWidthStrip_res_all[2*i+1] = LineWidthStrip_Rsh_W[i];
   }


   TCanvas *LineWidthStrip_cc_final[4];
   TMultiGraph *LineWidthStrip_mg[4];
   for (int i=0;i<2;i++) {
     LineWidthStrip_cc_final[i] = new TCanvas(Form("LineWidthStrip_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     LineWidthStrip_mg[i] = new TMultiGraph(Form("LineWidthStrip_mg_%d",i),"LineWidthStrip");
     if(i == 0) LineWidthStrip_mg[i]->SetTitle("LineWidthStrip_E");
     if(i == 1) LineWidthStrip_mg[i]->SetTitle("LineWidthStrip_W");

    // LineWidthStrip_iv[4*i+]->Draw("ALP");
     for(int id=0;id<nFiles;id++) {

        LineWidthStrip_mg[i]->Add(LineWidthStrip_iv[2*id+i]);

     }
     LineWidthStrip_mg[i]->Draw("LPsame");
     LineWidthStrip_mg[i]->GetXaxis()->SetTitle("Current [A]");
     LineWidthStrip_mg[i]->GetYaxis()->SetTitle("Voltage [V]");
     LineWidthStrip_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     LineWidthStrip_cc_final[i]->BuildLegend(0.15,0.5,0.5,0.85);


     gPad->Modified();
     gPad->Update();
   }


   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;

   cout << "      LineWidthStrip   Id    |  LineWidthStrip_s_Rsh [um]   | LineWidthStrip_r_Rsh [um]|" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(18) <<  setprecision(4) << LineWidthStrip_s_Rsh_Rgeom_Corrected_E[id]  <<  setw(30)  << LineWidthStrip_r_Rsh_Rgeom_Corrected_E[id]
	       << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(18) <<  setprecision(4) << LineWidthStrip_s_Rsh_Rgeom_Corrected_W[id]  <<  setw(30)  << LineWidthStrip_r_Rsh_Rgeom_Corrected_W[id]
	       << endl;
   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

}
void LineWidthStrip_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  cout << "LineWidthStrip Xml production is on" << endl;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<=1; j++){
          //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = LineWidthStrip_DataFileName[2*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = LineWidthStrip_DataFileName[2*i+j].Last(delimeter2);
         xml_filename = LineWidthStrip_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[4*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);

        // First create engine
        TXMLEngine LineWidthStrip_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = LineWidthStrip_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = LineWidthStrip_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = LineWidthStrip_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = LineWidthStrip_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = LineWidthStrip_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = LineWidthStrip_xml.NewChild(HEADER, 0, "RUN");
              LineWidthStrip_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              LineWidthStrip_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              LineWidthStrip_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              LineWidthStrip_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              LineWidthStrip_xml.NewChild(RUN, 0, "INITIATED_BY_USER", LineWidthStrip_Operators[2*i+j].c_str());
              LineWidthStrip_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", LineWidthStrip_Begin_Timestamp[2*i+j].c_str());
              LineWidthStrip_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = LineWidthStrip_xml.NewChild(ROOT, 0, "DATA_SET");
            LineWidthStrip_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            LineWidthStrip_xml.NewChild(DATA_SET, 0, "VERSION", "IV_measurement");
            XMLNodePointer_t PART = LineWidthStrip_xml.NewChild(DATA_SET, 0, "PART");
              LineWidthStrip_xml.NewChild(PART, 0, "NAME_LABEL", LineWidthStrip_name_labels[2*i+j].c_str());
              LineWidthStrip_xml.NewChild(PART, 0, "KIND_OF_PART", LineWidthStrip_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = LineWidthStrip_xml.NewChild(DATA_SET, 0, "DATA");
              LineWidthStrip_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC2");
              LineWidthStrip_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", "LINEWIDTH_STRIP");
              LineWidthStrip_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", LineWidthStrip_set_id[2*i+j].c_str());
              LineWidthStrip_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", LineWidthStrip_config_id[2*i+j].c_str());
              LineWidthStrip_xml.NewChild(DATA, 0, "EQUIPMENT", LineWidthStrip_equipment[2*i+j].c_str());
              LineWidthStrip_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(LineWidthStrip_waiting_time[2*i+j], 3).c_str());
              LineWidthStrip_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(LineWidthStrip_temp[2*i+j], 3).c_str());
              LineWidthStrip_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(LineWidthStrip_av_temp[2*i+j], 3).c_str());
              LineWidthStrip_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", LineWidthStrip_struct_id[2*i+j].c_str());
              LineWidthStrip_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = LineWidthStrip_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = LineWidthStrip_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = LineWidthStrip_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  LineWidthStrip_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  LineWidthStrip_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = LineWidthStrip_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                LineWidthStrip_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                LineWidthStrip_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "IV_measurement");
                XMLNodePointer_t PART_CDS1 = LineWidthStrip_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  LineWidthStrip_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", LineWidthStrip_name_labels[2*i+j].c_str());
                  LineWidthStrip_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", LineWidthStrip_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<LineWidthStrip_Nmeas[2*i+j]-1; k++){
                    XMLNodePointer_t DATA_CDS1 = LineWidthStrip_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                      LineWidthStrip_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string_scientific_format(LineWidthStrip_voltage[2*i+j][k], 3).c_str());
                      LineWidthStrip_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(LineWidthStrip_current[2*i+j][k]*1E+9, 3).c_str());
                      LineWidthStrip_xml.NewChild(DATA_CDS1, 0, "TIME", LineWidthStrip_timestamp_meas[2*i+j][k].c_str());
                      LineWidthStrip_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(LineWidthStrip_temp_meas[2*i+j][k], 3).c_str());
                      LineWidthStrip_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(LineWidthStrip_air_temp_meas[2*i+j][k], 3).c_str());
                      LineWidthStrip_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(LineWidthStrip_rh_prcnt_meas[2*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = LineWidthStrip_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = LineWidthStrip_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = LineWidthStrip_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      LineWidthStrip_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      LineWidthStrip_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = LineWidthStrip_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    LineWidthStrip_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    LineWidthStrip_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement");
                    XMLNodePointer_t PART_CDS2 = LineWidthStrip_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      LineWidthStrip_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", LineWidthStrip_name_labels[2*i+j].c_str());
                      LineWidthStrip_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", LineWidthStrip_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = LineWidthStrip_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      LineWidthStrip_xml.NewChild(DATA_CDS2, 0, "T_UM", Conv_float_to_string(LineWidthStrip_all[2*i+j], 3).c_str());
                      LineWidthStrip_xml.NewChild(DATA_CDS2, 0, "R_OHM", Conv_float_to_string_scientific_format(LineWidthStrip_res_all[2*i+j], 3).c_str());

        XMLDocPointer_t LineWidthStrip_xmldoc = LineWidthStrip_xml.NewDoc();
        LineWidthStrip_xml.DocSetRootElement(LineWidthStrip_xmldoc, ROOT);
        // Save document to file
        LineWidthStrip_xml.SaveDoc(LineWidthStrip_xmldoc, Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        LineWidthStrip_xml.FreeDoc(LineWidthStrip_xmldoc);
      }
    }
  }
}
// ==================================================================================================//
double RPoly_Read_Info_from_File(TString RPoly_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(RPoly_DataFileName);
   char delimeter1 = '/';
   int pos1 = RPoly_DataFileName.Last(delimeter1);
   char delimeter2 = '.';
   int pos2 = RPoly_DataFileName.Last(delimeter2);
   TString RPoly_filename = RPoly_DataFileName(pos1+1, ((pos2-pos1)-1));
   TObjArray *tx = RPoly_filename .Tokenize("_");
   //tx->Print();
   vector<TString> name_segments;
   for (Int_t i = 0; i < tx->GetEntries(); i++){
     name_segments.push_back(((TObjString *)(tx->At(i)))->String());
     //std::cout << ((TObjString *)(tx->At(i)))->String() << std::endl;
   }
   //cout << A[5] << endl;

   cout << endl << "=============>>>>>> Analyze File : " << RPoly_DataFileName << " wit Id = " << id << endl;
   int nlines = 0;
   int N_meas = 0;
   string line;
   while (in.good()) {
     getline (in,line);
      if(nlines < 18) {
         if(nlines == 1) {
         //cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         RPoly_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         RPoly_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         RPoly_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         RPoly_kind_of_parts.push_back(Kind_of_part);
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
         RPoly_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         //RPoly_struct_id.push_back(name_segments[4]);
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         RPoly_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         //RPoly_set_id.push_back(name_segments[5]);
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         RPoly_config_id.push_back(Kind_of_HM_config_id);
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
         RPoly_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         RPoly_waiting_time.push_back(Waiting_time_s);
         }
         if(nlines == 12) {
         cout << line << '\n';
         double Bias_V;
         sscanf(line.c_str(), "%*s %lf", &Bias_V);
         cout<<Bias_V<<endl;
         RPoly_bias.push_back(Bias_V);
         }
         if(nlines == 13) {
         cout << line << '\n';
         double Temp_set_degC;
         sscanf(line.c_str(), "%*s %lf", &Temp_set_degC);
         cout<<Temp_set_degC<<endl;
         RPoly_temp.push_back(Temp_set_degC);
        }
        if(nlines == 14) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        RPoly_av_temp.push_back(Av_temp_degC);
       }
        if(nlines == 15) {
        cout << line << '\n';
        double ComplianceA;
        sscanf(line.c_str(), "%*s %lf", &ComplianceA);
        cout<<ComplianceA<<endl;
        RPoly_compliance.push_back(ComplianceA);
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

    	RPoly_voltage_d[id][N_meas] = voltage;
    	RPoly_current_d[id][N_meas] = current_namp;
      RPoly_temp_meas.push_back(std::vector<double>());
      RPoly_temp_meas[id].push_back(temp_degC);
      RPoly_air_temp_meas.push_back(std::vector<double>());
      RPoly_air_temp_meas[id].push_back(air_temp_degC);
      RPoly_rh_prcnt_meas.push_back(std::vector<double>());
      RPoly_rh_prcnt_meas[id].push_back(rh_prcnt);
      RPoly_timestamp_meas.push_back(std::vector<string>());
      RPoly_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
    	RPoly_voltage[id][N_meas] = (float) RPoly_voltage_d[id][N_meas];
      RPoly_current[id][N_meas] = (float) RPoly_current_d[id][N_meas];
      //cout<< RPoly_vv[N_meas] << "," << RPoly_Curr[N_meas] << endl;
      N_meas++;
      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   RPoly_Nmeas.push_back(N_meas);

   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;

   TGaxis::SetMaxDigits(3);

   //   TCanvas *cc_spline = new TCanvas();
   RPoly_iv[id] =  new TGraph(N_meas, RPoly_current_d[id], RPoly_voltage_d[id]);

   if((id % 2) == 0) RPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC2_RPoly",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) RPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC2_RPoly",HM_Id,Structure_Id[file_id].c_str()));

   RPoly_iv[id]->GetXaxis()->SetTitle("RPoly Current [A]");
   RPoly_iv[id]->GetYaxis()->SetTitle("RPoly Voltage [V]");
   RPoly_iv[id]->SetDrawOption("AP");
   RPoly_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   RPoly_iv[id]->SetMarkerColor(colorid);
   //   RPoly_iv[id]->Draw();
   // Set linear fit ranges
   // -------------------------
   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   RPoly_cv_low[id] = RPoly_current_d[id][0];
   RPoly_cv_high[id] = RPoly_current_d[id][N_meas-1];
   // ---------------------------------------------------------------------------------------------------------------
   TF1 *f_cv = new TF1("f_cv","pol1", RPoly_cv_low[id], RPoly_cv_high[id]);
   //   RPoly_iv[id]->Fit("f_cv","0R+");
   RPoly_iv[id]->Fit("f_cv","0R+");
   // Accumulation Region
   // ---------------------
   RPoly_A_constant[id] = RPoly_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   RPoly_A_slope[id]  = RPoly_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double R_sh = RPoly_A_slope[id]*1000; // RPoly R_sh

   cout << "----->>>   R_sh Voltage = " << R_sh  << " [MOhm]" << endl;

   return R_sh;

}
void RPoly_IV_PQC2_Analysis_Final(int nFiles)
{
   //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1
   //   gROOT->Reset();
   gStyle->SetOptStat(0);

   for(int i=0;i<nFiles;i++) {
    RPoly_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Left_PQC2_RPoly.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    RPoly_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Left_PQC2_RPoly.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }
    //   float s = 187e-6; // in m
    //   s *= 100; // in cm
    //   float F1 = 1.089;
    //   float F2 = 1.218;
   for (int i=0;i<nFiles;i++) {
     cout << "RPoly Structure : Analyze file : " << RPoly_DataFileName[2*i] << endl;
     RPoly_E[i] = RPoly_Read_Info_from_File(RPoly_DataFileName[2*i],(2*i),i);
     //     RPoly_Resistivity_s_E[i] = 2*TMath::Pi()*s*F1*RPoly_E[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());
     cout << "RPoly Structure :Analyze file : " << RPoly_DataFileName[2*i+1] << endl;
     RPoly_W[i] = RPoly_Read_Info_from_File(RPoly_DataFileName[2*i+1],(2*i+1),i);
     //     RPoly_Resistivity_s_W[i] = 2*TMath::Pi()*s*F1*RPoly_W[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());
     RPoly_all[2*i+0] = RPoly_E[i];
     RPoly_all[2*i+1] = RPoly_W[i];
   }
   TCanvas *RPoly_cc_final[2];
   TMultiGraph *RPoly_mg[2];
   for (int i=0;i<2;i++) {
     RPoly_cc_final[i] = new TCanvas(Form("RPoly_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     RPoly_mg[i] = new TMultiGraph(Form("RPoly_mg_%d",i),"RPoly R_sh");
     if(i == 0) RPoly_mg[i]->SetTitle("RPoly_E");
     if(i == 1) RPoly_mg[i]->SetTitle("RPoly_W");
    // RPoly_iv[4*i+]->Draw("ALP");
     for(int id=0;id<nFiles;id++) {
        RPoly_mg[i]->Add(RPoly_iv[2*id+i]);
     }
     RPoly_mg[i]->Draw("LPsame");
     RPoly_mg[i]->GetXaxis()->SetTitle("Current [A]");
     RPoly_mg[i]->GetYaxis()->SetTitle("Voltage [V]");
     RPoly_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();
     RPoly_cc_final[i]->BuildLegend(0.15,0.5,0.5,0.85);
     gPad->Modified();
     gPad->Update();
   }
   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "      RPoly   Id                 |  RPoly [MOhm] " <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(20) <<  setprecision(4) << RPoly_E[id]
	       << endl;
         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(20) <<  setprecision(4) << RPoly_W[id]
	       << endl;

   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;


}
//--------------------------------------------------------------------------------------------------//
// RPoly xml production
// --------------------------------------------------------------------------------------------------//
void RPoly_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    cout << "RPoly Xml production is on" << endl;
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<=1; j++){
        //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = RPoly_DataFileName[2*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = RPoly_DataFileName[2*i+j].Last(delimeter2);
         xml_filename = RPoly_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[2*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);

         // First create engine
        TXMLEngine RPoly_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = RPoly_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = RPoly_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = RPoly_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = RPoly_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = RPoly_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = RPoly_xml.NewChild(HEADER, 0, "RUN");
              RPoly_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              RPoly_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              RPoly_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              RPoly_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              RPoly_xml.NewChild(RUN, 0, "INITIATED_BY_USER", RPoly_Operators[2*i+j].c_str());
              RPoly_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", RPoly_Begin_Timestamp[2*i+j].c_str());
              RPoly_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = RPoly_xml.NewChild(ROOT, 0, "DATA_SET");
            RPoly_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            RPoly_xml.NewChild(DATA_SET, 0, "VERSION", "v2");
            XMLNodePointer_t PART = RPoly_xml.NewChild(DATA_SET, 0, "PART");
              RPoly_xml.NewChild(PART, 0, "NAME_LABEL", RPoly_name_labels[2*i+j].c_str());
              RPoly_xml.NewChild(PART, 0, "KIND_OF_PART", RPoly_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = RPoly_xml.NewChild(DATA_SET, 0, "DATA");
              RPoly_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC2");
              RPoly_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", "R_POLY");
              RPoly_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", RPoly_set_id[2*i+j]);
              RPoly_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", RPoly_config_id[2*i+j].c_str());
              RPoly_xml.NewChild(DATA, 0, "EQUIPMENT", RPoly_equipment[2*i+j].c_str());
              RPoly_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(RPoly_waiting_time[2*i+j],3).c_str());
              RPoly_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(RPoly_temp[2*i+j],3).c_str());
              RPoly_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(RPoly_av_temp[2*i+j],3).c_str());
              RPoly_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "RPoly1");
              RPoly_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = RPoly_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = RPoly_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = RPoly_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  RPoly_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  RPoly_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = RPoly_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                RPoly_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                RPoly_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "V2");
                XMLNodePointer_t PART_CDS1 = RPoly_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  RPoly_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", RPoly_name_labels[2*i+j].c_str());
                  RPoly_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", RPoly_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<RPoly_Nmeas[2*i+j]-1; k++){
                XMLNodePointer_t DATA_CDS1 = RPoly_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  RPoly_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(RPoly_voltage[2*i+j][k],3).c_str());
                  //RPoly_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", to_string(RPoly_current[k]).c_str());
                  RPoly_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string(RPoly_current[2*i+j][k],3).c_str());
                  RPoly_xml.NewChild(DATA_CDS1, 0, "TIME", RPoly_timestamp_meas[2*i+j][k].c_str());
                  RPoly_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(RPoly_temp_meas[2*i+j][k], 3).c_str());
                  RPoly_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(RPoly_air_temp_meas[2*i+j][k], 3).c_str());
                  RPoly_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(RPoly_rh_prcnt_meas[2*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = RPoly_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = RPoly_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = RPoly_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      RPoly_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      RPoly_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = RPoly_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    RPoly_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    RPoly_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement-004");
                    XMLNodePointer_t PART_CDS2 = RPoly_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      RPoly_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", RPoly_name_labels[2*i+j].c_str());
                      RPoly_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", RPoly_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = RPoly_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      RPoly_xml.NewChild(DATA_CDS2, 0, "RPOLY_MOHM", Conv_float_to_string(RPoly_all[2*i+j], 3).c_str());
        XMLDocPointer_t RPoly_xmldoc = RPoly_xml.NewDoc();
        RPoly_xml.DocSetRootElement(RPoly_xmldoc, ROOT);
        // Save document to file
        RPoly_xml.SaveDoc(RPoly_xmldoc, Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        RPoly_xml.FreeDoc(RPoly_xmldoc);
      }
    }
  }
}


void GCD_Read_IV_from_File(TString GCD_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(GCD_DataFileName);

   epsilon_0 = 0.0885418782; // [pF/cm]
   epsilon_Si = 11.68*epsilon_0; // [pF/cm]
   epsilon_SiO2 = 3.9*epsilon_0; // [pF/cm]
   GCD_area = 0.505e-2; // [cm^2]
   //GCD_area = 0.723e-2; // [cm^2]
   q = 1.60219E-7; // [pF.V]
   kB = 8.6173324E-5; // [eV/K]
   //   TKelvin = 273 + Temp_measurement[id]; // [K]
   //   Nintrisic = 1.45e10; // [cm^-3]
   cout << endl << "=============>>>>>> Analyze File : " << GCD_DataFileName << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;
   while (in.good()) {
     getline (in,line);
      if(nlines < 18) {
         if(nlines == 1) {
         //cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         GCD_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         GCD_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         GCD_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         GCD_kind_of_parts.push_back(Kind_of_part);
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
         GCD_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         //GCD_struct_id.push_back(name_segments[4]);
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         GCD_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         //GCD_set_id.push_back(name_segments[5]);
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         GCD_config_id.push_back(Kind_of_HM_config_id);
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
         GCD_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         GCD_waiting_time.push_back(Waiting_time_s);
         }
         if(nlines == 12) {
         cout << line << '\n';
         double Bias_V;
         sscanf(line.c_str(), "%*s %lf", &Bias_V);
         cout<<Bias_V<<endl;
         GCD_bias.push_back(Bias_V);
         }
         if(nlines == 13) {
         cout << line << '\n';
         double Temp_set_degC;
         sscanf(line.c_str(), "%*s %lf", &Temp_set_degC);
         cout<<Temp_set_degC<<endl;
         TKelvin = 273 + Temp_set_degC; // [K]
         GCD_temp.push_back(Temp_set_degC);
        }
        if(nlines == 14) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        GCD_av_temp.push_back(Av_temp_degC);
       }
        if(nlines == 15) {
        cout << line << '\n';
        double ComplianceA;
        sscanf(line.c_str(), "%*s %lf", &ComplianceA);
        cout<<ComplianceA<<endl;
        GCD_compliance.push_back(ComplianceA);
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
    	GCD_vv_d[id][N_meas] = voltage;
    	GCD_current_d[id][N_meas] = current_namp*1E+3;
      GCD_temp_meas.push_back(std::vector<double>());
      GCD_temp_meas[id].push_back(temp_degC);
      GCD_air_temp_meas.push_back(std::vector<double>());
      GCD_air_temp_meas[id].push_back(air_temp_degC);
      GCD_rh_prcnt_meas.push_back(std::vector<double>());
      GCD_rh_prcnt_meas[id].push_back(rh_prcnt);
      GCD_timestamp_meas.push_back(std::vector<string>());
      GCD_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
    	GCD_vv[id][N_meas] = (float) GCD_vv[id][N_meas];
      GCD_current[id][N_meas] = (float) GCD_current_d[id][N_meas];
      //cout<< GCD_vv[N_meas] << "," << GCD_Curr[N_meas] << endl;
      N_meas++;
      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   GCD_Nmeas.push_back(N_meas);
   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl << endl;

   //TKelvin = 273 + GCD_Temp_measurement[id]; // [K]
   Nintrisic = 5.29E+19*pow((TKelvin/300),2.54)*exp(-6726/TKelvin);

   cout << endl << ">>>>>>>>>>> Nintrisic  = " << Nintrisic << " cm^{-3} and area = " << GCD_area <<  " cm^{-2} " << endl;
   int  N_meas_reduced = N_meas-1;

   if(N_meas_reduced > 300 ) {
      cout << "Serious Problem : Measurements > 300 that TSpline3 can handle...  " << endl;
      cout << " Solution: Reduce number of measurements for TSpline3 rootine.... " << endl;
      cout << " Now Abort..." << endl;
      abort();
   }
   double GCD_vv_d_reduced[300];
   double GCD_current_d_reduced[300];
   for (Int_t i=0; i<N_meas_reduced; i++) {
      GCD_vv_d_reduced[i] = GCD_vv_d[id][i];
      GCD_current_d_reduced[i] = GCD_current_d[id][i];
      //      cout << i << ")  "  <<  GCD_vv_d[i] << "   " << GCD_current_d[i] << endl;
   }
    //   CSpline_Interpolator[id] = new ROOT::Math::Interpolator();
    //   CSpline_Interpolator[id]->SetData(N_meas_reduced,vv_d_reduced,current_d_inversed);
    GCD_CuSpl_iv[id] = new TSpline3("Cubic Spline", GCD_vv_d_reduced, GCD_current_d_reduced, N_meas_reduced, "cb2e2", 0, 0);  // Seems that works well for less than 300 points
    //GCD_cc_spline[id] = new TCanvas(Form("GCD_cc_spline_%d",id),"", 700, 500, 700,700);
    GCD_iv[id] =  new TGraph(N_meas-2,GCD_vv_d_reduced,GCD_current_d_reduced);
    if((id % 2) == 0) GCD_iv[id]->SetTitle(Form("GCD%d_0%s E IV",HM_Id,Structure_Id[file_id].c_str()));
    if((id % 2) == 1) GCD_iv[id]->SetTitle(Form("GCD%d_0%s W IV",HM_Id,Structure_Id[file_id].c_str()));
    GCD_iv[id]->GetXaxis()->SetTitle("Gate Voltage [V]");
    GCD_iv[id]->GetYaxis()->SetTitle("GCD Current [pA]");
    GCD_iv[id]->GetYaxis()->SetRangeUser(-1,15);
    GCD_iv[id]->SetDrawOption("AP");
    GCD_iv[id]->SetMarkerStyle(20+file_id);
    int colorid = 1+file_id;
    if(colorid == 5) colorid = 46;
    GCD_iv[id]->SetMarkerColor(colorid);
    //GCD_iv[id]->Draw();
    GCD_CuSpl_iv[id]->SetLineColor(kBlue);
    //GCD_CuSpl_iv[id]->Draw("same");


    // Get V_Fb grosso - modo

    // --------------------------

    for(int i = 0;i<N_meas_reduced; i++) {
       GCD_dev_spline_simple[i] = GCD_CuSpl_iv[id]->Derivative(GCD_vv_d_reduced[i]);
    }

    //GCD_cc_spline_der[id] = new TCanvas(Form("GCD_cc_spline_der_%d",id),"", 50, 200, 500,500);
    TGraph *GCD_iv_deriv =  new TGraph(N_meas_reduced,GCD_vv_d_reduced,GCD_dev_spline_simple);
    if((id % 2) == 0) GCD_iv_deriv->SetTitle(Form("GCD%d_0%s E IV derivarive",HM_Id,Structure_Id[file_id].c_str()));
    if((id % 2) == 1) GCD_iv_deriv->SetTitle(Form("GCD%d_0%s W IV derivarive",HM_Id,Structure_Id[file_id].c_str()));
    //GCD_iv_deriv->Draw();

    /*
        for (Int_t i=0; i<N_meas_reduced; i++) {

          cout << "Data points : " << i << ")  "  <<  GCD_vv_d_reduced[i] << "   " << "  "
               << GCD_current_d_reduced[i] <<  "  " << GCD_CuSpl_iv[id]->Eval(GCD_vv_d_reduced[i]) << "   " << GCD_dev_spline_simple[i] << endl;

        }
    */
    // get the minimum of the derivative
    int maxdev_index_simple = 0;
    float max_GCD_dev_spline_simple = 0;
    for (int i = 1;i<N_meas_reduced; i++) {
       if(GCD_dev_spline_simple[i] >  max_GCD_dev_spline_simple && GCD_vv_d_reduced[i] > -9 &&  GCD_vv_d_reduced[i] < -7) {
         max_GCD_dev_spline_simple = GCD_dev_spline_simple[i];
	 maxdev_index_simple = i;
       }
    }
    cout  << "---------------- >>>>> Max derivative Simple at i = " << maxdev_index_simple
          << " with value = " << GCD_dev_spline_simple[maxdev_index_simple] << endl;

    // Get V_Fb fine tuned
    // -------------------------
    int Nub_steps = 1000;

    float vv_initial = GCD_vv_d_reduced[0];
    float vv_final = GCD_vv_d_reduced[N_meas-1];

    Double_t dev_spline[2000];
    Float_t dev_spline_f[2000];
    Float_t voltage_extented[2000];
    for(int i = 0;i<Nub_steps; i++) {
      voltage_extented[i] =  vv_initial + i*(vv_final - vv_initial)/(Nub_steps-1);
      dev_spline[i] = GCD_CuSpl_iv[id]->Derivative(voltage_extented[i]);
      dev_spline_f[i] = float(dev_spline[i]);
    }

    // get the minimum of the derivative
    int min_index = 0;
    float min_dev_spline_f = 0.0;
    int max_index = 0;
    float max_dev_spline_f = 0.0;
    for (int i = 1;i<Nub_steps; i++) {
       if(dev_spline_f[i] < min_dev_spline_f && voltage_extented[i] > -4 &&  voltage_extented[i] < 0) {
         min_dev_spline_f = dev_spline_f[i];
	        min_index = i;
       }
       if(dev_spline_f[i] > max_dev_spline_f && voltage_extented[i] > -9 &&  voltage_extented[i] < -7) {
         max_dev_spline_f = dev_spline_f[i];
	        max_index = i;
       }
    }
    double V_FB_derivative = double(voltage_extented[max_index]) - GCD_bias[id];
    double I_FB_derivative =GCD_CuSpl_iv[id]->Eval(V_FB_derivative);
    cout  << "---------------- >>>>> Max derivative accurate at i = " << max_index
          << "  with value = " << dev_spline_f[max_index]
	  << " V_flatband = " <<  V_FB_derivative
	  << " and I_flatband = " << I_FB_derivative << endl;

    double V_FB_plus_V_Bias_derivative = double(voltage_extented[min_index]);
    double I_FB_plus_V_Bias_derivative =GCD_CuSpl_iv[id]->Eval(V_FB_plus_V_Bias_derivative);
    cout  << "---------------- >>>>> Min derivative accurate at i = " << min_index
          << "  with value = " << dev_spline_f[min_index]
	  << " V_flatband = " <<  V_FB_plus_V_Bias_derivative
	  << " and I_flatband_plus_B_Bias = " << I_FB_plus_V_Bias_derivative << endl << endl;

   GCD_VFB_Derivative[id] = V_FB_derivative; // in V
   cout << endl<<  "======>>> V_flatband = " <<  V_FB_derivative << endl;


   // find the maximum of the IV curve in range maxdev_index_simple + 16
   int IV_max_index = 0;
   double IV_Max = 0.0;
   for (int i = maxdev_index_simple;i<maxdev_index_simple+16; i++) {
       if(GCD_current_d[id][i] > IV_Max ) {
        IV_Max  = GCD_current_d[id][i];
	IV_max_index  = i;
       }
    }
    cout  << "---------------- >>>>> Max IV at i = " <<  IV_max_index
          << " with value = " << IV_Max << endl;


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------

   GCD_iv_low[id] = GCD_vv_d_reduced[IV_max_index-20];
   GCD_iv_high[id] = GCD_vv_d_reduced[IV_max_index-30];

   // Depletion Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   int Dep_Central_Point = IV_max_index;

   GCD_iv1_low[id] = GCD_vv_d_reduced[Dep_Central_Point];
   GCD_iv1_high[id] = GCD_vv_d_reduced[Dep_Central_Point+5];


   // Inversion Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   GCD_iv2_low[id] = GCD_vv_d_reduced[N_meas-60];;
   GCD_iv2_high[id] = GCD_vv_d_reduced[N_meas-45];


// Get tree branches to find the fit points
// -------------------------------------------

//   Int_t nentries = (Int_t)T->GetEntries();
//   cout << "Number of Measurements = " << nentries << endl;

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol0", GCD_iv_low[id], GCD_iv_high[id]);
//   f_cv->SetParameter(1,0.0);

   GCD_iv[id]->Fit("f_cv","R+");

   TF1 *f_cv1 = new TF1("f_cv1","pol0",GCD_iv1_low[id], GCD_iv1_high[id]);

   GCD_iv[id]->Fit("f_cv1","R+");

   TF1 *f_cv2 = new TF1("f_cv2","pol0", GCD_iv2_low[id], GCD_iv2_high[id]);

   GCD_iv[id]->Fit("f_cv2","R+");

   // Accumulation Region
   // ---------------------
   GCD_A_constant[id] = GCD_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   GCD_A_slope[id]  = 0.0; // Accumulation line slope

   GCD_I_accumulation[id] = GCD_A_constant[id]; // GCD Accumulation Current
   float I_accumulation_Error = GCD_iv[id]->GetFunction("f_cv")->GetParError(0);

   cout << endl << "======>>>   Accumulation Current = " <<  GCD_I_accumulation[id] << " +/- " << I_accumulation_Error << " [pA]" << endl;

   // Depletion Region
   // ---------------------
   GCD_D_constant[id] = GCD_iv[id]->GetFunction("f_cv1")->GetParameter(0); // Depletion Region constant term
   GCD_D_slope[id] = 0.0; // Depletion Region slope

   GCD_I_depletion[id] = GCD_D_constant[id]; // GCD Depletion Current
   float GCD_I_depletion_Error = GCD_iv[id]->GetFunction("f_cv1")->GetParError(0);

   cout << "======>>>   Depletion Current = " << GCD_I_depletion[id] << " +/- " << GCD_I_depletion_Error << " [pA]" << endl;

   // Inversion Region
   // ---------------------
   GCD_I_constant[id] = GCD_iv[id]->GetFunction("f_cv2")->GetParameter(0); // Inversion constant term
   GCD_I_slope[id]  = 0.0; // Inversion line slope

   GCD_I_inversion[id] = GCD_I_constant[id]; // GCD High Frequency Current - Inversion Current
   float GCD_I_inversion_Error = GCD_iv[id]->GetFunction("f_cv2")->GetParError(0);

   cout << "======>>>   Inversion Current = " << GCD_I_inversion[id] << " +/- " << GCD_I_inversion_Error << " [pA]" << endl;

   GCD_I_surface[id] = GCD_I_depletion[id] - GCD_I_inversion[id];

   cout << "======>>>   I_surface = " << GCD_I_surface[id] << " +/- "
        << sqrt(GCD_I_depletion_Error*GCD_I_depletion_Error + GCD_I_inversion_Error*GCD_I_inversion_Error) << " [pA]" << endl;

   GCD_S_interface_recombination_Velocity[id] = GCD_I_surface[id]/(q*Nintrisic*GCD_area);

   cout << "======>>>   Interface Recombination Velocity (S_{0}) = " << GCD_S_interface_recombination_Velocity[id] << " [cm/sec]" <<endl;
   cout << "======>>>   V_flatband = " <<  V_FB_derivative << endl;

   float sigma = 2.5e-16; // [cm^2]
   float u_thermal = 2.0e7; // [cm/sec]
   GCD_D_it[id] = GCD_S_interface_recombination_Velocity[id]/(3.14159 * sigma * u_thermal *  TKelvin * kB);
   GCD_N_it[id] = 1.12*GCD_D_it[id]/2;;


}
void GCD_PQC_Analysis_Final(int nFiles)
{

   gStyle->SetOptStat(0);


   //TString GCD_DataFileName[100];
   for(int i=0;i<nFiles;i++) {
      GCD_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Left_PQC2_GCD_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
      GCD_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Left_PQC2_GCD_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }


   for (int i=0;i<nFiles;i++) {

     GCD_Read_IV_from_File(GCD_DataFileName[2*i],2*i,i);
     GCD_Read_IV_from_File(GCD_DataFileName[2*i+1],2*i+1,i);
   }

   TCanvas *GCD_cc_final[2];
   TMultiGraph *GCD_mg[2];
   for (int i=0;i<2;i++) {
     GCD_cc_final[i] = new TCanvas(Form("cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     GCD_mg[i] = new TMultiGraph(Form("GCD_mg_%d",i),"GCD IV");
     if(i == 0) GCD_mg[i]->SetTitle("GCD_E IV ");
     if(i == 1) GCD_mg[i]->SetTitle("GCD_W IV ");

     for(int id=0;id<nFiles;id++) {

        GCD_mg[i]->Add(GCD_iv[2*id+i]);

     }
     GCD_mg[i]->Draw("LPsame");
     GCD_mg[i]->GetXaxis()->SetTitle("Gate GCD_voltage [V]");
     GCD_mg[i]->GetYaxis()->SetTitle("GCD_current [pF]");
     GCD_mg[i]->GetYaxis()->SetRangeUser(-1,15);
     GCD_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     GCD_cc_final[i]->BuildLegend(0.6,0.6,0.85,0.85);
     gPad->Modified();
     gPad->Update();
   }

   cout <<endl;
   cout << "=====================================================================================" << endl;
   cout << "------------ VPX28442: HMSet_Left GCD Rectangular Irradiated  -------------------------" << endl;
   cout << "=====================================================================================" << endl;

   cout << "         GCD Id                 | I_acc [pA]   | I_dep [pA]   |  I_inv [pA]  | I_Sur [pA]   | S_0 [cm/sec]    | V_FB [V]  |     N_it " << endl;
   cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

   for(int id=0;id<nFiles;id++) {

         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       <<  setw(15) <<  setprecision(4) << GCD_I_depletion[2*id]
	       <<  setw(15) <<  setprecision(4) << GCD_I_depletion[2*id]
	       <<  setw(15) <<  setprecision(2) << GCD_I_inversion[2*id]
	       <<  setw(15) <<  setprecision(3)<< GCD_I_surface[2*id]
	       <<  setw(15) << GCD_S_interface_recombination_Velocity[2*id]
	       <<  setw(15) << GCD_VFB_Derivative[2*id]
	       <<  setw(15) << GCD_N_it[2*id]    << endl;


         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       <<  setw(15) <<  setprecision(4) << GCD_I_depletion[2*id+1]
	       <<  setw(15) <<  setprecision(4) << GCD_I_depletion[2*id+1]
	       <<  setw(15) <<  setprecision(2) << GCD_I_inversion[2*id+1]
	       <<  setw(15) <<  setprecision(3)<< GCD_I_surface[2*id+1]
	       <<  setw(15) << GCD_S_interface_recombination_Velocity[2*id+1]
	       <<  setw(15) << GCD_VFB_Derivative[2*id+1]
	       <<  setw(15) << GCD_N_it[2*id+1]    << endl;


   }
   cout << "-------------------------------------------------------------------------------------" << endl;

}
void GCD_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    cout << "GCD Xml production is on" << endl;
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<=1; j++){
        //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = GCD_DataFileName[2*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = GCD_DataFileName[2*i+j].Last(delimeter2);
         xml_filename = GCD_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[2*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);

         // First create engine
        TXMLEngine GCD_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = GCD_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = GCD_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = GCD_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = GCD_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = GCD_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = GCD_xml.NewChild(HEADER, 0, "RUN");
              GCD_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              GCD_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              GCD_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              GCD_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              GCD_xml.NewChild(RUN, 0, "INITIATED_BY_USER", GCD_Operators[2*i+j].c_str());
              GCD_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", GCD_Begin_Timestamp[2*i+j].c_str());
              GCD_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = GCD_xml.NewChild(ROOT, 0, "DATA_SET");
            GCD_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            GCD_xml.NewChild(DATA_SET, 0, "VERSION", "v2");
            XMLNodePointer_t PART = GCD_xml.NewChild(DATA_SET, 0, "PART");
              GCD_xml.NewChild(PART, 0, "NAME_LABEL", GCD_name_labels[2*i+j].c_str());
              GCD_xml.NewChild(PART, 0, "KIND_OF_PART", GCD_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = GCD_xml.NewChild(DATA_SET, 0, "DATA");
              GCD_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC2");
              GCD_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", "GCD");
              GCD_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", GCD_set_id[2*i+j]);
              GCD_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", GCD_config_id[2*i+j].c_str());
              GCD_xml.NewChild(DATA, 0, "EQUIPMENT", GCD_equipment[2*i+j].c_str());
              GCD_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(GCD_waiting_time[2*i+j],3).c_str());
              GCD_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(GCD_temp[2*i+j],3).c_str());
              GCD_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(GCD_av_temp[2*i+j],3).c_str());
              GCD_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "GCD1");
              GCD_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = GCD_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = GCD_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = GCD_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  GCD_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  GCD_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = GCD_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                GCD_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                GCD_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "V2");
                XMLNodePointer_t PART_CDS1 = GCD_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  GCD_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", GCD_name_labels[2*i+j].c_str());
                  GCD_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", GCD_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<GCD_Nmeas[2*i+j]; k++){
                XMLNodePointer_t DATA_CDS1 = GCD_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  GCD_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(GCD_vv_d[2*i+j][k],3).c_str());
                  //GCD_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", to_string(GCD_current[k]).c_str());
                  GCD_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(GCD_current[2*i+j][k]*1E-3,3).c_str());
                  GCD_xml.NewChild(DATA_CDS1, 0, "TIME", GCD_timestamp_meas[2*i+j][k].c_str());
                  GCD_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(GCD_temp_meas[2*i+j][k], 3).c_str());
                  GCD_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(GCD_air_temp_meas[2*i+j][k], 3).c_str());
                  GCD_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(GCD_rh_prcnt_meas[2*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = GCD_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = GCD_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = GCD_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      GCD_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      GCD_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = GCD_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    GCD_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    GCD_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement-004");
                    XMLNodePointer_t PART_CDS2 = GCD_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      GCD_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", GCD_name_labels[2*i+j].c_str());
                      GCD_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", GCD_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = GCD_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      GCD_xml.NewChild(DATA_CDS2, 0, "ISURF_PAMPR", Conv_float_to_string(GCD_I_surface[2*i+j], 3).c_str());
                      GCD_xml.NewChild(DATA_CDS2, 0, "S0_CMSEC", Conv_float_to_string(GCD_S_interface_recombination_Velocity[2*i+j], 3).c_str());
                      GCD_xml.NewChild(DATA_CDS2, 0, "VFB_ACC_V", Conv_float_to_string(GCD_VFB_Derivative[2*i+j], 3).c_str());
        XMLDocPointer_t GCD_xmldoc = GCD_xml.NewDoc();
        GCD_xml.DocSetRootElement(GCD_xmldoc, ROOT);
        // Save document to file
        GCD_xml.SaveDoc(GCD_xmldoc, Form("./XML_Info_Flute2/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        GCD_xml.FreeDoc(GCD_xmldoc);
      }
    }
  }
}

void Flute2_2_S_Quick_Characterization_with_xml(int nFiles)
{
   Diel_IV_PQC2_Analysis_Final(nFiles);
   Diel_IV_PQC_xml_production(xml_on, nFiles);

   LineWidthStop_PQC_Flute2_Analysis_Final(nFiles);
   LineWidthStop_PQC_xml_production(xml_on, nFiles);

   LineWidthStrip_PQC_Flute2_Analysis_Final(nFiles);
   LineWidthStrip_PQC_xml_production(xml_on, nFiles);

   RPoly_IV_PQC2_Analysis_Final(nFiles);
   RPoly_PQC_xml_production(xml_on, nFiles);

   GCD_PQC_Analysis_Final(nFiles);
   GCD_PQC_xml_production(xml_on, nFiles);

   CVS_Output_File.open(Form("CSV_Info/VPX%d/Flute2_VPX%d_0xx_2-S_HM_E_W_Left_Quick_Info.csv",HM_Id,HM_Id));

   cout << "=========================================================================================================================================" << endl;
   cout << "-----------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  ------------------------------------------------------" << endl;
   cout << "=========================================================================================================================================" << endl;
   cout << "     Flute2  Id            | DielNE | DielNW | DielSW | LineW_Pstop | LineW_Strip  | RPoly   | GCD I_surf | GCD S_0  | GCD Vfb_Acc_Dep |" <<endl;
   cout << "                           |   [V]  |  [V]   |   [V]  |    [um]     |     [um]     | [MOhm]  |    [pA]    | [cm/sec] |       [V]       |" <<endl;
   cout << "                           |        |        |        | stad | rot  | stad  | rot  |         |            |          |                 |" <<endl ;
   cout << "     Spec Limit            |  >150  |  >150  |  >150  |             |              |         |            |          |                 |" <<endl ;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(8) <<  setprecision(3) << Diel_BreakDown_Voltage[6*id]
	       << setw(9) << Diel_BreakDown_Voltage[6*id+1] <<  setw(9)  <<  Diel_BreakDown_Voltage[6*id+2]
	       << setw(8) <<  setprecision(4) << LineWidthStop_s_Rsh_Rgeom_Corrected_E[id]  <<  setw(7) << LineWidthStop_r_Rsh_Rgeom_Corrected_E[id]
	       << setw(8) <<  setprecision(4) << LineWidthStrip_s_Rsh_Rgeom_Corrected_E[id] <<  setw(7) << LineWidthStrip_r_Rsh_Rgeom_Corrected_E[id]
	       << setw(9) <<  setprecision(4) << RPoly_E[id]
	       <<  setw(11) << setprecision(3)<< GCD_I_surface[2*id]
	       <<  setw(12) << GCD_S_interface_recombination_Velocity[2*id]
	       <<  setw(14) << GCD_VFB_Derivative[2*id]
	       << endl;


         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(8) <<  setprecision(3) << Diel_BreakDown_Voltage[6*id+3]
	       << setw(9) << Diel_BreakDown_Voltage[6*id+4] <<  setw(9)  <<  Diel_BreakDown_Voltage[6*id+5]
	       << setw(8) <<  setprecision(4) << LineWidthStop_s_Rsh_Rgeom_Corrected_W[id]  <<  setw(7) << LineWidthStop_r_Rsh_Rgeom_Corrected_W[id]
	       << setw(8) <<  setprecision(4) << LineWidthStrip_s_Rsh_Rgeom_Corrected_W[id] <<  setw(7) << LineWidthStrip_r_Rsh_Rgeom_Corrected_W[id]
	       << setw(9) <<  setprecision(4) << RPoly_W[id]
	       << setw(11) << setprecision(3)<< GCD_I_surface[2*id+1]
	       << setw(12) << GCD_S_interface_recombination_Velocity[2*id+1]
	       << setw(14) << GCD_VFB_Derivative[2*id+1]
	       << endl;
   }

   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   CVS_Output_File<< "     Flute2  Id            ;Mode ; DielNE ; DielNW ; DielSW  ;  LineW_Pstop ;  ;LineW_Strip ;  ; RPoly ; GCD I_surf   ;  GCD S_0   ; GCD Vfb_Acc_Dep" << endl;
   CVS_Output_File<< "                           ;     ;   [V]    ; [V]    ; [V]     ;  [um]        ;  ;[um]        ;   ;[MOhm] ;   [pA]       ;  [cm/sec]  ;    [V]      "<< endl;
   CVS_Output_File<< "                           ;     ;        ;        ;         ; stad    ; rot     ; stad  ;    rot ;   ;       ;              ;            ; " << endl;
   CVS_Output_File<< "     Spec Limit            ;     ;>150    ;>150    ; >150    ;         ;       ;        ;      ;       ;              ;            ;  " << endl;
   for(int id=0;id<nFiles;id++) {
         CVS_Output_File  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) << ";" << "Quick" << ";"
	                  << Diel_BreakDown_Voltage[6*id] << ";"  << Diel_BreakDown_Voltage[6*id+1] <<  ";"  <<  Diel_BreakDown_Voltage[6*id+2] << ";"
	                  << LineWidthStop_s_Rsh_Rgeom_Corrected_E[id]  << ";" << LineWidthStop_r_Rsh_Rgeom_Corrected_E[id] << ";"
	                  << LineWidthStrip_s_Rsh_Rgeom_Corrected_E[id] << ";" << LineWidthStrip_r_Rsh_Rgeom_Corrected_E[id] << ";"
	                  << RPoly_E[id] << ";" << GCD_I_surface[2*id]  << ";" << GCD_S_interface_recombination_Velocity[2*id]  << ";" << GCD_VFB_Derivative[2*id]
	                  << endl;


          CVS_Output_File <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << ";" << "Quick" << ";"
	                  << Diel_BreakDown_Voltage[6*id+3] << ";"  << Diel_BreakDown_Voltage[6*id+4] <<  ";"  <<  Diel_BreakDown_Voltage[6*id+5] << ";"
	                  << LineWidthStop_s_Rsh_Rgeom_Corrected_W[id]  << ";" << LineWidthStop_r_Rsh_Rgeom_Corrected_W[id] << ";"
	                  << LineWidthStrip_s_Rsh_Rgeom_Corrected_W[id] << ";" << LineWidthStrip_r_Rsh_Rgeom_Corrected_W[id] << ";"
	                  << RPoly_W[id] << ";" << GCD_I_surface[2*id+1]  << ";" << GCD_S_interface_recombination_Velocity[2*id+1]  << ";" << GCD_VFB_Derivative[2*id+1]
  	                  << endl;
   }

}
