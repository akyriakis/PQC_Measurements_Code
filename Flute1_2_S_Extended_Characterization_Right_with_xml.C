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

// ------------------------------------------------------------------------------------
// VDP arrays and parameters declartion
//--------------------------------------------------------------------------------
TGraph *VDP_cv[600];
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
//int HM_Id = 37400; // batch number
//string Structure_Id[40] = {"06","19","49"}; //wafers to be analyzed
//int HM_Id = 37401;
//string Structure_Id[40] = {"09","26","41"}; 
//int HM_Id = 37408;
//string Structure_Id[40] = {"10","21","46"}; 
int HM_Id = 37897;
string Structure_Id[40] = {"05","18","43"}; 


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
   if((id % 12) == 0) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPPoly_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 1) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPPoly_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 2) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPPoly_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 3) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPPoly_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 4) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPStop_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 5) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPStop_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 6) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPStop_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 7) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPStop_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 8) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPStrip_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 9) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPStrip_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 10) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPStrip_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 12) == 11) VDP_cv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPStrip_r",HM_Id,Structure_Id[file_id].c_str()));
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
    VDP_DataFileName[12*i+0] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPPoly_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPPoly_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPPoly_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPPoly_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+4] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPStop_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+5] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPStop_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+6] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPStop_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+7] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPStop_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+8] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPStrip_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+9] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPStrip_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_DataFileName[12*i+10] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Right_PQC1_VDPStrip_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_DataFileName[12*i+11] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Right_PQC1_VDPStrip_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
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

   if((id % 2) == 0) FET_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Right_PQC1_FET",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) FET_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Right_PQC1_FET",HM_Id,Structure_Id[file_id].c_str()));

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
    Analysis_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Right_PQC1_FETPSs_trans.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    Analysis_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Right_PQC1_FETPSs_trans.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
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

     FET_cc_final[i]->BuildLegend(0.15,0.5,0.65,0.85);


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
// ========================================================================================================================================

void Flute1_2_S_Extended_Characterization_Right_with_xml(int nFiles)
{
   VDP_PQC_Flute1_Analysis_Final(nFiles);
   VDP_PQC_xml_production(xml_on, nFiles);

   FET_PQC_Analysis_Final(nFiles);
   FET_PQC_xml_production(xml_on, nFiles);

   CVS_Output_File.open(Form("CSV_Info/VPX%d/Flute1_VPX%d_0xx_2-S_HM_E_W_Right_Extended_Info.csv",HM_Id,HM_Id));

   
   cout << endl << endl;
   cout << "=========================================================================================================" << endl;
   cout << "---------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Right",HM_Id) << "  --------------------------------------------------" << endl;
   cout << "==============================================================================================================================================" << endl;
   cout << "     Flute1  Id            |   VdPPoly   |   VdPStop   |  VdPStrip   | FET_Vth | " << endl; 
   cout << "                           |  [kOhm/sq]  |  [kOhm/sq]  |  [Ohm/sq]   |   [V]   |" << endl; 
   cout << "                           |  stad | rot | stad |  rot | stad |  rot |         | " << endl; 
   cout << "     Spec Limit            |  2.2- |2.2- |18-20 |18-20 |32-35 |32-35 |         |" << endl; 
   cout << "     Spec Limit            |   2.4 | 2.4 |      |      |      |      |         |" << endl; 
   cout << "------------------------------------------------------------------------------------------ "<< endl;
   cout << showpoint;  
   for(int id=0;id<nFiles;id++) {  
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Right",HM_Id,Structure_Id[id].c_str())	       
	       << setw(7) <<  setprecision(3) << VDPoly_s_Rsh_E[id]/1000.  <<  setw(7) 
	       << VDPoly_r_Rsh_E[id]/1000. <<  setw(7) <<  setprecision(4) << VDPstop_s_Rsh_E[id]/1000.  <<  setw(7) 
	       << VDPstop_r_Rsh_E[id]/1000.  <<  setw(7)  << VDPStrip_s_Rsh_E[id] <<  setw(7)  << VDPStrip_r_Rsh_E[id]
	       << setw(7)  <<  setprecision(3) << FET_Gate_Voltage_Threshold[2*id]
	       << endl; 
	          
   
         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Right",HM_Id,Structure_Id[id].c_str())
	       << setw(7) <<  setprecision(3) << VDPoly_s_Rsh_W[id]/1000.  <<  setw(7) 
	       << VDPoly_r_Rsh_W[id]/1000. <<  setw(7) <<  setprecision(4) << VDPstop_s_Rsh_W[id]/1000.  <<  setw(7) 
	       << VDPstop_r_Rsh_W[id]/1000.  <<  setw(7)  << VDPStrip_s_Rsh_W[id] <<  setw(7)  << VDPStrip_r_Rsh_W[id] 
	       << setw(7)  <<  setprecision(3) << FET_Gate_Voltage_Threshold[2*id+1] 
	       << endl; 
   }   
 
   CVS_Output_File<< "     Flute1  Id            ;Mode ;   VdPPoly;   ;   VdPStop;   ;  VdPStrip;   ; FET_Vth" << endl; 
   CVS_Output_File<< "                           ; ; [kOhm/sq];  ;  [kOhm/sq];  ;  [Ohm/sq];   ;   [V]  " << endl; 
   CVS_Output_File<< "                           ; ; stad ; rot ; stad ;  rot ; stad ;  rot ;        " << endl; 
   CVS_Output_File<< "     Spec Limit            ; ; 2.2-2.4 ;2.2-2.4 ;18- 20 ;18 - 20 ;32 - 35 ;32 - 35 ;        " << endl; 
   
   for(int id=0;id<nFiles;id++) {  
         CVS_Output_File  <<  Form("VPX%d_0%s_2-S_HM_E_Right",HM_Id,Structure_Id[id].c_str())<< ";" << "Extended" << ";"	       
                          << VDPoly_s_Rsh_E[id]/1000.   <<  ";" << VDPoly_r_Rsh_E[id]/1000.  << ";"
	                  << VDPstop_s_Rsh_E[id]/1000.  <<  ";" << VDPstop_r_Rsh_E[id]/1000. << ";"
	                  << VDPStrip_s_Rsh_E[id] <<  ";"  << VDPStrip_r_Rsh_E[id] << ";" << FET_Gate_Voltage_Threshold[2*id] 
                	  << endl; 
	          
   
         CVS_Output_File  <<  Form("VPX%d_0%s_2-S_HM_W_Right",HM_Id,Structure_Id[id].c_str()) << ";" << "Extended" << ";"
	                  << VDPoly_s_Rsh_W[id]/1000.   <<  ";" << VDPoly_r_Rsh_W[id]/1000.  << ";"
	                  << VDPstop_s_Rsh_W[id]/1000.  <<  ";" << VDPstop_r_Rsh_W[id]/1000. << ";"
	                  << VDPStrip_s_Rsh_W[id] <<  ";"  << VDPStrip_r_Rsh_W[id] << ";" << FET_Gate_Voltage_Threshold[2*id+1] 
	                  << endl; 
   }   

    CVS_Output_File.close();

}
