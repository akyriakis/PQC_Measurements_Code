
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



double epsilon_0 = 0.0885418782; // [pF/cm]
double epsilon_Si = 11.68*epsilon_0; // [pF/cm]
double epsilon_SiO2 = 3.9*epsilon_0; // [pF/cm]
//double MOS_area;
double q= 1.60219E-7; // [pF.V]
double kB = 8.6173324E-5 * q; // [pF*V^2/K]
double TKelvin;
double Nintrisic;

// ------------------------------------------------------

float VDPBulk_current_d[500][500];
float VDPBulk_voltage_d[500][500];
float VDPBulk_Temp_measurement[600];
float VDPBulk_Humid_measurement[600];
double VDPBulk_A_constant[600];
double VDPBulk_A_slope[600];
float VDPBulk_cv_low[600];
float VDPBulk_cv_high[600];
double VDPBulk_s_Rsh_E[100];
double VDPBulk_r_Rsh_E[100];
double VDPBulk_s_Rsh_W[100];
double VDPBulk_r_Rsh_W[100];
double VDPBulk_Rsh_all[100];
double VDPBulk_Resistivity_all[100];
double VDPBulk_Resistivity_s_F1_E[100];
double VDPBulk_Resistivity_s_F2_E[100];
double VDPBulk_Resistivity_s_F1_W[100];
double VDPBulk_Resistivity_s_F2_W[100];
double VDPBulk_Resistivity_r_F1_E[100];
double VDPBulk_Resistivity_r_F2_E[100];
double VDPBulk_Resistivity_r_F1_W[100];
double VDPBulk_Resistivity_r_F2_W[100];
TGraph *VDPBulk_iv[600];
vector<vector <string> > VDPBulk_timestamp_meas;
vector<vector <double> > VDPBulk_temp_meas;
vector<vector <double> > VDPBulk_air_temp_meas;
vector<vector <double> > VDPBulk_rh_prcnt_meas;
vector<std::string> VDPBulk_Operators;
vector<std::string> VDPBulk_Begin_Timestamp;
vector<std::string> VDPBulk_name_labels;
vector<std::string> VDPBulk_kind_of_parts;
vector<std::string> VDPBulk_kind_of_HM_flute_id;
vector<std::string> VDPBulk_struct_id;
vector<std::string> VDPBulk_set_id;
vector<std::string> VDPBulk_config_id;
vector<std::string> VDPBulk_equipment;
vector<double> VDPBulk_waiting_time;
vector<double> VDPBulk_temp;
vector<double> VDPBulk_av_temp;
vector<double> VDPBulk_Nmeas;


// ------------------------------------------------------

float Diode_IV_vv[500][500];
float Diode_IV_Curr[500][500];
double Diode_IV_vv_d[500][500];
double Diode_IV_Curr_d[500][500];
float Diode_IV_Temp_measurement[500];
float Diode_IV_Humid_measurement[500];
TGraph *Diode_IV[100];
float Breakdown_Voltage[100];
float Diode_Current_at_600V[100];
float Diode_Current_at_600V_over_Volume[100];
float Diode_Current_at_600V_times_volume[100];
vector<vector <string> > Diode_IV_timestamp_meas;
vector<vector <double> > Diode_IV_temp_meas;
vector<vector <double> > Diode_IV_air_temp_meas;
vector<vector <double> > Diode_IV_rh_prcnt_meas;
vector<std::string> Diode_IV_Operators;
vector<std::string> Diode_IV_Begin_Timestamp;
vector<std::string> Diode_IV_name_labels;
vector<std::string> Diode_IV_kind_of_parts;
vector<std::string> Diode_IV_struct_id;
vector<std::string> Diode_IV_set_id;
vector<std::string> Diode_IV_config_id;
vector<std::string> Diode_IV_equipment;
vector<double> Diode_IV_waiting_time;
vector<double> Diode_IV_bias;
vector<double> Diode_IV_temp;
vector<double> Diode_IV_av_temp;
vector<double> Diode_IV_compliance;
vector<double> Diode_IV_Nmeas;

// ------------------------------------------------------

float Diode_CV_vv[500][500];
float Diode_CV_Cap[500][500];
double Diode_CV_vv_d[500][500];
double Diode_CV_Cap_d[500][500];
float Inv_Diode_CV_Cap[500][500];
double Inv_Diode_CV_Cap_d[500][500];
double Diode_CV_dev_spline_InvCV2[500];
float Diode_CV_Temp_measurement[500];
float Diode_CV_Humid_measurement[500];
TGraph *Diode_CV_cv[40];
TGraph *Diode_CV_invC2_vs_V[40];
TSpline3 * Diode_CV_CuSpl_InvCV2[40];
float Diode_CV_FullDepletion_Voltage[40];
float Diode_CV_Minimum_Capacitance[40];
float Diode_CV_Bulk_Resistivity[40];
float Diode_CV_Bulk_Concentration[40];
float Diode_CV_active_thickness[40];
vector<vector <string> > Diode_CV_timestamp_meas;
vector<vector <double> > Diode_CV_rsstnc_mhom;
vector<vector <double> > Diode_CV_temp_meas;
vector<vector <double> > Diode_CV_air_temp_meas;
vector<vector <double> > Diode_CV_rh_prcnt_meas;
vector<std::string> Diode_CV_Operators;
vector<std::string> Diode_CV_Begin_Timestamp;
vector<std::string> Diode_CV_name_labels;
vector<std::string> Diode_CV_kind_of_parts;
vector<std::string> Diode_CV_kind_of_HM_flute_id;
vector<std::string> Diode_CV_struct_id;
vector<std::string> Diode_CV_set_id;
vector<std::string> Diode_CV_config_id;
vector<std::string> Diode_CV_equipment;
vector<double> Diode_CV_waiting_time;
vector<double> Diode_CV_ac_freq;
vector<double> Diode_CV_ac_ampl;
vector<double> Diode_CV_temp;
vector<double> Diode_CV_av_temp;
vector<double> Diode_CV_compliance;
vector<double> Diode_CV_Nmeas;


ofstream CVS_Output_File;

TString Data_Dir = "/home/akyriakis/MOS-Measurements_Analysis/PQC_Measurements/Data_New_Format/";

TString Meander_DataFileName[600];
TString CloverMetal_DataFileName[600];
TString VDPEdge_DataFileName[600];
TString LineWidthEdge_DataFileName[600];
TString VDPBulk_DataFileName[600];
TString Diode_IV_DataFileName[100];
TString Diode_CV_DataFileName[100];


//int HM_Id = 34355;
//string Structure_Id[100] = {"02","03","04","05","06","10","12","14","16","22","26","30","31","41","42","43","45","46","47","48"};
//int HM_Id = 34356;
//string Structure_Id[100] = {"01","02","03","06","14","15","20","22","23","24","30","31","34","46","47"};
//int HM_Id = 35715;
//string Structure_Id[40] = {"10","21","34","45"};
//int HM_Id = 35716;
//string Structure_Id[40] = {"04","17","31","45"};
//int HM_Id = 35720;
//string Structure_Id[40] = {"08","15","34","49"};
//int HM_Id = 35721;
//string Structure_Id[40] = {"02","12","20","33"};
//int HM_Id = 36243;
//string Structure_Id[40] = {"05","15","26","35"};
//int HM_Id = 36244;
//string Structure_Id[40] = {"03","16","28","42"};
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

// ------------------------------------------------------------------------
double VDPBulk_Read_Info_from_File(TString VDPBulk_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(VDPBulk_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << VDPBulk_DataFileName << " wit Id = " << id << endl;

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
         VDPBulk_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         VDPBulk_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         VDPBulk_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         VDPBulk_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         VDPBulk_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         VDPBulk_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         VDPBulk_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20];
         sscanf(line.c_str(), "%*s %s", str1);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1);
         cout<<Kind_of_HM_config_id<<endl;
         VDPBulk_config_id.push_back(Kind_of_HM_config_id);
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
         VDPBulk_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         VDPBulk_waiting_time.push_back(Waiting_time_s);
         }
        if(nlines == 12) {
        cout << line << '\n';
        double temp_degC;
        sscanf(line.c_str(), "%*s %lf", &temp_degC);
        cout<<temp_degC<<endl;
        VDPBulk_temp.push_back(temp_degC);
       }
        if(nlines == 13) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        VDPBulk_av_temp.push_back(Av_temp_degC);
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
     VDPBulk_voltage_d[id][N_meas] = voltage*1E-9;  // to 1E-9 exei mpei gia na diorthwsei to lathos pou yparxei sta txt. tha prepei na aferethei an diorthwthei
     VDPBulk_current_d[id][N_meas] = current_namp;
      VDPBulk_temp_meas.push_back(std::vector<double>());
      VDPBulk_temp_meas[id].push_back(temp_degC);
      VDPBulk_air_temp_meas.push_back(std::vector<double>());
      VDPBulk_air_temp_meas[id].push_back(air_temp_degC);
      VDPBulk_rh_prcnt_meas.push_back(std::vector<double>());
      VDPBulk_rh_prcnt_meas[id].push_back(rh_prcnt);
      VDPBulk_timestamp_meas.push_back(std::vector<string>());
      VDPBulk_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
     //VDPBulk_voltage[id][N_meas] = float(VDPBulk_voltage_d[id][N_meas]);
      //VDPBulk_current[id][N_meas] = float(VDPBulk_current_d[id][N_meas]);


    //        cout << N_meas << ")  "  <<  VDP_vv[N_meas] << "   " << VDP_cap[N_meas] << endl;
    //  cout << N_meas << ")  "  <<  VDP_voltage[id][N_meas] << "   " << VDP_current[id][N_meas]  << endl;
     N_meas++;

      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   //cout  << "----->>>>>>   Number of Measurements = " << N_meas << endl;
   VDPBulk_Nmeas.push_back(N_meas);

   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;

   TGaxis::SetMaxDigits(3);

//   TCanvas *cc_spline = new TCanvas();
   VDPBulk_iv[id] =  new TGraph(N_meas,VDPBulk_current_d[id], VDPBulk_voltage_d[id]);
   if((id % 4) == 0) VDPBulk_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC3_VDPBulk_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 1) VDPBulk_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC3_VDPBulk_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 2) VDPBulk_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC3_VDPBulk_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 3) VDPBulk_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC3_VDPBulk_r",HM_Id,Structure_Id[file_id].c_str()));
   VDPBulk_iv[id]->GetXaxis()->SetTitle("VDPBulk Current [A]");
   VDPBulk_iv[id]->GetYaxis()->SetTitle("VDPBulk Voltage [V]");
   VDPBulk_iv[id]->SetDrawOption("AP");
   VDPBulk_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   VDPBulk_iv[id]->SetMarkerColor(colorid);
//   VDPBulk_iv[id]->Draw();


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   VDPBulk_cv_low[id] = VDPBulk_current_d[id][3];
   VDPBulk_cv_high[id] = VDPBulk_current_d[id][38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", VDPBulk_cv_low[id], VDPBulk_cv_high[id]);
//   VDPBulk_iv[id]->Fit("f_cv","0R+");
   VDPBulk_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   VDPBulk_A_constant[id] = VDPBulk_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   VDPBulk_A_slope[id]  = VDPBulk_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double R_sh = (TMath::Pi()/TMath::Log(2))*VDPBulk_A_slope[id]; // VDPBulk R_sh
   double R_sh_error = (TMath::Pi()/TMath::Log(2))*(VDPBulk_iv[id]->GetFunction("f_cv")->GetParError(1));

   cout << "----->>>   R_sh Voltage = " << R_sh << " +/- " <<  R_sh_error<< " [Ohm/sq]" << endl;

   return R_sh;

}
void VDPBulk_PQC_Flute3_Analysis_Final(int nFiles)
{

 //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1

 //   gROOT->Reset();
   gStyle->SetOptStat(0);



   //TString VDPBulk_DataFileName[600];
   for(int i=0;i<nFiles;i++) {
    VDPBulk_DataFileName[4*i+0] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC3_VDPBulk_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDPBulk_DataFileName[4*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC3_VDPBulk_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDPBulk_DataFileName[4*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC3_VDPBulk_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDPBulk_DataFileName[4*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC3_VDPBulk_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

   }

   float s = 187e-6; // in m
   s *= 100; // in cm
   float F1 = 1.089;
   float F2 = 1.218;

   for (int i=0;i<nFiles;i++) {

     cout << "VDPBulk Structure : Analyze file : " << VDPBulk_DataFileName[4*i+0] << endl;
     VDPBulk_s_Rsh_E[i] = VDPBulk_Read_Info_from_File(VDPBulk_DataFileName[4*i+0],(4*i+0),i);
     VDPBulk_Resistivity_s_F1_E[i] = 2*TMath::Pi()*s*F1*VDPBulk_s_Rsh_E[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());
     VDPBulk_Resistivity_s_F2_E[i] = 2*TMath::Pi()*s*F2*VDPBulk_s_Rsh_E[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());

     cout << "VDPBulk Structure :Analyze file : " << VDPBulk_DataFileName[4*i+1] << endl;
     VDPBulk_s_Rsh_W[i] = VDPBulk_Read_Info_from_File(VDPBulk_DataFileName[4*i+1],(4*i+1),i);
     VDPBulk_Resistivity_s_F1_W[i] = 2*TMath::Pi()*s*F1*VDPBulk_s_Rsh_W[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());
     VDPBulk_Resistivity_s_F2_W[i] = 2*TMath::Pi()*s*F2*VDPBulk_s_Rsh_W[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());

     cout << "VDPBulk Structure :Analyze file : " << VDPBulk_DataFileName[4*i+2] << endl;
     VDPBulk_r_Rsh_E[i] = VDPBulk_Read_Info_from_File(VDPBulk_DataFileName[4*i+2],(4*i+2),i);
     VDPBulk_Resistivity_r_F1_E[i] = 2*TMath::Pi()*s*F1*VDPBulk_s_Rsh_W[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());
     VDPBulk_Resistivity_r_F2_E[i] = 2*TMath::Pi()*s*F2*VDPBulk_s_Rsh_W[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());

     cout << "VDPBulk Structure :Analyze file : " << VDPBulk_DataFileName[4*i+3] << endl;
     VDPBulk_r_Rsh_W[i] = VDPBulk_Read_Info_from_File(VDPBulk_DataFileName[4*i+3],(4*i+3),i);
     VDPBulk_Resistivity_r_F1_W[i] = 2*TMath::Pi()*s*F1*VDPBulk_s_Rsh_W[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());
     VDPBulk_Resistivity_r_F2_W[i] = 2*TMath::Pi()*s*F2*VDPBulk_s_Rsh_W[i]/(2-sqrt(2))*(TMath::Log(2)/TMath::Pi());

     VDPBulk_Rsh_all[4*i+0] = VDPBulk_s_Rsh_E[i];
     VDPBulk_Rsh_all[4*i+1] = VDPBulk_s_Rsh_W[i];
     VDPBulk_Rsh_all[4*i+2] = VDPBulk_r_Rsh_E[i];
     VDPBulk_Rsh_all[4*i+3] = VDPBulk_r_Rsh_W[i];

     VDPBulk_Resistivity_all[4*i+0] = VDPBulk_Resistivity_s_F1_E[i];
     VDPBulk_Resistivity_all[4*i+1] = VDPBulk_Resistivity_s_F1_W[i];
     VDPBulk_Resistivity_all[4*i+2] = VDPBulk_Resistivity_r_F1_E[i];
     VDPBulk_Resistivity_all[4*i+3] = VDPBulk_Resistivity_r_F1_W[i];

   }

   TCanvas *cc_final[4];
   TMultiGraph *mg[4];
   for (int i=0;i<4;i++) {
     cc_final[i] = new TCanvas(Form("cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     mg[i] = new TMultiGraph(Form("mg_%d",i),"VDPBulk R_sh");
     if(i == 0) mg[i]->SetTitle("VDPBulk_s_E R_sh ");
     if(i == 1) mg[i]->SetTitle("VDPBulk_s_W R_sh ");
     if(i == 2) mg[i]->SetTitle("VDPBulk_r_E R_sh ");
     if(i == 3) mg[i]->SetTitle("VDPBulk_r_W R_sh ");

    // VDPBulk_iv[4*i+]->Draw("ALP");
     for(int id=0;id<nFiles;id++) {

        mg[i]->Add(VDPBulk_iv[4*id+i]);

     }
     mg[i]->Draw("LPsame");
     mg[i]->GetXaxis()->SetTitle("Current [A]");
     mg[i]->GetYaxis()->SetTitle("Voltage [V]");
     mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     cc_final[i]->BuildLegend(0.15,0.5,0.5,0.85);


     gPad->Modified();
     gPad->Update();
   }


   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;

   cout << "      VDPBulk   Id                 |  VDPBulk_s_Rsh [kOhm/sq]   | VDPBulk_r_Rsh [kOhm/sq] |  VDPBulk_Resistivity_s_F1 [kOhmcm]|  VDPBulk_Resistivity_s_F2[kOhmcm]" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(20) <<  setprecision(4) << VDPBulk_s_Rsh_E[id]/1000.  <<  setw(30)  << VDPBulk_r_Rsh_E[id]/1000.
	       << setw(20) <<  setprecision(3) << VDPBulk_Resistivity_s_F1_E[id]/1000.  <<  setw(30)  << VDPBulk_Resistivity_s_F2_E[id]/1000.
	       << endl;
         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(20) <<  setprecision(4) << VDPBulk_s_Rsh_W[id]/1000.  <<  setw(30)  << VDPBulk_r_Rsh_W[id]/1000.
	       << setw(20) <<  setprecision(3) << VDPBulk_Resistivity_s_F1_W[id]/1000.  <<  setw(30)  << VDPBulk_Resistivity_s_F2_W[id]/1000.
	       << endl;

   }
  cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
}

//--------------------------------------------------------------------------------------------------//
// VDPBulk xml production
// --------------------------------------------------------------------------------------------------//
void VDPBulk_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  cout << "VDPBulk Xml production is on" << endl;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<4; j++){
          //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = VDPBulk_DataFileName[4*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = VDPBulk_DataFileName[4*i+j].Last(delimeter2);
         xml_filename = VDPBulk_DataFileName[4*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << Diode_DataFileName[4*i+j](0,pos1) << endl;
         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute3/VPX%d/",HM_Id) +xml_filename);

        // First create engine
        TXMLEngine VDPBulk_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = VDPBulk_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = VDPBulk_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = VDPBulk_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = VDPBulk_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = VDPBulk_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = VDPBulk_xml.NewChild(HEADER, 0, "RUN");
              VDPBulk_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              VDPBulk_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              VDPBulk_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              VDPBulk_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              VDPBulk_xml.NewChild(RUN, 0, "INITIATED_BY_USER", VDPBulk_Operators[4*i+j].c_str());
              VDPBulk_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", VDPBulk_Begin_Timestamp[4*i+j].c_str());
              VDPBulk_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = VDPBulk_xml.NewChild(ROOT, 0, "DATA_SET");
            VDPBulk_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            VDPBulk_xml.NewChild(DATA_SET, 0, "VERSION", "IV_measurement");
            XMLNodePointer_t PART = VDPBulk_xml.NewChild(DATA_SET, 0, "PART");
              VDPBulk_xml.NewChild(PART, 0, "NAME_LABEL", VDPBulk_name_labels[4*i+j].c_str());
              VDPBulk_xml.NewChild(PART, 0, "KIND_OF_PART", VDPBulk_kind_of_parts[4*i+j].c_str());
            XMLNodePointer_t DATA = VDPBulk_xml.NewChild(DATA_SET, 0, "DATA");
              VDPBulk_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC3");
              VDPBulk_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", VDPBulk_struct_id[4*i+j].c_str());
              VDPBulk_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", VDPBulk_set_id[4*i+j].c_str());
              VDPBulk_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", VDPBulk_config_id[4*i+j].c_str());
              VDPBulk_xml.NewChild(DATA, 0, "EQUIPMENT", VDPBulk_equipment[4*i+j].c_str());
              VDPBulk_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(VDPBulk_waiting_time[4*i+j], 3).c_str());
              VDPBulk_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(VDPBulk_temp[4*i+j], 3).c_str());
              VDPBulk_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(VDPBulk_av_temp[4*i+j], 3).c_str());
              VDPBulk_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", VDPBulk_struct_id[4*i+j].c_str());
              VDPBulk_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = VDPBulk_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = VDPBulk_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = VDPBulk_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  VDPBulk_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  VDPBulk_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = VDPBulk_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                VDPBulk_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                VDPBulk_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "IV_measurement");
                XMLNodePointer_t PART_CDS1 = VDPBulk_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  VDPBulk_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", VDPBulk_name_labels[4*i+j].c_str());
                  VDPBulk_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", VDPBulk_kind_of_parts[4*i+j].c_str());
                for (int k=0; k<VDPBulk_Nmeas[4*i+j]-1; k++){
                    XMLNodePointer_t DATA_CDS1 = VDPBulk_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                      VDPBulk_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string_scientific_format(VDPBulk_voltage_d[4*i+j][k], 3).c_str());
                      VDPBulk_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(VDPBulk_current_d[4*i+j][k]*1E+9, 3).c_str());
                      VDPBulk_xml.NewChild(DATA_CDS1, 0, "TIME", VDPBulk_timestamp_meas[4*i+j][k].c_str());
                      VDPBulk_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(VDPBulk_temp_meas[4*i+j][k], 3).c_str());
                      VDPBulk_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(VDPBulk_air_temp_meas[4*i+j][k], 3).c_str());
                      VDPBulk_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(VDPBulk_rh_prcnt_meas[4*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = VDPBulk_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = VDPBulk_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = VDPBulk_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      VDPBulk_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      VDPBulk_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = VDPBulk_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    VDPBulk_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    VDPBulk_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement");
                    XMLNodePointer_t PART_CDS2 = VDPBulk_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      VDPBulk_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", VDPBulk_name_labels[4*i+j].c_str());
                      VDPBulk_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", VDPBulk_kind_of_parts[4*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = VDPBulk_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      VDPBulk_xml.NewChild(DATA_CDS2, 0, "RSH_OHMSQR", Conv_float_to_string_scientific_format(VDPBulk_Rsh_all[4*i+j], 3).c_str());
                      VDPBulk_xml.NewChild(DATA_CDS2, 0, "R_OHM", Conv_float_to_string_scientific_format(VDPBulk_Rsh_all[4*i+j]*(TMath::Log(2)/TMath::Pi()), 3).c_str());
                      VDPBulk_xml.NewChild(DATA_CDS2, 0, "RHO_KOHMCM", Conv_float_to_string(VDPBulk_Resistivity_all[4*i+j]*1E-3, 3).c_str());

        XMLDocPointer_t VDPBulk_xmldoc = VDPBulk_xml.NewDoc();
        VDPBulk_xml.DocSetRootElement(VDPBulk_xmldoc, ROOT);
        // Save document to file
        VDPBulk_xml.SaveDoc(VDPBulk_xmldoc, Form("./XML_Info_Flute3/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        VDPBulk_xml.FreeDoc(VDPBulk_xmldoc);
      }
    }
  }
}

// ------------------------------------------------------------------------

void Diode_IV_Read_Info_from_File(TString Diode_IV_DataFileName, int id, int file_id)
{

   float Diode_half_area = 0.25*0.25; // in cm^2
   float Diode_width = 0.03; // in cm = 300um
   float Diode_Volume = Diode_half_area*Diode_width; //in cm^3
   float Diode_Volume_in_mm = Diode_half_area*Diode_width*1000; //in mm^3

   float Diode_Current_600V = 0;
   ifstream in;
   in.open(Diode_IV_DataFileName);

   int nlines = 0;
   int N_meas = 0;

   string line;

   //char date[20],time[20];
   //float Diode_Voltage, Diode_Current, Diode_Temperature, Diode_AirTemperature, Diode_Hymidity;
   while (in.good()) {
     getline (in,line);
      if(nlines < 18) {
         if(nlines == 1) {
         //cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         Diode_IV_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         Diode_IV_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         Diode_IV_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         Diode_IV_kind_of_parts.push_back(Kind_of_part);
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
         Diode_IV_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         Diode_IV_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         Diode_IV_config_id.push_back(Kind_of_HM_config_id);
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
         Diode_IV_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         Diode_IV_waiting_time.push_back(Waiting_time_s);
         }
         if(nlines == 12) {
         cout << line << '\n';
         double Bias_V;
         sscanf(line.c_str(), "%*s %lf", &Bias_V);
         cout<<Bias_V<<endl;
         Diode_IV_bias.push_back(Bias_V);
         }
         if(nlines == 13) {
         cout << line << '\n';
         double Temp_set_degC;
         sscanf(line.c_str(), "%*s %lf", &Temp_set_degC);
         cout<<Temp_set_degC<<endl;
         Diode_IV_temp.push_back(Temp_set_degC);
        }
        if(nlines == 14) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        Diode_IV_av_temp.push_back(Av_temp_degC);
       }
        if(nlines == 15) {
        cout << line << '\n';
        double ComplianceA;
        sscanf(line.c_str(), "%*s %lf", &ComplianceA);
        cout<<ComplianceA<<endl;
        Diode_IV_compliance.push_back(ComplianceA);
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

    	Diode_IV_vv_d[id][N_meas] = -voltage;
    	Diode_IV_Curr_d[id][N_meas] = current_namp;
      Diode_IV_temp_meas.push_back(std::vector<double>());
      Diode_IV_temp_meas[id].push_back(temp_degC);
      Diode_IV_air_temp_meas.push_back(std::vector<double>());
      Diode_IV_air_temp_meas[id].push_back(air_temp_degC);
      Diode_IV_rh_prcnt_meas.push_back(std::vector<double>());
      Diode_IV_rh_prcnt_meas[id].push_back(rh_prcnt);
      Diode_IV_timestamp_meas.push_back(std::vector<string>());
      Diode_IV_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
    	Diode_IV_vv[id][N_meas] = (float) Diode_IV_vv_d[id][N_meas];
      Diode_IV_Curr[id][N_meas] = (float) Diode_IV_Curr_d[id][N_meas];

      if(Diode_IV_vv_d[id][N_meas] == 600) Diode_Current_600V = Diode_IV_Curr_d[id][N_meas];
      if((Diode_IV_Curr_d[id][N_meas]-Diode_IV_Curr_d[id][N_meas-1])> 500){
        Breakdown_Voltage[id] =  Diode_IV_vv_d[id][N_meas];
      }else{
        Breakdown_Voltage[id] =  1000;
      }
      //cout<< Diode_IV_vv[N_meas] << "," << Diode_IV_Curr[N_meas] << endl;
      N_meas++;
      }
      if (!in.good()) break;
      nlines++;
   }
   //cout << N_meas << ")  "  <<  Diode_IV_vv[N_meas] << "   " << Diode_IV_Curr[N_meas]  << endl;
   printf(" Found %d lines\n",nlines);
   in.close();

   cout  << "----->>>>>>   Number of Measurements = " << N_meas << endl;
   Diode_IV_Nmeas.push_back(N_meas);


   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;


   //TCanvas *cc = new TCanvas("cc","Diode Current vs Diode_Voltage",500,50,500,500);
   Diode_IV[id] =  new TGraph(N_meas, Diode_IV_vv_d[id], Diode_IV_Curr_d[id]);

   if((id % 2) == 0) Diode_IV[id]->SetTitle(Form("DiodeE%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) Diode_IV[id]->SetTitle(Form("DiodeW%d_0%s IV",HM_Id,Structure_Id[file_id].c_str()));

   Diode_IV[id]->GetXaxis()->SetTitle("Diode_Voltage [V]");
   Diode_IV[id]->GetYaxis()->SetTitle("Diode_Current [nA]");
   Diode_IV[id]->SetDrawOption("AP");
   Diode_IV[id]->SetMarkerStyle(20+id);
   int colorid = 1+id;
   if(colorid == 5) colorid = 46;
   Diode_IV[id]->SetMarkerColor(colorid);
   //Diode_IV[id]->Draw();

   Diode_Current_at_600V[id] = Diode_Current_600V*1E+3;
   cout << "----->>>  Diode_Current at 600V = " << Diode_Current_600V*1E+3  << " pA"<<  endl;
   Diode_Current_at_600V_over_Volume[id] = Diode_Current_600V/Diode_Volume_in_mm;
   cout << "----->>>  Diode_Current at 600V per Volume = " << Diode_Current_600V / Diode_Volume_in_mm << " nA/mm^3"<<  endl;


}
void Diode_IV_PQC3_Analysis_Final(int nFiles)
{

 //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1

 //   gROOT->Reset();
   gStyle->SetOptStat(0);

   //TString Diode_IV_DataFileName[100];
   for(int i=0;i<nFiles;i++) {
    Diode_IV_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Left_PQC3_DiodeHalf_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    Diode_IV_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Left_PQC3_DiodeHalf_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {
     cout << endl;
     cout << "Analyze file : " << Diode_IV_DataFileName[i] << endl;
     Diode_IV_Read_Info_from_File(Diode_IV_DataFileName[2*i],2*i,i);
     Diode_IV_Read_Info_from_File(Diode_IV_DataFileName[2*i+1],2*i+1,i);
  }

   TCanvas *cc_final[2];
   TMultiGraph *mg[2];
   for (int i=0;i<2;i++) {
     cc_final[i] = new TCanvas(Form("cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     mg[i] = new TMultiGraph(Form("mg_%d",i),"Diode  Current vs V");
     if(i == 0) mg[i]->SetTitle("Diode E Current vs V");
     if(i == 1) mg[i]->SetTitle("Diode W Current vs V");

     for(int id=0;id<nFiles;id++) {

        mg[i]->Add(Diode_IV[2*id+i]);

     }
     mg[i]->Draw("LPsame");
     mg[i]->GetXaxis()->SetTitle("Diode_Voltage [V]");
     mg[i]->GetYaxis()->SetTitle("Diode_Current [A]");
     mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     cc_final[i]->BuildLegend(0.55,0.15,0.85,0.55);


     gPad->Modified();
     gPad->Update();
   }

   cout <<endl;
   cout << "===========================+==================================================================" << endl;
   cout << "---------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  ---------------------------------" << endl;
   cout << "==============================================================================================" << endl;

   cout << "       PQC3  Diode Id                 |  Diode_Current [pA] " <<endl;
   cout << "----------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
      cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
            << setw(30)  <<  setprecision(3) <<  Diode_Current_at_600V_over_Volume[2*id] << endl;

      cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
            << setw(30)  <<  setprecision(3) <<  Diode_Current_at_600V_over_Volume[2*id+1] << endl;

   }
   cout << "----------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------------------------------------------//
// Diode_IV xml production
// --------------------------------------------------------------------------------------------------//
void Diode_IV_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    cout << "Diode IV Xml production is on" << endl;
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<=1; j++){
  //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = Diode_IV_DataFileName[2*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = Diode_IV_DataFileName[2*i+j].Last(delimeter2);
         xml_filename = Diode_IV_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
 //        cout << Diode_DataFileName[2*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute3/VPX%d/",HM_Id) +xml_filename);


         // First create engine
        TXMLEngine Diode_IV_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = Diode_IV_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = Diode_IV_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = Diode_IV_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = Diode_IV_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = Diode_IV_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = Diode_IV_xml.NewChild(HEADER, 0, "RUN");
              Diode_IV_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              Diode_IV_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              Diode_IV_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              Diode_IV_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              Diode_IV_xml.NewChild(RUN, 0, "INITIATED_BY_USER", Diode_IV_Operators[2*i+j].c_str());
              Diode_IV_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", Diode_IV_Begin_Timestamp[2*i+j].c_str());
              Diode_IV_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = Diode_IV_xml.NewChild(ROOT, 0, "DATA_SET");
            Diode_IV_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            Diode_IV_xml.NewChild(DATA_SET, 0, "VERSION", "v2");
            XMLNodePointer_t PART = Diode_IV_xml.NewChild(DATA_SET, 0, "PART");
              Diode_IV_xml.NewChild(PART, 0, "NAME_LABEL", Diode_IV_name_labels[2*i+j].c_str());
              Diode_IV_xml.NewChild(PART, 0, "KIND_OF_PART", Diode_IV_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = Diode_IV_xml.NewChild(DATA_SET, 0, "DATA");
              Diode_IV_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC3");
              Diode_IV_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", Diode_IV_struct_id[2*i+j].c_str());
              Diode_IV_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", Diode_IV_set_id[2*i+j].c_str());
              Diode_IV_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", Diode_IV_config_id[2*i+j].c_str());
              Diode_IV_xml.NewChild(DATA, 0, "EQUIPMENT", Diode_IV_equipment[2*i+j].c_str());
              Diode_IV_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(Diode_IV_waiting_time[2*i+j],3).c_str());
              Diode_IV_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(Diode_IV_temp[2*i+j],3).c_str());
              Diode_IV_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(Diode_IV_av_temp[2*i+j],3).c_str());
              Diode_IV_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "Diode_IV1");
              Diode_IV_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = Diode_IV_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = Diode_IV_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = Diode_IV_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  Diode_IV_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  Diode_IV_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = Diode_IV_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                Diode_IV_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                Diode_IV_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "V2");
                XMLNodePointer_t PART_CDS1 = Diode_IV_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  Diode_IV_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", Diode_IV_name_labels[2*i+j].c_str());
                  Diode_IV_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", Diode_IV_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<Diode_IV_Nmeas[2*i+j]-1; k++){
                XMLNodePointer_t DATA_CDS1 = Diode_IV_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  Diode_IV_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(Diode_IV_vv[2*i+j][k],3).c_str());
                  //Diode_IV_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", to_string(Diode_IV_Curr[k]).c_str());
                  Diode_IV_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(Diode_IV_Curr[2*i+j][k],3).c_str());
                  Diode_IV_xml.NewChild(DATA_CDS1, 0, "TIME", Diode_IV_timestamp_meas[2*i+j][k].c_str());
                  Diode_IV_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(Diode_IV_temp_meas[2*i+j][k], 3).c_str());
                  Diode_IV_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(Diode_IV_air_temp_meas[2*i+j][k], 3).c_str());
                  Diode_IV_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(Diode_IV_rh_prcnt_meas[2*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = Diode_IV_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = Diode_IV_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = Diode_IV_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      Diode_IV_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      Diode_IV_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = Diode_IV_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    Diode_IV_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    Diode_IV_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement-004");
                    XMLNodePointer_t PART_CDS2 = Diode_IV_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      Diode_IV_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", Diode_IV_name_labels[2*i+j].c_str());
                      Diode_IV_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", Diode_IV_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = Diode_IV_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      Diode_IV_xml.NewChild(DATA_CDS2, 0, "VBD_V", Conv_float_to_string(Breakdown_Voltage[2*i+j], 3).c_str());
                      Diode_IV_xml.NewChild(DATA_CDS2, 0, "I600_PAMPR", Conv_float_to_string(Diode_Current_at_600V[2*i+j], 3).c_str());
                      Diode_IV_xml.NewChild(DATA_CDS2, 0, "VI600_NAMPRMM3", Conv_float_to_string(Diode_Current_at_600V_over_Volume[2*i+j], 3).c_str());
        XMLDocPointer_t Diode_IV_xmldoc = Diode_IV_xml.NewDoc();
        Diode_IV_xml.DocSetRootElement(Diode_IV_xmldoc, ROOT);
        // Save document to file
        Diode_IV_xml.SaveDoc(Diode_IV_xmldoc, Form("./XML_Info_Flute3/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        Diode_IV_xml.FreeDoc(Diode_IV_xmldoc);
      }
    }
  }
}
// ------------------------------------------------------------------------
void Read_CV_from_File(TString Diode_CV_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(Diode_CV_DataFileName);
   //double TKelvin, Nintrisic;

   float Diode_half_area = 0.25*0.25; // in cm^2    Area=0.250*0.250 -math.pi*(math.pow((0.05/2.0),2)) #cm^2 Diode Area of flute1
   float Diode_width = 0.0290; // in cm = 300um
   float Diode_Volume = Diode_half_area*Diode_width; //in cm^3
   float Diode_Volume_in_mm = Diode_half_area*Diode_width*1000; //in mm^3

   //epsilon_0 = 8.85418782E-14; // [F/cm]
   //epsilon_Si = 11.68;
   double sensor_width = 290E-4; // [cm]
   double hole_mobility = 484;  // [cm^2/(V*sec)]

   int nlines = 0;
   int N_meas = 0;

   string line;
   //char date[20],time[20];
   //float Diode_Voltage, Diode_Capacitance, Diode_Resistance, Diode_Temperature, Diode_AirTemperature, Diode_Hymidity;
   while (in.good()) {
     getline (in,line);
      if(nlines < 17) {
         if(nlines == 1) {
           cout << line << '\n';
           char First_Name[20], Last_Name[20];
           sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
           cout<<First_Name<<Last_Name<<endl;
           Diode_CV_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
           cout << line << '\n';
           char Date[20], Time[20];
           sscanf(line.c_str(), "%*s %s %s", Date, Time);
           cout<<Date<<Time<<endl;
           Diode_CV_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
           cout << line << '\n';
           char Name_Label[25];
           sscanf(line.c_str(), "%*s %s", Name_Label);
           cout<<Name_Label<<endl;
           Diode_CV_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
           cout << line << '\n';
           char batch_type[20], loc[20];
           sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
           string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
           cout<<Kind_of_part<<endl;
           Diode_CV_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
           cout << line << '\n';
           char Kind_of_HM_flute_id[20];
           sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
           cout<<Kind_of_HM_flute_id<<endl;
           Diode_CV_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
           cout << line << '\n';
           char Kind_of_HM_struct_id[20];
           sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
           cout<<Kind_of_HM_struct_id<<endl;
           Diode_CV_struct_id.push_back(Conv_Char_To_String(Kind_of_HM_struct_id));
         }
         if(nlines == 7) {
           cout << line << '\n';
           char Kind_of_HM_set_id[20];
           sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
           cout<<Kind_of_HM_set_id<<endl;
           Diode_CV_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
           cout << line << '\n';
           char str1[20], str2[20];
           sscanf(line.c_str(), "%*s %s %s", str1, str2);
           string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
           cout<<Kind_of_HM_config_id<<endl;
           Diode_CV_config_id.push_back(Kind_of_HM_config_id);
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
           Diode_CV_waiting_time.push_back(Waiting_time_s);
         }
         if(nlines == 11) {
           cout << line << '\n';
           double ac_freq;
           sscanf(line.c_str(), "%*s %lf", &ac_freq);
           cout<<ac_freq<<endl;
           Diode_CV_ac_freq.push_back(ac_freq);
         }
         if(nlines == 12) {
           cout << line << '\n';
           double ac_ampl;
           sscanf(line.c_str(), "%*s %lf", &ac_ampl);
           cout<<ac_ampl<<endl;
           Diode_CV_ac_ampl.push_back(ac_ampl);
        }
        if(nlines == 13) {
          cout << line << '\n';
          double temp_degC;
          sscanf(line.c_str(), "%*s %lf", &temp_degC);
          cout<<temp_degC<<endl;
          Diode_CV_temp.push_back(temp_degC);
          TKelvin = 273 + temp_degC; // [K]
        }
        if(nlines == 14) {
          cout << line << '\n';
          double Av_temp_degC;
          sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
          cout<<Av_temp_degC<<endl;
          Diode_CV_av_temp.push_back(Av_temp_degC);
        }
      }
      else {

      Nintrisic= (5.29*1E+19)*(TMath::Power((TKelvin/300), 2.5))*exp(-6726/TKelvin); // [cm^-3]
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double capctnc_pfrd, rsstnc_mhom;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf %lf", date, time, &voltage, &capctnc_pfrd, &rsstnc_mhom, &temp_degC, &air_temp_degC, &rh_prcnt);
      cout<<date<<","<<time<<","<<voltage<<","<<(float) capctnc_pfrd<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
      //        in >>  Diode_CV_voltage >> Diode_CV_capacitance >> Diode_CV_conductance;
      //  in >>  Diode_CV_voltage >> Diode_CV_capacitance;
    	Diode_CV_vv_d[id][N_meas] =  -voltage;
    	Diode_CV_Cap_d[id][N_meas] = (capctnc_pfrd);

      Diode_CV_temp_meas.push_back(std::vector<double>());
      Diode_CV_temp_meas[id].push_back(temp_degC);
      Diode_CV_air_temp_meas.push_back(std::vector<double>());
      Diode_CV_air_temp_meas[id].push_back(air_temp_degC);
      Diode_CV_rh_prcnt_meas.push_back(std::vector<double>());
      Diode_CV_rh_prcnt_meas[id].push_back(rh_prcnt);
      Diode_CV_timestamp_meas.push_back(std::vector<string>());
      Diode_CV_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      Diode_CV_rsstnc_mhom.push_back(std::vector<double>());
      Diode_CV_rsstnc_mhom[id].push_back(rsstnc_mhom);

    	Diode_CV_vv[id][N_meas] = float(Diode_CV_vv_d[id][N_meas]);
      Diode_CV_Cap[id][N_meas] = float(Diode_CV_Cap_d[id][N_meas]);

      Inv_Diode_CV_Cap[id][N_meas] = 1/(TMath::Power(Diode_CV_Cap[id][N_meas],2));
      Inv_Diode_CV_Cap_d[id][N_meas] = 1/(TMath::Power(Diode_CV_Cap_d[id][N_meas],2));
      //        cout << N_meas << ")  "  <<  Diode_CV_vv[N_meas] << "   " << Diode_CV_cap[N_meas] << endl;
      //cout << N_meas << ")  "  <<  Diode_CV_vv[N_meas] << "   " << Diode_CV_cap[N_meas] << "   " << inv_Diode_CV_cap[N_meas] << endl;
    	N_meas++;
      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   cout  << "----->>>>>>   Number of Measurements = " << N_meas << endl;
   Diode_CV_Nmeas.push_back(N_meas); // insert the number of measurments to an array

   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;


//   TCanvas *cc = new TCanvas("cc","Diode Cap vs Diode_CV_Voltage",500,50,500,500);
   Diode_CV_cv[id] =  new TGraph(N_meas, Diode_CV_vv_d[id], Diode_CV_Cap_d[id]);

   if((id % 2) == 0) Diode_CV_cv[id]->SetTitle(Form("Diode%d_0%s CV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) Diode_CV_cv[id]->SetTitle(Form("Diode%d_0%s CV",HM_Id,Structure_Id[file_id].c_str()));

   Diode_CV_cv[id]->GetXaxis()->SetTitle("Diode_CV_Voltage [V]");
   Diode_CV_cv[id]->GetYaxis()->SetTitle("Diode_CV_Capacitance [pF]");
   Diode_CV_cv[id]->SetDrawOption("AP");
   Diode_CV_cv[id]->SetMarkerStyle(20+id);
   int colorid = 1+id;
   if(colorid == 5) colorid = 46;
   Diode_CV_cv[id]->SetMarkerColor(colorid);
//   Diode_CV_cv[id]->Draw();

   TF1 *f_cv_s = new TF1("f_cv_s","pol0", 300, 500);
   Diode_CV_cv[id]->Fit("f_cv_s","0R+");
   float fix_term_s = Diode_CV_cv[id]->GetFunction("f_cv_s")->GetParameter(0);
   float slop_term_s = 0;

   TF1 *func_A_s = new TF1("func_A_s","[0] + [1]*x",0,500);
   func_A_s->SetParameters(fix_term_s, slop_term_s);
//   func_A_s->Draw("same");

   Diode_CV_Minimum_Capacitance[id] = fix_term_s;
   cout << "----->>> Minimun Diode_CV_Capacitance = " << Diode_CV_Minimum_Capacitance[id]  <<  endl;


// -----------------------------------------------------------------------------------------------------------


//   TCanvas *cc_inv = new TCanvas("cc_inv","Diode 1/Cap2 vs Diode_CV_Voltage",1000,50,500,500);
   Diode_CV_invC2_vs_V[id] =  new TGraph(N_meas, Diode_CV_vv_d[id], Inv_Diode_CV_Cap_d[id]);
   Diode_CV_invC2_vs_V[id]->Sort();
   if((id % 2) == 0) Diode_CV_invC2_vs_V[id]->SetTitle(Form("Diode%d_0%s E,1/C^{2} vs V",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 2) == 1) Diode_CV_invC2_vs_V[id]->SetTitle(Form("Diode%d_0%s W,1/C^{2} vs V",HM_Id,Structure_Id[file_id].c_str()));

   Diode_CV_invC2_vs_V[id]->GetXaxis()->SetTitle("Diode_CV_Voltage [V]");
   Diode_CV_invC2_vs_V[id]->GetYaxis()->SetTitle("Diode_CV_Capacitance [pF]");
   Diode_CV_invC2_vs_V[id]->SetDrawOption("AP");
   Diode_CV_invC2_vs_V[id]->SetMarkerStyle(20+id);
   colorid = 1+id;
   if(colorid == 5) colorid = 46;
   Diode_CV_invC2_vs_V[id]->SetMarkerColor(colorid);
//   Diode_CV_invC2_vs_V[id]->Draw();


   TF1 *f_cv = new TF1("f_cv","pol0", 300, 400);
   Diode_CV_invC2_vs_V[id]->Fit("f_cv","0R+");
   float fix_term = Diode_CV_invC2_vs_V[id]->GetFunction("f_cv")->GetParameter(0);
   float slop_term = 0;

   TF1 *func_A = new TF1("func_A","[0] + [1]*x",0,500);
   func_A->SetParameters(fix_term, slop_term);
//   func_A->Draw("same");

   Diode_CV_CuSpl_InvCV2[id] = new TSpline3("Cubic Spline", Diode_CV_vv_d[id],Inv_Diode_CV_Cap_d[id], N_meas-2, "cb2e2", 0, 0);
//   Diode_CV_CuSpl_InvCV2[id] = new TSpline3("CuSpl_InvCV2", Diode_CV_invC2_vs_V[id]);
//   Diode_CV_CuSpl_InvCV2[id]->SetLineColor(kBlue);
//   Diode_CV_CuSpl_InvCV2[id]->Draw("same");


    // Find the Best Fit point for 1/CV^2
    // --------------------------------
//   for(int i = 0;i<N_meas; i++) {
//       Diode_CV_dev_spline_InvCV2[i] = Diode_CV_CuSpl_InvCV2[id]->Derivative(Diode_CV_vv_d[i]);
//       cout << i << ")  "  <<  Diode_CV_vv_d[i] << "   " << Inv_Diode_CV_Cap_d[i]
//            << "   " << Diode_CV_CuSpl_InvCV2[id]->Eval(Diode_CV_vv_d[i]) << "   " << Diode_CV_dev_spline_InvCV2[i] << endl;
//
//   }

    // get the minimum of the derivative
   int max_index_InvCV2 = 0;
   float max_dev_spline_InvCV2 = Diode_CV_dev_spline_InvCV2[0];
   for (int i = 1;i<N_meas; i++) {
       if(Diode_CV_dev_spline_InvCV2[i] >  max_dev_spline_InvCV2 && i<30) {
         max_dev_spline_InvCV2 = Diode_CV_dev_spline_InvCV2[i];
	        max_index_InvCV2 = i;
       }
   }
   cout  << "---------------- >>>>> Max derivative  InvCV2at i = " << max_index_InvCV2
          << " with InvCV2 value = " << Inv_Diode_CV_Cap[max_index_InvCV2] << endl;

//   TGraph *cv_deriv =  new TGraph(N_meas-2,vv_d,dev_spline_InvCV2);
//   cv_deriv->SetTitle(Form("First Derivative of Diode CV"));
//   cv_deriv->GetXaxis()->SetTitle("Gate Diode_CV_Voltage [V]");
//   cv_deriv->Draw("same");

   int id_fit_low = max_index_InvCV2 - 4;
   int id_fit_high = max_index_InvCV2 + 14;

   float invc2v_low;
   float invc2v_high;
   invc2v_low = Diode_CV_vv_d[id][id_fit_low];
   invc2v_high = Diode_CV_vv_d[id][id_fit_high];

   cout << "----->>> id_fit_low = " << id_fit_low << " id_fit_high = " <<  id_fit_high <<  endl;
   cout << "----->>> invc2v_low = " << invc2v_low << " invc2v_high = " << invc2v_high <<  endl;


   TF1 *f_invc2v = new TF1("f_invc2v","pol1",invc2v_low , invc2v_high);

   Diode_CV_invC2_vs_V[id]->Fit("f_invc2v","0R+");

   double invc2v_fit_fix = Diode_CV_invC2_vs_V[id]->GetFunction("f_invc2v")->GetParameter(0);
   double invc2v_fit_slope = Diode_CV_invC2_vs_V[id]->GetFunction("f_invc2v")->GetParameter(1);

   TF1 *func_B = new TF1("func_B","[0] + [1]*x",100,300);
   func_B->SetParameters(invc2v_fit_fix, invc2v_fit_slope);
//   func_B->Draw("same");

   Diode_CV_FullDepletion_Voltage[id] = (fix_term - invc2v_fit_fix)/invc2v_fit_slope; //
   cout << "----->>> Depletion Diode_CV_Voltage  = " << Diode_CV_FullDepletion_Voltage[id]  <<  endl;
// Calculate hole mobility with Arora model
   float Amin = 54.3;
   float amin = -0.57;
   float Ad = 407;
   float ad = -2.23;
   float AN = 2.35E17;
   float aN = 2.4;
   float Aa = 0.88;
   float aa = -0.146;

//   float N_acceptors =  4E12;
   float N_acceptors =  3E13;
   cout <<"invc2v_fit_slope="<< invc2v_fit_slope << endl;
   Diode_CV_Bulk_Concentration[id] = 2/(q*epsilon_Si*Diode_half_area*Diode_half_area*invc2v_fit_slope);


   //float T_in_Kelvin = Diode_CV_Temp_measurement[id] + 273;

   float mu_min =  Amin*TMath::Power( (TKelvin/300.) , amin);
   float mu_d = Ad*TMath::Power( (TKelvin/300.) , ad);
   float N_0 = AN*TMath::Power( (TKelvin/300.) ,aN );
   float A_star = Aa*TMath::Power( (TKelvin/300.) ,aa );

   float mu_hole = mu_min + mu_d/( 1 + TMath::Power( (N_acceptors/N_0) , A_star) );

   cout << "----->>> Hole Mobility in Diode for this tempetature  = " << mu_hole <<  " cm^2/(V * sec)" << endl;
   cout << "----->>> Standarrd Hole Mobility in Diode  = " << hole_mobility << " cm^2/(V * sec)" <<  endl;

// Calculate Bulk Resistivity taking into account that sensor width: d = 290um = 290 x 10^-4 cm and hole mobility:  mu_h = 450 cm^2/(V * sec)

//   Diode_CV_Bulk_Resistivity[id] = (sensor_width*sensor_width)/(2*epsilon_0*epsilon_Si*hole_mobility*Diode_CV_FullDepletion_Voltage[id]);
   Diode_CV_Bulk_Resistivity[id] = (sensor_width*sensor_width)/(2*epsilon_Si*1E-12*mu_hole*Diode_CV_FullDepletion_Voltage[id]);
   cout << "----->>> Diode Bulk Resistivity  = " << Diode_CV_Bulk_Resistivity[id]/1000.  <<  endl;

   Diode_CV_active_thickness[id] = epsilon_Si*Diode_half_area/Diode_CV_Minimum_Capacitance[id]; // active_thickness=((e0*esi*Area)/(Cmin)) # cm

}
void Diode_CV_PQC3_Analysis_Final(int nFiles)
{

 //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1

 //   gROOT->Reset();
   gStyle->SetOptStat(0);


   //TString Diode_CV_DataFileName[100];
   for(int i=0;i<nFiles;i++) {
    Diode_CV_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/CV/VPX%d_0%s_2-S_HM_E_Left_PQC3_DiodeHalf_CV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    Diode_CV_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/CV/VPX%d_0%s_2-S_HM_W_Left_PQC3_DiodeHalf_CV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {
     cout << endl;
     cout << "Analyze file : " << Diode_CV_DataFileName[i] << endl;
     Read_CV_from_File(Diode_CV_DataFileName[2*i],2*i,i);
     Read_CV_from_File(Diode_CV_DataFileName[2*i+1],2*i+1,i);
  }

   TCanvas *cc_final[2];
   TMultiGraph *mg[2];
   for (int i=0;i<2;i++) {
     cc_final[i] = new TCanvas(Form("cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     mg[i] = new TMultiGraph(Form("mg_%d",i),"Diode  1/C2 vs V");
     if(i == 0) mg[i]->SetTitle("Diode E 1/C^{2} vs V");
     if(i == 1) mg[i]->SetTitle("Diode W 1/C^{2} vs V");

     for(int id=0;id<nFiles;id++) {

        mg[i]->Add(Diode_CV_invC2_vs_V[2*id+i]);

     }
     mg[i]->Draw("LPsame");
     mg[i]->GetXaxis()->SetTitle("Diode_CV_Voltage [V]");
     mg[i]->GetYaxis()->SetTitle("1/C^{2} [pF]^{-2}");
     mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     cc_final[i]->BuildLegend(0.55,0.15,0.85,0.55);


     gPad->Modified();
     gPad->Update();
   }

   cout <<endl;
   cout << "===========================+==================================================================" << endl;
   cout << "---------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  ---------------------------------" << endl;
   cout << "==============================================================================================" << endl;

   cout << "       PQC3  Diode Id                 |  Min Diode Capacitance [pF] |   Diode Depletion Voltage [V]" <<endl;
   cout << "----------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
      cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) << setw(20)  <<  setprecision(3) << Diode_CV_Minimum_Capacitance[2*id] << setw(30)  <<  setprecision(5) <<  Diode_CV_FullDepletion_Voltage[2*id] << endl;
      cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << setw(20)  <<  setprecision(3) << Diode_CV_Minimum_Capacitance[2*id+1] << setw(30)  <<  setprecision(5) <<  Diode_CV_FullDepletion_Voltage[2*id+1] << endl;

   }
   cout << "----------------------------------------------------------------------------------------------" << endl;


}

// Diode_CV XML production (set xml_on to 1)
// --------------------------------------------------------------------------------------------------//
void Diode_CV_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      cout << "Diode CV Xml production is on" << endl;
      for(int j=0; j<=1; j++){
 //       TString xml_filename[nFiles];
        char delimeter1 = '/';
        int pos1 = Diode_CV_DataFileName[2*i+j].Last(delimeter1);
        char delimeter2 = '.';
        int pos2 = Diode_CV_DataFileName[2*i+j].Last(delimeter2);
        xml_filename = Diode_CV_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1))+".xml";
//        cout << Diode_CV_DataFileName[2*i+j](0,pos1) << endl;

        cout << xml_filename << endl;
        ofstream xml_file;

        xml_file.open(Form("./XML_Info_Flute3/VPX%d/",HM_Id)+xml_filename);

        // First create engine
        TXMLEngine Diode_CV_xml;
        // Create main node of document tree
        XMLNodePointer_t ROOT = Diode_CV_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = Diode_CV_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = Diode_CV_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = Diode_CV_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = Diode_CV_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = Diode_CV_xml.NewChild(HEADER, 0, "RUN");
              Diode_CV_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              Diode_CV_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              Diode_CV_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              Diode_CV_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              Diode_CV_xml.NewChild(RUN, 0, "INITIATED_BY_USER", Diode_CV_Operators[2*i+j].c_str());
              Diode_CV_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", Diode_CV_Begin_Timestamp[2*i+j].c_str());
              Diode_CV_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = Diode_CV_xml.NewChild(ROOT, 0, "DATA_SET");
            Diode_CV_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            Diode_CV_xml.NewChild(DATA_SET, 0, "VERSION", "v2");
            XMLNodePointer_t PART = Diode_CV_xml.NewChild(DATA_SET, 0, "PART");
              Diode_CV_xml.NewChild(PART, 0, "NAME_LABEL", Diode_CV_name_labels[2*i+j].c_str());
              Diode_CV_xml.NewChild(PART, 0, "KIND_OF_PART", Diode_CV_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = Diode_CV_xml.NewChild(DATA_SET, 0, "DATA");
              Diode_CV_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC1");
              Diode_CV_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", Diode_CV_struct_id[2*i+j].c_str());
              Diode_CV_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", Diode_CV_set_id[2*i+j].c_str());
              Diode_CV_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", Diode_CV_config_id[2*i+j].c_str());
              Diode_CV_xml.NewChild(DATA, 0, "EQUIPMENT", "HP 4192A");
              Diode_CV_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(Diode_CV_waiting_time[2*i+j],3).c_str());
              Diode_CV_xml.NewChild(DATA, 0, "AC_FREQ_HZ", Conv_float_to_string_scientific_format(Diode_CV_ac_freq[2*i+j],2).c_str());
              Diode_CV_xml.NewChild(DATA, 0, "AC_AMPL_V", Conv_float_to_string(Diode_CV_ac_ampl[2*i+j],3).c_str());
              Diode_CV_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(Diode_CV_temp[2*i+j],3).c_str());
              Diode_CV_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(Diode_CV_av_temp[2*i+j],3).c_str());
              Diode_CV_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "Diode_CV1");
//              Diode_CV_xml.NewChild(DATA, 0, "FILE_NAME", xml_filename[2*i+j]);
              Diode_CV_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = Diode_CV_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = Diode_CV_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = Diode_CV_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  Diode_CV_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_CV");
                  Diode_CV_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon CV Test");
              XMLNodePointer_t DATA_SET_CDS1 = Diode_CV_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                Diode_CV_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                Diode_CV_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "V2");
                XMLNodePointer_t PART_CDS1 = Diode_CV_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  Diode_CV_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", Diode_CV_name_labels[2*i+j].c_str());
                  Diode_CV_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", Diode_CV_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<Diode_CV_Nmeas[2*i+j]-1; k++){
                XMLNodePointer_t DATA_CDS1 = Diode_CV_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  Diode_CV_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(Diode_CV_vv[2*i+j][k],3).c_str());
                  Diode_CV_xml.NewChild(DATA_CDS1, 0, "CAPCTNC_PFRD", Conv_float_to_string(Diode_CV_Cap[2*i+j][k],3).c_str());
                  Diode_CV_xml.NewChild(DATA_CDS1, 0, "RESSTNC_MOHM", Conv_float_to_string(Diode_CV_rsstnc_mhom[2*i+j][k],3).c_str());
                  Diode_CV_xml.NewChild(DATA_CDS1, 0, "TIME", Diode_CV_timestamp_meas[2*i+j][k].c_str());
                  Diode_CV_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(Diode_CV_temp_meas[2*i+j][k],3).c_str());
                  Diode_CV_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(Diode_CV_air_temp_meas[2*i+j][k],3).c_str());
                  Diode_CV_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(Diode_CV_rh_prcnt_meas[2*i+j][k],3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = Diode_CV_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = Diode_CV_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = Diode_CV_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      Diode_CV_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_CV_PAR");
                      Diode_CV_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon CV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = Diode_CV_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    Diode_CV_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    Diode_CV_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "CV_measurement-004");
                    XMLNodePointer_t PART_CDS2 = Diode_CV_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      Diode_CV_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", Diode_CV_name_labels[2*i+j].c_str());
                      Diode_CV_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", Diode_CV_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = Diode_CV_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      Diode_CV_xml.NewChild(DATA_CDS2, 0, "VD_V", Conv_float_to_string(Diode_CV_FullDepletion_Voltage[2*i+j],3).c_str());
                      Diode_CV_xml.NewChild(DATA_CDS2, 0, "CMIN_PFRD", Conv_float_to_string(Diode_CV_Minimum_Capacitance[2*i+j],3).c_str());
                      Diode_CV_xml.NewChild(DATA_CDS2, 0, "RHO_KOHMCM", Conv_float_to_string(Diode_CV_Bulk_Resistivity[2*i+j]/1000,3).c_str());
                      Diode_CV_xml.NewChild(DATA_CDS2, 0, "NA", Conv_float_to_string_scientific_format(Diode_CV_Bulk_Concentration[2*i+j],3).c_str());
                      Diode_CV_xml.NewChild(DATA_CDS2, 0, "D_UM", Conv_float_to_string_scientific_format(Diode_CV_active_thickness[2*i+j]*1e+4,3).c_str());

        XMLDocPointer_t Diode_CV_xmldoc = Diode_CV_xml.NewDoc();
        Diode_CV_xml.DocSetRootElement(Diode_CV_xmldoc, ROOT);

        // Save document to file
        Diode_CV_xml.SaveDoc(Diode_CV_xmldoc, Form("./XML_Info_Flute3/VPX%d/",HM_Id)+xml_filename);

        // Release memory before exit
        Diode_CV_xml.FreeDoc(Diode_CV_xmldoc);
      }
    }
  }
}
// ----

void Flute3_2_S_Quick_Characterization_with_xml(int nFiles)
{
   VDPBulk_PQC_Flute3_Analysis_Final(nFiles);
   VDPBulk_PQC_xml_production(xml_on, nFiles);

   Diode_IV_PQC3_Analysis_Final(nFiles);
   Diode_IV_PQC_xml_production(xml_on, nFiles);

   Diode_CV_PQC3_Analysis_Final(nFiles);
   Diode_CV_PQC_xml_production(xml_on, nFiles);

   CVS_Output_File.open(Form("CSV_Info/VPX%d/Flute3_VPX%d_0xx_2-S_HM_E_W_Left_Quick_Info.csv",HM_Id,HM_Id));


   cout << endl << endl;
   cout << "====================================================================================================" << endl;
   cout << "-------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E_Left",HM_Id) << "  ---------------------------------" << endl;
   cout << "===================================================================================================================" << endl;
   cout << "     Flute1  Id            |   VDPBulk   |    VDPBulk_Res    | Diode I600 | Diode Cap | Diode Vdep | Diode_Rho_Bulk" <<endl; 
   cout << "                           |   [kOhm/sq] |     [kOhmcm]      |    [nA]    |    [pf]   |     [V]    |   [kOhm cm]    " <<endl; 
   cout << "                           | stad | rot  | F=1.089 | F=1.218 |            |           |            |  mu_h = arora    " <<endl ; 
   cout << "     Spec Limit            |             | 3.5-8.0 | 3.5-8.0 |    <2.5    |           |   < 350    |                " <<endl ; 

   cout << "----------------------------------------------------------------------------------------------------------------------" << endl;

   cout << showpoint;  
   for(int id=0;id<nFiles;id++) {  
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())  
	       <<  setw(7) <<  setprecision(4) <<  VDPBulk_s_Rsh_E[id]/1000.  <<  setw(7)  << VDPBulk_r_Rsh_E[id]/1000. 
	       <<  setw(8) <<  setprecision(4) << VDPBulk_Resistivity_s_F1_E[id]/1000.  <<  setw(10)  << VDPBulk_Resistivity_s_F2_E[id]/1000.
	       <<  setw(11)  <<  setprecision(2) <<  Diode_Current_at_600V_over_Volume[2*id]
	       <<  setw(13)  <<  setprecision(3) << Diode_CV_Minimum_Capacitance[2*id] << setw(14)  <<  setprecision(5) <<  Diode_CV_FullDepletion_Voltage[2*id]
	       <<  setw(14)  <<  setprecision(5) <<  Diode_CV_Bulk_Resistivity[2*id]/1000.<< endl; 
	          
   
         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())  
	       <<  setw(7) <<  setprecision(4) <<  VDPBulk_s_Rsh_W[id]/1000.  <<  setw(7)  << VDPBulk_r_Rsh_W[id]/1000. 
	       <<  setw(8) <<  setprecision(4) << VDPBulk_Resistivity_s_F1_W[id]/1000.  <<  setw(10)  << VDPBulk_Resistivity_s_F2_W[id]/1000.
	       <<  setw(11)  <<  setprecision(2) <<  Diode_Current_at_600V_over_Volume[2*id+1]  
	       <<  setw(13)  <<  setprecision(3) << Diode_CV_Minimum_Capacitance[2*id+1] << setw(14)  <<  setprecision(5) <<  Diode_CV_FullDepletion_Voltage[2*id+1]
	       << setw(14)  <<  setprecision(5) <<  Diode_CV_Bulk_Resistivity[2*id+1]/1000.  << endl; 
   }   


   cout << "----------------------------------------------------------------------------------------------------------------------" << endl;
   CVS_Output_File << "     Flute3  Id            ;Mode ;   VDPBulk  ;  ; VDPBulk_Res;     ; Diode I600  ; Diode Cap    ; Diode Vdep; Diode Rho_Bulk " << endl; 
   CVS_Output_File << "                           ;     ;  [kOhm/sq] ;  ;  [kOhmcm] ;      ;    [nA]     ;  [pF]        ;    [V]    ;   [kOhm cm]    "<< endl; 
   CVS_Output_File << "                           ;     ;  stad  ;  rot ; F=1.089 ; F=1.218;             ;              ;           ;   mu_h = arora   " << endl; 
   CVS_Output_File << "     Spec Limit            ;     ;        ;      ; 3.5-8.0 ; 3.5-8.0;     <2.5    ;             ; < 350      ;                " << endl; 


   for(int id=0;id<nFiles;id++) {  
        CVS_Output_File   <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) << ";" << "Quick" << ";"  
	                  << VDPBulk_s_Rsh_E[id]/1000.  << ";"   << VDPBulk_r_Rsh_E[id]/1000.<< ";"  
	                  << VDPBulk_Resistivity_s_F1_E[id]/1000.  <<  ";"  << VDPBulk_Resistivity_s_F2_E[id]/1000.<< ";" 
	                  << Diode_Current_at_600V_over_Volume[2*id] <<  ";" << Diode_CV_Minimum_Capacitance[2*id] << ";"
			  <<  Diode_CV_FullDepletion_Voltage[2*id] << ";" << Diode_CV_Bulk_Resistivity[2*id]/1000.
		          << endl; 
	          
   
         CVS_Output_File  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << ";" << "Quick" << ";"	 
	                  << VDPBulk_s_Rsh_W[id]/1000.  << ";"   << VDPBulk_r_Rsh_W[id]/1000.<< ";"  
	                  << VDPBulk_Resistivity_s_F1_W[id]/1000.  <<  ";"  << VDPBulk_Resistivity_s_F2_W[id]/1000.<< ";" 
	                  << Diode_Current_at_600V_over_Volume[2*id+1] <<  ";" << Diode_CV_Minimum_Capacitance[2*id+1] << ";"
			  <<  Diode_CV_FullDepletion_Voltage[2*id+1] << ";" << Diode_CV_Bulk_Resistivity[2*id+1]/1000.
		          << endl; 
 
   } 
}
