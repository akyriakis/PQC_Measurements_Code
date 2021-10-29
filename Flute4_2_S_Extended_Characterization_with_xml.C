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
//#include <boost/algorithm/string.hpp>

double epsilon_0;
double epsilon_Si;
double epsilon_SiO2;
double GCD_area;
double q;
double kB;
double TKelvin;
double Nintrisic;


float VDP_CCStructures_current_d[500][500];
float VDP_CCStructures_voltage_d[500][500];
float VDP_CCStructures_current[500][500];
float VDP_CCStructures_voltage[500][500];
float VDP_CCStructures_Temp_measurement[600];
float VDP_CCStructures_Humid_measurement[600];
double VDP_CCStructures_A_constant[600];
double VDP_CCStructures_A_slope[600];
float VDP_CCStructures_cv_low[600];
float VDP_CCStructures_cv_high[600];
double VDP_CCEdge_Rc_E[100];
double VDP_CCPoly_Rc_E[100];
double VDP_CCStrip_Rc_E[100];
double VDP_CCEdge_Rc_W[100];
double VDP_CCPoly_Rc_W[100];
double VDP_CCStrip_Rc_W[100];
double VDP_CC_Rc_all[100];
TGraph *VDP_CCStructures_iv[600];
vector<vector <string> > VDP_CCStructures_timestamp_meas;
vector<vector <double> > VDP_CCStructures_temp_meas;
vector<vector <double> > VDP_CCStructures_air_temp_meas;
vector<vector <double> > VDP_CCStructures_rh_prcnt_meas;

vector<std::string> VDP_CCStructures_Operators;
vector<std::string> VDP_CCStructures_Begin_Timestamp;
vector<std::string> VDP_CCStructures_name_labels;
vector<std::string> VDP_CCStructures_kind_of_parts;
vector<std::string> VDP_CCStructures_kind_of_HM_flute_id;
vector<std::string> VDP_CCStructures_struct_id;
vector<std::string> VDP_CCStructures_set_id;
vector<std::string> VDP_CCStructures_config_id;
vector<std::string> VDP_CCStructures_equipment;
vector<double> VDP_CCStructures_waiting_time;
vector<double> VDP_CCStructures_temp;
vector<double> VDP_CCStructures_av_temp;
vector<double> VDP_CCStructures_Nmeas;


// --------------------------------------------------------

float VDPoly_current_d[500];
float VDPoly_voltage_d[500];
float VDPoly_Temp_measurement[600];
float VDPoly_Humid_measurement[600];
double VDPoly_A_constant[600];
double VDPoly_A_slope[600];
float VDPoly_cv_low[600];
float VDPoly_cv_high[600];
double VDPoly_s_Rsh_E[100];
double VDPoly_r_Rsh_E[100];
double VDPoly_s_Rsh_W[100];
double VDPoly_r_Rsh_W[100];
TGraph *VDPoly_iv[600];


float CBKRPoly_current_d[500][500];
float CBKRPoly_voltage_d[500][500];
float CBKRPoly_current[500][500];
float CBKRPoly_voltage[500][500];
float CBKRPoly_Temp_measurement[600];
float CBKRPoly_Humid_measurement[600];
double CBKRPoly_A_constant[600];
double CBKRPoly_A_slope[600];
float CBKRPoly_cv_low[600];
float CBKRPoly_cv_high[600];
double CBKRPoly_s_Rsh_E[100];
double CBKRPoly_r_Rsh_E[100];
double CBKRPoly_s_Rsh_W[100];
double CBKRPoly_r_Rsh_W[100];
double CBKRPoly_RC_all[100];
double CBKRPoly_R_all[100];
double CBKRPoly_s_Rsh_Rgeom_Corrected_E[100];
double CBKRPoly_r_Rsh_Rgeom_Corrected_E[100];
double CBKRPoly_s_Rsh_Rgeom_Corrected_W[100];
double CBKRPoly_r_Rsh_Rgeom_Corrected_W[100];
TGraph *CBKRPoly_iv[600];
vector<vector <string> > CBKRPoly_timestamp_meas;
vector<vector <double> > CBKRPoly_temp_meas;
vector<vector <double> > CBKRPoly_air_temp_meas;
vector<vector <double> > CBKRPoly_rh_prcnt_meas;
vector<std::string> CBKRPoly_Operators;
vector<std::string> CBKRPoly_Begin_Timestamp;
vector<std::string> CBKRPoly_name_labels;
vector<std::string> CBKRPoly_kind_of_parts;
vector<std::string> CBKRPoly_kind_of_HM_flute_id;
vector<std::string> CBKRPoly_struct_id;
vector<std::string> CBKRPoly_set_id;
vector<std::string> CBKRPoly_config_id;
vector<std::string> CBKRPoly_equipment;
vector<double> CBKRPoly_waiting_time;
vector<double> CBKRPoly_temp;
vector<double> CBKRPoly_av_temp;
vector<double> CBKRPoly_Nmeas;

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


float CBKRStrip_voltage_d[500][500];
float CBKRStrip_current_d[500][500];
float CBKRStrip_voltage[500][500];
float CBKRStrip_current[500][500];
float CBKRStrip_Temp_measurement[600];
float CBKRStrip_Humid_measurement[600];
double CBKRStrip_A_constant[600];
double CBKRStrip_A_slope[600];
float CBKRStrip_cv_low[600];
float CBKRStrip_cv_high[600];
double CBKRStrip_s_Rsh_E[100];
double CBKRStrip_r_Rsh_E[100];
double CBKRStrip_s_Rsh_W[100];
double CBKRStrip_r_Rsh_W[100];
double CBKRStrip_RC_all[100];
double CBKRStrip_R_all[100];
double CBKRStrip_s_Rsh_Rgeom_Corrected_E[100];
double CBKRStrip_r_Rsh_Rgeom_Corrected_E[100];
double CBKRStrip_s_Rsh_Rgeom_Corrected_W[100];
double CBKRStrip_r_Rsh_Rgeom_Corrected_W[100];
TGraph *CBKRStrip_iv[600];
vector<vector <string> > CBKRStrip_timestamp_meas;
vector<vector <double> > CBKRStrip_temp_meas;
vector<vector <double> > CBKRStrip_air_temp_meas;
vector<vector <double> > CBKRStrip_rh_prcnt_meas;
vector<std::string> CBKRStrip_Operators;
vector<std::string> CBKRStrip_Begin_Timestamp;
vector<std::string> CBKRStrip_name_labels;
vector<std::string> CBKRStrip_kind_of_parts;
vector<std::string> CBKRStrip_kind_of_HM_flute_id;
vector<std::string> CBKRStrip_struct_id;
vector<std::string> CBKRStrip_set_id;
vector<std::string> CBKRStrip_config_id;
vector<std::string> CBKRStrip_equipment;
vector<double> CBKRStrip_waiting_time;
vector<double> CBKRStrip_temp;
vector<double> CBKRStrip_av_temp;
vector<double> CBKRStrip_Nmeas;
// --------------------------------------------------------

float GCD05_Temp_measurement[100];
float GCD05_Humid_measurement[100];
float GCD05_Vbias_measurement[100];
double GCD05_A_constant[100];
double GCD05_A_slope[100];
double GCD05_D_constant[100];
double GCD05_D_slope[100];
double GCD05_I_constant[100];
double GCD05_I_slope[100];
float GCD05_iv_low[100];
float GCD05_iv_high[100];
float GCD05_iv1_low[100];
float GCD05_iv1_high[100];
float GCD05_iv2_low[100];
float GCD05_iv2_high[100];
float GCD05_I_accumulation[100];
float GCD05_I_depletion[100];
float GCD05_I_inversion[100];
float GCD05_VFB_Derivative[100];
float GCD05_I_surface[100];
float GCD05_S_interface_recombination_Velocity[100];
float GCD05_D_it[100];
float GCD05_N_it[100];
float GCD05_area;
double GCD05_vv[500][500];
double GCD05_current[500][500];
double GCD05_vv_d[500][500];
double GCD05_current_d[500][500];
double GCD05_vv_d_reduced[500];
double GCD05_current_d_reduced[500];
double GCD05_dev_spline_simple[500];
TCanvas *GCD05_cc_spline[100];
TCanvas *GCD05_cc_spline_der[100];
TGraph *GCD05_iv[100];
TSpline3 * GCD05_CuSpl_iv[100];
vector<vector <string> > GCD05_timestamp_meas;
vector<vector <double> > GCD05_temp_meas;
vector<vector <double> > GCD05_air_temp_meas;
vector<vector <double> > GCD05_rh_prcnt_meas;
vector<std::string> GCD05_Operators;
vector<std::string> GCD05_Begin_Timestamp;
vector<std::string> GCD05_name_labels;
vector<std::string> GCD05_kind_of_parts;
vector<std::string> GCD05_kind_of_HM_flute_id;
vector<std::string> GCD05_struct_id;
vector<TString> GCD05_set_id;
vector<std::string> GCD05_config_id;
vector<std::string> GCD05_equipment;
vector<double> GCD05_waiting_time;
vector<double> GCD05_Bias;
vector<double> GCD05_temp;
vector<double> GCD05_av_temp;
vector<double> GCD05_compliance;
vector<double> GCD05_Nmeas;


ROOT::Math::Interpolator * CSpline_Interpolator[100];

ofstream CVS_Output_File;

TString Data_Dir = "/home/akyriakis/MOS-Measurements_Analysis/PQC_Measurements/Data_New_Format/";

TString CBKRPoly_DataFileName[600];
TString CBKRStrip_DataFileName[600];
TString LineWidthStrip_DataFileName[600];
TString VDP_DataFileName[100];
TString VDP_CCStructures_DataFileName[600];
TString GCD05_DataFileName[100];


int xml_on = 1;

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
//string Structure_Id[40] = {"06","19","49"};
//int HM_Id = 37401;
//string Structure_Id[40] = {"09","26","41"}; 
//int HM_Id = 37408;
//string Structure_Id[40] = {"10","21","46"}; 
int HM_Id = 37897;
string Structure_Id[40] = {"05","18","43"}; 

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
double VDP_CCStructures_Read_Info_from_File(TString VDP_CCStructures_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(VDP_CCStructures_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << VDP_CCStructures_DataFileName << " with Id = " << id << endl;

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
         VDP_CCStructures_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         VDP_CCStructures_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         VDP_CCStructures_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         VDP_CCStructures_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         VDP_CCStructures_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         //std::string data = Conv_Char_To_String(Kind_of_HM_struct_id);
         //boost::to_upper(data);
         // convert string to upper case
         char uppercase[100];
         int i=0;
         for (int i=0; i<=20; i++){
           uppercase[i]=toupper(Kind_of_HM_struct_id[i]);
         }
         VDP_CCStructures_struct_id.push_back(Conv_Char_To_String(uppercase));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         VDP_CCStructures_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         VDP_CCStructures_config_id.push_back(Kind_of_HM_config_id);
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
         VDP_CCStructures_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         VDP_CCStructures_waiting_time.push_back(Waiting_time_s);
         }
        if(nlines == 12) {
        cout << line << '\n';
        double temp_degC;
        sscanf(line.c_str(), "%*s %lf", &temp_degC);
        cout<<temp_degC<<endl;
        VDP_CCStructures_temp.push_back(temp_degC);
       }
        if(nlines == 13) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        VDP_CCStructures_av_temp.push_back(Av_temp_degC);
       }
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double current_namp;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf ", date, time, &current_namp, &voltage, &temp_degC, &air_temp_degC, &rh_prcnt);
      //cout<<date<<","<<time<<","<<voltage<<","<<(float) current_namp<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
   //        in >>  VDP_CCStructures_voltage >> VDP_CCStructures_capacitance >> VDP_CCStructures_conductance;
      //  in >>  VDP_CCStructures_voltage >> VDP_CCStructures_capacitance;
      VDP_CCStructures_voltage_d[id][N_meas] = voltage*1E-9;
      VDP_CCStructures_current_d[id][N_meas] = current_namp; // this is not current in namp
      VDP_CCStructures_temp_meas.push_back(std::vector<double>());
      VDP_CCStructures_temp_meas[id].push_back(temp_degC);
      VDP_CCStructures_air_temp_meas.push_back(std::vector<double>());
      VDP_CCStructures_air_temp_meas[id].push_back(air_temp_degC);
      VDP_CCStructures_rh_prcnt_meas.push_back(std::vector<double>());
      VDP_CCStructures_rh_prcnt_meas[id].push_back(rh_prcnt);
      VDP_CCStructures_timestamp_meas.push_back(std::vector<string>());
      VDP_CCStructures_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      VDP_CCStructures_voltage[id][N_meas] = float(VDP_CCStructures_voltage_d[id][N_meas]);
      VDP_CCStructures_current[id][N_meas] = float(VDP_CCStructures_current_d[id][N_meas]);

      VDP_CCStructures_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));


    //        cout << N_meas << ")  "  <<  VDP_CCStructures_vv[N_meas] << "   " << VDP_CCStructures_cap[N_meas] << endl;
      cout << N_meas << ")  "  <<  VDP_CCStructures_voltage[id][N_meas] << "   " << VDP_CCStructures_current[id][N_meas]  << endl;
      N_meas++;

      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   VDP_CCStructures_Nmeas.push_back(N_meas);

   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl;

   TGaxis::SetMaxDigits(3);

//   TCanvas *cc_spline = new TCanvas();
   VDP_CCStructures_iv[id] =  new TGraph(N_meas,VDP_CCStructures_current_d[id], VDP_CCStructures_voltage_d[id]);
   if((id % 6) == 0) VDP_CCStructures_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC4_VDP_CCEdge_IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 6) == 1) VDP_CCStructures_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC4_VDP_CCEdge_IV",HM_Id,Structure_Id[file_id].c_str()));

   if((id % 6) == 2) VDP_CCStructures_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC4_VDP_CCCCPoly_IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 6) == 3) VDP_CCStructures_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC4_VDP_CCCCPoly_IV",HM_Id,Structure_Id[file_id].c_str()));

   if((id % 6) == 4) VDP_CCStructures_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC4_VDP_CCStrip_IV",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 6) == 5) VDP_CCStructures_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC4_VDP_CCStrip_IV",HM_Id,Structure_Id[file_id].c_str()));
   VDP_CCStructures_iv[id]->GetXaxis()->SetTitle("VDP_CCStructures Current [A]");
   VDP_CCStructures_iv[id]->GetYaxis()->SetTitle("VDP_CCStructures Voltage [V]");
   VDP_CCStructures_iv[id]->SetDrawOption("AP");
   VDP_CCStructures_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   VDP_CCStructures_iv[id]->SetMarkerColor(colorid);
//   VDP_CCStructures_iv[id]->Draw();


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   VDP_CCStructures_cv_low[id] = VDP_CCStructures_current_d[id][3];
   VDP_CCStructures_cv_high[id] = VDP_CCStructures_current_d[id][38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", VDP_CCStructures_cv_low[id], VDP_CCStructures_cv_high[id]);
//   VDP_CCStructures_iv[id]->Fit("f_cv","0R+");
   VDP_CCStructures_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   VDP_CCStructures_A_constant[id] = VDP_CCStructures_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   VDP_CCStructures_A_slope[id]  = VDP_CCStructures_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double R_sh = VDP_CCStructures_A_slope[id]; // VDP_CCStructures R_sh
   double R_sh_error = (VDP_CCStructures_iv[id]->GetFunction("f_cv")->GetParError(1));

   cout << "----->>>   R_sh Voltage = " << R_sh << " +/- " <<  R_sh_error<< " [Ohm/sq]" << endl;

   return R_sh;

}
void VDP_CCStructures_PQC_Flute4_Analysis_Final(int nFiles)
{

 //  AC_ampl(V)=0.250        AC_freq(Hz)=10.000E+3   Settling_time(s)=1000.0 , Temperature(C)=22.1

 //   gROOT->Reset();
   gStyle->SetOptStat(0);


   //TString VDP_CCStructures_DataFileName[600];
   for(int i=0;i<nFiles;i++) {
    VDP_CCStructures_DataFileName[6*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC4_CCEdge_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_CCStructures_DataFileName[6*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC4_CCEdge_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDP_CCStructures_DataFileName[6*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC4_CCPoly_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_CCStructures_DataFileName[6*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC4_CCPoly_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());


    VDP_CCStructures_DataFileName[6*i+4] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC4_CCStrip_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDP_CCStructures_DataFileName[6*i+5] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC4_CCStrip_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {

     cout << endl << endl << "CCEdge Structure E: Analyze file : " << VDP_CCStructures_DataFileName[6*i] << endl;
     VDP_CCEdge_Rc_E[i] = VDP_CCStructures_Read_Info_from_File(VDP_CCStructures_DataFileName[6*i],(6*i),i);

     cout << endl << endl << "CCEdge Structure W: Analyze file : " << VDP_CCStructures_DataFileName[6*i+1] << endl;
     VDP_CCEdge_Rc_W[i] = VDP_CCStructures_Read_Info_from_File(VDP_CCStructures_DataFileName[6*i+1],(6*i+1),i);


     cout << endl << endl << "CCPoly Structure W: Analyze file : " << VDP_CCStructures_DataFileName[6*i+2] << endl;
     VDP_CCPoly_Rc_E[i] = VDP_CCStructures_Read_Info_from_File(VDP_CCStructures_DataFileName[6*i+2],(6*i+2),i);

     cout << endl << endl << "CCPoly Structure W: Analyze file : " << VDP_CCStructures_DataFileName[6*i+3] << endl;
     VDP_CCPoly_Rc_W[i] = VDP_CCStructures_Read_Info_from_File(VDP_CCStructures_DataFileName[6*i+3],(6*i+3),i);


     cout <<endl << endl << "CCStrip Structure W: Analyze file : " << VDP_CCStructures_DataFileName[6*i+4] << endl;
     VDP_CCStrip_Rc_E[i] = VDP_CCStructures_Read_Info_from_File(VDP_CCStructures_DataFileName[6*i+4],(6*i+4),i);

     cout << endl << endl << "CCStrip Structure W: Analyze file : " << VDP_CCStructures_DataFileName[6*i+5] << endl;
     VDP_CCStrip_Rc_W[i] = VDP_CCStructures_Read_Info_from_File(VDP_CCStructures_DataFileName[6*i+5],(6*i+5),i);

     VDP_CC_Rc_all[6*i+0] = VDP_CCEdge_Rc_E[i];
     VDP_CC_Rc_all[6*i+1] = VDP_CCEdge_Rc_W[i];
     VDP_CC_Rc_all[6*i+2] = VDP_CCPoly_Rc_E[i];
     VDP_CC_Rc_all[6*i+3] = VDP_CCPoly_Rc_W[i];
     VDP_CC_Rc_all[6*i+4] = VDP_CCStrip_Rc_E[i];
     VDP_CC_Rc_all[6*i+5] = VDP_CCStrip_Rc_W[i];

   }

   TCanvas *VDP_CC_cc_final[6];
   TMultiGraph *VDP_CC_mg[6];
   for (int i=0;i<6;i++) {
     VDP_CC_cc_final[i] = new TCanvas(Form("VDP_CC_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     VDP_CC_mg[i] = new TMultiGraph(Form("VDP_CC_mg_%d",i),"VDP_CCStructures R_sh");
     if(i == 0) VDP_CC_mg[i]->SetTitle("VDP_CCEdge_E R_sh ");
     if(i == 1) VDP_CC_mg[i]->SetTitle("VDP_CCEdge_W R_sh ");

     if(i == 2) VDP_CC_mg[i]->SetTitle("VDP_CCPoly_Rc_E R_sh ");
     if(i == 3) VDP_CC_mg[i]->SetTitle("VDP_CCPoly_Rc_W R_sh ");


     if(i == 4) VDP_CC_mg[i]->SetTitle("VDP_CCStrip_Rc_E R_sh ");
     if(i == 5) VDP_CC_mg[i]->SetTitle("VDP_CCStrip_Rc_W R_sh ");

    // VDP_CCStructures_iv[4*i+]->Draw("ALP");
     for(int id=0;id<nFiles;id++) {

        VDP_CC_mg[i]->Add(VDP_CCStructures_iv[6*id+i]);

     }
     VDP_CC_mg[i]->Draw("LPsame");
     VDP_CC_mg[i]->GetXaxis()->SetTitle("Current [A]");
     VDP_CC_mg[i]->GetYaxis()->SetTitle("Voltage [V]");
     VDP_CC_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     VDP_CC_cc_final[i]->BuildLegend(0.15,0.5,0.5,0.85);


     gPad->Modified();
     gPad->Update();
   }


   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;

   cout << "      VDP_CCStructures   Id   |  VDP_CCEdge_Rc [kOhm/sq]   | VDP_CCPoly_Rc [MOhm/sq] | VDP_CCStrip_Rc [KOhm/sq] " <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(20) <<  setprecision(4) << VDP_CCEdge_Rc_E [id]*1E-3  <<  setw(30) << VDP_CCPoly_Rc_E[id]*1E-6   << setw(30) << VDP_CCStrip_Rc_E[id]*1E-3   << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(20) <<  setprecision(4) << VDP_CCEdge_Rc_W [id]*1E-3   <<  setw(30) << VDP_CCPoly_Rc_W[id]*1E-6 << setw(30) << VDP_CCStrip_Rc_W[id]*1E-3  << endl;

   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
}

void VDP_CCStructures_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  cout << "VDP Xml production is on" << endl;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<6; j++){
          //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = VDP_CCStructures_DataFileName[6*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = VDP_CCStructures_DataFileName[6*i+j].Last(delimeter2);
         xml_filename = VDP_CCStructures_DataFileName[6*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[4*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute4/VPX%d/",HM_Id) +xml_filename);

        // First create engine
        TXMLEngine VDP_CCStructures_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = VDP_CCStructures_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = VDP_CCStructures_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = VDP_CCStructures_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = VDP_CCStructures_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = VDP_CCStructures_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = VDP_CCStructures_xml.NewChild(HEADER, 0, "RUN");
              VDP_CCStructures_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              VDP_CCStructures_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              VDP_CCStructures_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              VDP_CCStructures_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              VDP_CCStructures_xml.NewChild(RUN, 0, "INITIATED_BY_USER", VDP_CCStructures_Operators[6*i+j].c_str());
              VDP_CCStructures_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", VDP_CCStructures_Begin_Timestamp[6*i+j].c_str());
              VDP_CCStructures_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = VDP_CCStructures_xml.NewChild(ROOT, 0, "DATA_SET");
            VDP_CCStructures_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            VDP_CCStructures_xml.NewChild(DATA_SET, 0, "VERSION", "IV_measurement");
            XMLNodePointer_t PART = VDP_CCStructures_xml.NewChild(DATA_SET, 0, "PART");
              VDP_CCStructures_xml.NewChild(PART, 0, "NAME_LABEL", VDP_CCStructures_name_labels[6*i+j].c_str());
              VDP_CCStructures_xml.NewChild(PART, 0, "KIND_OF_PART", VDP_CCStructures_kind_of_parts[6*i+j].c_str());
            XMLNodePointer_t DATA = VDP_CCStructures_xml.NewChild(DATA_SET, 0, "DATA");
              VDP_CCStructures_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC4");
              VDP_CCStructures_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", VDP_CCStructures_struct_id[6*i+j].c_str());
              VDP_CCStructures_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", VDP_CCStructures_set_id[6*i+j].c_str());
              VDP_CCStructures_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", VDP_CCStructures_config_id[6*i+j].c_str());
              VDP_CCStructures_xml.NewChild(DATA, 0, "EQUIPMENT", VDP_CCStructures_equipment[6*i+j].c_str());
              VDP_CCStructures_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(VDP_CCStructures_waiting_time[6*i+j], 5).c_str());
              VDP_CCStructures_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(VDP_CCStructures_temp[6*i+j], 3).c_str());
              VDP_CCStructures_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(VDP_CCStructures_av_temp[6*i+j], 3).c_str());
              VDP_CCStructures_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", VDP_CCStructures_struct_id[6*i+j].c_str());
              VDP_CCStructures_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = VDP_CCStructures_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = VDP_CCStructures_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = VDP_CCStructures_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  VDP_CCStructures_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  VDP_CCStructures_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = VDP_CCStructures_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                VDP_CCStructures_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                VDP_CCStructures_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "IV_measurement");
                XMLNodePointer_t PART_CDS1 = VDP_CCStructures_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  VDP_CCStructures_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", VDP_CCStructures_name_labels[6*i+j].c_str());
                  VDP_CCStructures_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", VDP_CCStructures_kind_of_parts[6*i+j].c_str());
                for (int k=0; k<VDP_CCStructures_Nmeas[6*i+j]-1; k++){
                    XMLNodePointer_t DATA_CDS1 = VDP_CCStructures_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                    cout<< Conv_float_to_string(VDP_CCStructures_current[6*i+j][k], 4).c_str() << endl;
                      VDP_CCStructures_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string_scientific_format(VDP_CCStructures_voltage[6*i+j][k], 3).c_str());
                      VDP_CCStructures_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(VDP_CCStructures_current[6*i+j][k]*1E+9, 3).c_str());
                      VDP_CCStructures_xml.NewChild(DATA_CDS1, 0, "TIME", VDP_CCStructures_timestamp_meas[6*i+j][k].c_str());
                      VDP_CCStructures_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(VDP_CCStructures_temp_meas[6*i+j][k], 3).c_str());
                      VDP_CCStructures_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(VDP_CCStructures_air_temp_meas[6*i+j][k], 3).c_str());
                      VDP_CCStructures_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(VDP_CCStructures_rh_prcnt_meas[6*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = VDP_CCStructures_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = VDP_CCStructures_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = VDP_CCStructures_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      VDP_CCStructures_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      VDP_CCStructures_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = VDP_CCStructures_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    VDP_CCStructures_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    VDP_CCStructures_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement");
                    XMLNodePointer_t PART_CDS2 = VDP_CCStructures_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      VDP_CCStructures_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", VDP_CCStructures_name_labels[6*i+j].c_str());
                      VDP_CCStructures_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", VDP_CCStructures_kind_of_parts[6*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = VDP_CCStructures_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      VDP_CCStructures_xml.NewChild(DATA_CDS2, 0, "R_OHM", Conv_float_to_string_scientific_format(VDP_CC_Rc_all[6*i+j], 3).c_str());

        XMLDocPointer_t VDP_CCStructures_xmldoc = VDP_CCStructures_xml.NewDoc();
        VDP_CCStructures_xml.DocSetRootElement(VDP_CCStructures_xmldoc, ROOT);
        // Save document to file
        VDP_CCStructures_xml.SaveDoc(VDP_CCStructures_xmldoc, Form("./XML_Info_Flute4/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        VDP_CCStructures_xml.FreeDoc(VDP_CCStructures_xmldoc);
      }
    }
  }
}

//==========================================================================================
double VDPoly_Read_Info_from_File(TString VDPoly_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(VDPoly_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << VDPoly_DataFileName << " wit Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;
   char date[20],time[20];
   float Voltage, Current, VDPoly_Temperature, VDPoly_AirTemperature, VDPoly_Hymidity;
   while (in.good()) {

     if(nlines < 16) {
         getline (in,line);
         if(nlines == 13) {
	   cout << "Line Number  = " << nlines << " with content : " <<  line << '\n';
	   string Temperature = line.substr(14, 7);
	   VDPoly_Temp_measurement[id] = std::stof(Temperature);
	   cout << "Temperature = " << VDPoly_Temp_measurement[id] << '\n';
        }
	} else {

          in  >> date >> time >> Current >> Voltage >> VDPoly_Temperature >> VDPoly_AirTemperature >> VDPoly_Hymidity;

	  VDPoly_current_d[N_meas] = Current;
	  VDPoly_voltage_d[N_meas] = Voltage*1E-9;

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
   VDPoly_iv[id] =  new TGraph(N_meas,VDPoly_current_d, VDPoly_voltage_d);
   if((id % 4) == 0) VDPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPoly_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 1) VDPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPoly_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 2) VDPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPoly_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 3) VDPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPoly_r",HM_Id,Structure_Id[file_id].c_str()));
   VDPoly_iv[id]->GetXaxis()->SetTitle("VdP Current [A]");
   VDPoly_iv[id]->GetYaxis()->SetTitle("VdP Voltage [V]");
   VDPoly_iv[id]->SetDrawOption("AP");
   VDPoly_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   VDPoly_iv[id]->SetMarkerColor(colorid);
   //   VDPoly_iv[id]->Draw();


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   VDPoly_cv_low[id] = VDPoly_current_d[3];
   VDPoly_cv_high[id] = VDPoly_current_d[38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", VDPoly_cv_low[id], VDPoly_cv_high[id]);
//   VDPoly_iv[id]->Fit("f_cv","0R+");
   VDPoly_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   VDPoly_A_constant[id] = VDPoly_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   VDPoly_A_slope[id]  = VDPoly_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double Resistance = (TMath::Pi()/TMath::Log(2))*VDPoly_A_slope[id]; // VdP Resistance
   double Resistance_error = (TMath::Pi()/TMath::Log(2))*(VDPoly_iv[id]->GetFunction("f_cv")->GetParError(1));

   cout << "----->>> VDPPoly PQC1 Structure Resistance  = " << Resistance << " +/- " <<  Resistance_error<< " [Ohm/sq]" << endl;

   return Resistance;

}
void VdPoly_PQC_Flute1_Analysis_Final(int nFiles)
{

   gStyle->SetOptStat(0);

   TString VDPoly_DataFileName[600];
   for(int i=0;i<nFiles;i++) {

    VDPoly_DataFileName[4*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPPoly_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDPoly_DataFileName[4*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPPoly_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());

    VDPoly_DataFileName[4*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC1_VDPPoly_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    VDPoly_DataFileName[4*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC1_VDPPoly_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }

   for (int i=0;i<nFiles;i++) {


     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPoly_DataFileName[4*i] << endl;
     VDPoly_s_Rsh_E[i] = VDPoly_Read_Info_from_File(VDPoly_DataFileName[4*i],(4*i),i);
     //VDPoly_s_Rsh_E[i] /=1000.;

     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPoly_DataFileName[4*i+1] << endl;
     VDPoly_s_Rsh_W[i] = VDPoly_Read_Info_from_File(VDPoly_DataFileName[4*i+1],(4*i+1),i);
     //VDPoly_s_Rsh_W[i] /=1000.;

     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPoly_DataFileName[4*i+2] << endl;
     VDPoly_r_Rsh_E[i] = VDPoly_Read_Info_from_File(VDPoly_DataFileName[4*i+2],(4*i+2),i);
     //VDPoly_r_Rsh_E[i] /=1000.;

     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPoly_DataFileName[4*i+3] << endl;
     VDPoly_r_Rsh_W[i] = VDPoly_Read_Info_from_File(VDPoly_DataFileName[4*i+3],(4*i+3),i);
     //VDPoly_r_Rsh_W[i] /=1000.;
   }

   cout <<endl;
   cout << "======================================================================================================================" << endl;
   cout << "---------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -----------------" << endl;
   cout << "==================================================================================================" << endl;

   cout << "      VDPoly   Id                  | VDPoly_s_Rsh [kOhm/sq] | VDPoly_r_Rsh [kOhm/sq]" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) <<  setw(25)  <<  setprecision(3)<< VDPoly_s_Rsh_E[id]*1E+3 <<  setw(25)  << VDPoly_r_Rsh_E[id]   << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) <<  setw(25) <<  setprecision(3) << VDPoly_s_Rsh_W[id]*1E+3 <<  setw(25)  << VDPoly_r_Rsh_W[id]   << endl;
   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;


}

// ------------------------------------------------------------------------
double CBKRPoly_Read_Info_from_File(TString CBKRPoly_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(CBKRPoly_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << CBKRPoly_DataFileName << " wit Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;
   string line = "";
   while (in.good()) {
     getline (in,line);
      if(nlines < 16) {
         if(nlines == 1) {
         cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         CBKRPoly_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         CBKRPoly_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         CBKRPoly_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         CBKRPoly_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         CBKRPoly_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         //std::string data = Conv_Char_To_String(Kind_of_HM_struct_id);
         //boost::to_upper(data);
         // convert string to upper case
         char uppercase[100];
         int i=0;
         for (int i=0; i<=20; i++){
           uppercase[i]=toupper(Kind_of_HM_struct_id[i]);
         }
         CBKRPoly_struct_id.push_back(Conv_Char_To_String(uppercase));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         CBKRPoly_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20];
         sscanf(line.c_str(), "%*s %s", str1);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1);
         cout<<Kind_of_HM_config_id<<endl;
         CBKRPoly_config_id.push_back(Kind_of_HM_config_id);
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
         CBKRPoly_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         CBKRPoly_waiting_time.push_back(Waiting_time_s);
         }
        if(nlines == 12) {
        cout << line << '\n';
        double temp_degC;
        sscanf(line.c_str(), "%*s %lf", &temp_degC);
        cout<<temp_degC<<endl;
        CBKRPoly_temp.push_back(temp_degC);
       }
        if(nlines == 13) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        CBKRPoly_av_temp.push_back(Av_temp_degC);
       }
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double current_namp;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf ", date, time, &current_namp, &voltage, &temp_degC, &air_temp_degC, &rh_prcnt);
      //cout<<date<<","<<time<<","<<voltage<<","<<(float) current_namp<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
   //        in >>  CBKRPoly_voltage >> CBKRPoly_capacitance >> CBKRPoly_conductance;
      //  in >>  CBKRPoly_voltage >> CBKRPoly_capacitance;
      CBKRPoly_voltage_d[id][N_meas] = voltage*1E-9;
      CBKRPoly_current_d[id][N_meas] = current_namp; // this is not current in namp
      CBKRPoly_temp_meas.push_back(std::vector<double>());
      CBKRPoly_temp_meas[id].push_back(temp_degC);
      CBKRPoly_air_temp_meas.push_back(std::vector<double>());
      CBKRPoly_air_temp_meas[id].push_back(air_temp_degC);
      CBKRPoly_rh_prcnt_meas.push_back(std::vector<double>());
      CBKRPoly_rh_prcnt_meas[id].push_back(rh_prcnt);
      CBKRPoly_timestamp_meas.push_back(std::vector<string>());
      CBKRPoly_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      CBKRPoly_voltage[id][N_meas] = float(CBKRPoly_voltage_d[id][N_meas]);
      CBKRPoly_current[id][N_meas] = float(CBKRPoly_current_d[id][N_meas]);

      CBKRPoly_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));


    //        cout << N_meas << ")  "  <<  CBKRPoly_vv[N_meas] << "   " << CBKRPoly_cap[N_meas] << endl;
      cout << N_meas << ")  "  <<  CBKRPoly_voltage[id][N_meas] << "   " << CBKRPoly_current[id][N_meas]  << endl;
      N_meas++;

      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   CBKRPoly_Nmeas.push_back(N_meas);
   TGaxis::SetMaxDigits(3);

//   TCanvas *cc_spline = new TCanvas();
   CBKRPoly_iv[id] =  new TGraph(N_meas,CBKRPoly_current_d[id], CBKRPoly_voltage_d[id]);
   if((id % 4) == 0) CBKRPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC4_CBKRPoly_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 1) CBKRPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC4_CBKRPoly_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 2) CBKRPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC4_CBKRPoly_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 3) CBKRPoly_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC4_CBKRPoly_r",HM_Id,Structure_Id[file_id].c_str()));
   CBKRPoly_iv[id]->GetXaxis()->SetTitle("CBKRPoly Current [A]");
   CBKRPoly_iv[id]->GetYaxis()->SetTitle("CBKRPoly Voltage [V]");
   CBKRPoly_iv[id]->SetDrawOption("AP");
   CBKRPoly_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   CBKRPoly_iv[id]->SetMarkerColor(colorid);
//   CBKRPoly_iv[id]->Draw();


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   CBKRPoly_cv_low[id] = CBKRPoly_current_d[id][3];
   CBKRPoly_cv_high[id] = CBKRPoly_current_d[id][38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", CBKRPoly_cv_low[id], CBKRPoly_cv_high[id]);
   //   CBKRPoly_iv[id]->Fit("f_cv","0R+");
   CBKRPoly_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   CBKRPoly_A_constant[id] = CBKRPoly_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   CBKRPoly_A_slope[id]  = CBKRPoly_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double R_sh = CBKRPoly_A_slope[id]; // CBKRPoly R_sh
   double R_sh_error = (CBKRPoly_iv[id]->GetFunction("f_cv")->GetParError(1));

   cout << "----->>>   R_sh Voltage = " << R_sh << " +/- " <<  R_sh_error<< " [Ohm/sq]" << endl;

   return R_sh;

}
void CBKRPoly_PQC_Flute4_Analysis_Final(int nFiles)
{
   gStyle->SetOptStat(0);

   VdPoly_PQC_Flute1_Analysis_Final(nFiles); // get the R_Geometry

   //TString CBKRPoly_DataFileName[600];
   for(int i=0;i<nFiles;i++) {
    CBKRPoly_DataFileName[4*i+0] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC4_CBKRPoly_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << CBKRPoly_DataFileName[4*i+0] << endl;
    CBKRPoly_DataFileName[4*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC4_CBKRPoly_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << CBKRPoly_DataFileName[4*i+1] << endl;

    CBKRPoly_DataFileName[4*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC4_CBKRPoly_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << CBKRPoly_DataFileName[4*i+2] << endl;
    CBKRPoly_DataFileName[4*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC4_CBKRPoly_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << CBKRPoly_DataFileName[4*i+3] << endl;

   }

   float d = 13; // in um
   float W = 33; // in um

   for (int i=0;i<nFiles;i++) {

     CBKRPoly_s_Rsh_E[i] = CBKRPoly_Read_Info_from_File(CBKRPoly_DataFileName[4*i+0],(4*i+0),i);
     CBKRPoly_s_Rsh_Rgeom_Corrected_E[i] = CBKRPoly_s_Rsh_E[i]  - 4*(VDPoly_s_Rsh_E[i])*d*d/(3*W*W)*(1+d/(2*(W-d)));

     CBKRPoly_s_Rsh_W[i] = CBKRPoly_Read_Info_from_File(CBKRPoly_DataFileName[4*i+1],(4*i+1),i);
     CBKRPoly_s_Rsh_Rgeom_Corrected_W[i] = CBKRPoly_s_Rsh_W[i]  - 4*(VDPoly_s_Rsh_W[i])*d*d/(3*W*W)*(1+d/(2*(W-d)));

     CBKRPoly_r_Rsh_E[i] = CBKRPoly_Read_Info_from_File(CBKRPoly_DataFileName[4*i+2],(4*i+2),i);
     CBKRPoly_r_Rsh_Rgeom_Corrected_E[i] = CBKRPoly_r_Rsh_E[i]  - 4*VDPoly_r_Rsh_E[i]*d*d/(3*W*W)*(1+d/(2*(W-d)));

     CBKRPoly_r_Rsh_W[i] = CBKRPoly_Read_Info_from_File(CBKRPoly_DataFileName[4*i+3],(4*i+3),i);
     CBKRPoly_r_Rsh_Rgeom_Corrected_W[i] = CBKRPoly_r_Rsh_W[i]  - 4*VDPoly_r_Rsh_W[i]*d*d/(3*W*W)*(1+d/(2*(W-d)));

     CBKRPoly_RC_all[4*i+0] =  CBKRPoly_s_Rsh_Rgeom_Corrected_E[i];
     CBKRPoly_RC_all[4*i+1] =  CBKRPoly_s_Rsh_Rgeom_Corrected_W[i];
     CBKRPoly_RC_all[4*i+2] =  CBKRPoly_r_Rsh_Rgeom_Corrected_E[i];
     CBKRPoly_RC_all[4*i+3] =  CBKRPoly_r_Rsh_Rgeom_Corrected_W[i];

     CBKRPoly_R_all[4*i+0] =  CBKRPoly_s_Rsh_E[i];
     CBKRPoly_R_all[4*i+1] =  CBKRPoly_s_Rsh_W[i];
     CBKRPoly_R_all[4*i+2] =  CBKRPoly_r_Rsh_E[i];
     CBKRPoly_R_all[4*i+3] =  CBKRPoly_r_Rsh_W[i];

   }

   TCanvas *CBKRPoly_cc_final[4];
   TMultiGraph *CBKRPoly_mg[4];
   for (int i=0;i<4;i++) {
     CBKRPoly_cc_final[i] = new TCanvas(Form("CBKRPoly_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     CBKRPoly_mg[i] = new TMultiGraph(Form("CBKRPoly_mg_%d",i),"CBKRPoly R_sh");
     if(i == 0) CBKRPoly_mg[i]->SetTitle("CBKRPoly_s_E R_sh ");
     if(i == 1) CBKRPoly_mg[i]->SetTitle("CBKRPoly_s_W R_sh ");
     if(i == 2) CBKRPoly_mg[i]->SetTitle("CBKRPoly_r_E R_sh ");
     if(i == 3) CBKRPoly_mg[i]->SetTitle("CBKRPoly_r_W R_sh ");

    // CBKRPoly_iv[4*i+]->Draw("ALP");
     for(int id=0;id<nFiles;id++) {

        CBKRPoly_mg[i]->Add(CBKRPoly_iv[4*id+i]);

     }
     CBKRPoly_mg[i]->Draw("LPsame");
     CBKRPoly_mg[i]->GetXaxis()->SetTitle("Current [A]");
     CBKRPoly_mg[i]->GetYaxis()->SetTitle("Voltage [V]");
     CBKRPoly_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     CBKRPoly_cc_final[i]->BuildLegend(0.15,0.5,0.5,0.85);


     gPad->Modified();
     gPad->Update();
   }


   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;

   cout << "      CBKRPoly   Id        |  CBKRPoly_s_Rsh [Ohm/sq]   | CBKRPoly_r_Rsh [Ohm/sq]|  CBKRPoly_s_Rsh_Corr [Ohm/sq]   | CBKRPoly_r_Rsh_Corr [Ohm/sq]|" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) << setw(18) <<  setprecision(4) << CBKRPoly_s_Rsh_E[id]  <<  setw(30)
	       << CBKRPoly_r_Rsh_E[id] << setw(30) <<  setprecision(4) << CBKRPoly_s_Rsh_Rgeom_Corrected_E[id]  <<  setw(30)
	       << CBKRPoly_r_Rsh_Rgeom_Corrected_E[id]   << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << setw(18) <<  setprecision(4) << CBKRPoly_s_Rsh_W[id]  <<  setw(30)
	       << CBKRPoly_r_Rsh_W[id]  << setw(30) <<  setprecision(4) << CBKRPoly_s_Rsh_Rgeom_Corrected_W[id]  <<  setw(30)
	       << CBKRPoly_r_Rsh_Rgeom_Corrected_W[id]  << endl;
   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------------------------------------------//
// CBKRPoly xml production
// --------------------------------------------------------------------------------------------------//
void CBKRPoly_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  cout << "CBKRPoly Xml production is on" << endl;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<4; j++){
          //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = CBKRPoly_DataFileName[4*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = CBKRPoly_DataFileName[4*i+j].Last(delimeter2);
         xml_filename = CBKRPoly_DataFileName[4*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[4*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute4/VPX%d/",HM_Id) +xml_filename);

        // First create engine
        TXMLEngine CBKRPoly_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = CBKRPoly_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = CBKRPoly_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = CBKRPoly_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = CBKRPoly_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = CBKRPoly_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = CBKRPoly_xml.NewChild(HEADER, 0, "RUN");
              CBKRPoly_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              CBKRPoly_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              CBKRPoly_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              CBKRPoly_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              CBKRPoly_xml.NewChild(RUN, 0, "INITIATED_BY_USER", CBKRPoly_Operators[4*i+j].c_str());
              CBKRPoly_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", CBKRPoly_Begin_Timestamp[4*i+j].c_str());
              CBKRPoly_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = CBKRPoly_xml.NewChild(ROOT, 0, "DATA_SET");
            CBKRPoly_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            CBKRPoly_xml.NewChild(DATA_SET, 0, "VERSION", "IV_measurement");
            XMLNodePointer_t PART = CBKRPoly_xml.NewChild(DATA_SET, 0, "PART");
              CBKRPoly_xml.NewChild(PART, 0, "NAME_LABEL", CBKRPoly_name_labels[4*i+j].c_str());
              CBKRPoly_xml.NewChild(PART, 0, "KIND_OF_PART", CBKRPoly_kind_of_parts[4*i+j].c_str());
            XMLNodePointer_t DATA = CBKRPoly_xml.NewChild(DATA_SET, 0, "DATA");
              CBKRPoly_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC4");
              CBKRPoly_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", CBKRPoly_struct_id[4*i+j].c_str());
              CBKRPoly_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", CBKRPoly_set_id[4*i+j].c_str());
              CBKRPoly_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", CBKRPoly_config_id[4*i+j].c_str());
              CBKRPoly_xml.NewChild(DATA, 0, "EQUIPMENT", CBKRPoly_equipment[4*i+j].c_str());
              CBKRPoly_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(CBKRPoly_waiting_time[4*i+j], 3).c_str());
              CBKRPoly_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(CBKRPoly_temp[4*i+j], 3).c_str());
              CBKRPoly_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(CBKRPoly_av_temp[4*i+j], 3).c_str());
              CBKRPoly_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", CBKRPoly_struct_id[4*i+j].c_str());
              CBKRPoly_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = CBKRPoly_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = CBKRPoly_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = CBKRPoly_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  CBKRPoly_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  CBKRPoly_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = CBKRPoly_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                CBKRPoly_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                CBKRPoly_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "IV_measurement");
                XMLNodePointer_t PART_CDS1 = CBKRPoly_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  CBKRPoly_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", CBKRPoly_name_labels[4*i+j].c_str());
                  CBKRPoly_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", CBKRPoly_kind_of_parts[4*i+j].c_str());
                  for (int k=0; k<CBKRPoly_Nmeas[4*i+j]-1; k++){
                      XMLNodePointer_t DATA_CDS1 = CBKRPoly_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                      cout<< Conv_float_to_string(CBKRPoly_current[4*i+j][k], 4).c_str() << endl;
                        CBKRPoly_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string_scientific_format(CBKRPoly_voltage[4*i+j][k], 3).c_str());
                        CBKRPoly_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string(CBKRPoly_current[4*i+j][k]*1E+9, 3).c_str());
                        CBKRPoly_xml.NewChild(DATA_CDS1, 0, "TIME", CBKRPoly_timestamp_meas[4*i+j][k].c_str());
                        CBKRPoly_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(CBKRPoly_temp_meas[4*i+j][k], 3).c_str());
                        CBKRPoly_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(CBKRPoly_air_temp_meas[4*i+j][k], 3).c_str());
                        CBKRPoly_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(CBKRPoly_rh_prcnt_meas[4*i+j][k], 3).c_str());
                  }
                XMLNodePointer_t CHILD_DATA_SET2 = CBKRPoly_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = CBKRPoly_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = CBKRPoly_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      CBKRPoly_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      CBKRPoly_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = CBKRPoly_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    CBKRPoly_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    CBKRPoly_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement");
                    XMLNodePointer_t PART_CDS2 = CBKRPoly_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      CBKRPoly_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", CBKRPoly_name_labels[4*i+j].c_str());
                      CBKRPoly_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", CBKRPoly_kind_of_parts[4*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = CBKRPoly_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      CBKRPoly_xml.NewChild(DATA_CDS2, 0, "RC_OHM", Conv_float_to_string_scientific_format(CBKRPoly_RC_all[4*i+j], 5).c_str());
                      CBKRPoly_xml.NewChild(DATA_CDS2, 0, "R_OHM", Conv_float_to_string_scientific_format(CBKRPoly_R_all[4*i+j], 5).c_str());

        XMLDocPointer_t CBKRPoly_xmldoc = CBKRPoly_xml.NewDoc();
        CBKRPoly_xml.DocSetRootElement(CBKRPoly_xmldoc, ROOT);
        // Save document to file
        CBKRPoly_xml.SaveDoc(CBKRPoly_xmldoc, Form("./XML_Info_Flute4/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        CBKRPoly_xml.FreeDoc(CBKRPoly_xmldoc);
      }
    }
  }
}
//==========================================================================================
double VdP_Read_Info_from_File(TString VDPStrip_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(VDPStrip_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << VDPStrip_DataFileName << " wit Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;
   char date[20],time[20];
   float Voltage, Current, VDPStrip_Temperature, VDPStrip_AirTemperature, VDPStrip_Hymidity;
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

          in  >> date >> time >> Current >> Voltage >> VDPStrip_Temperature >> VDPStrip_AirTemperature >> VDPStrip_Hymidity;

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
   double Resistance_error = (TMath::Pi()/TMath::Log(2))*(VDPStrip_iv[id]->GetFunction("f_cv")->GetParError(1));

   cout << "----->>> VDPStrip PQC1 Structure Resistance  = " << Resistance << " +/- " <<  Resistance_error<< " [Ohm/sq]" << endl;

   return Resistance;

}
void VdP_PQC_Flute1_Analysis_Final(int nFiles)
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
     VDPStrip_s_Rsh_E[i] = VdP_Read_Info_from_File(VDPStrip_DataFileName[4*i],(4*i),i);
     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPStrip_DataFileName[4*i+1] << endl;
     VDPStrip_s_Rsh_W[i] = VdP_Read_Info_from_File(VDPStrip_DataFileName[4*i+1],(4*i+1),i);

     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPStrip_DataFileName[4*i+2] << endl;
     VDPStrip_r_Rsh_E[i] = VdP_Read_Info_from_File(VDPStrip_DataFileName[4*i+2],(4*i+2),i);
     cout << "VDPStrip Flute1 Structure: Analyze file : " << VDPStrip_DataFileName[4*i+3] << endl;
     VDPStrip_r_Rsh_W[i] = VdP_Read_Info_from_File(VDPStrip_DataFileName[4*i+3],(4*i+3),i);
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
// ------------------------------------------------------------------------
double CBKRStrip_Read_Info_from_File(TString CBKRStrip_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(CBKRStrip_DataFileName);


   cout << endl << "=============>>>>>> Analyze File : " << CBKRStrip_DataFileName << " wit Id = " << id << endl;

   int nlines = 0;
   int N_meas = 0;
   string line = "";
   while (in.good()) {
     getline (in,line);
      if(nlines < 16) {
         if(nlines == 1) {
         cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         CBKRStrip_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         CBKRStrip_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         CBKRStrip_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         CBKRStrip_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         CBKRStrip_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         //std::string data = Conv_Char_To_String(Kind_of_HM_struct_id);
         //boost::to_upper(data);
         // convert string to upper case
         char uppercase[100];
         int i=0;
         for (int i=0; i<=20; i++){
           uppercase[i]=toupper(Kind_of_HM_struct_id[i]);
         }
         CBKRStrip_struct_id.push_back(Conv_Char_To_String(uppercase));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         CBKRStrip_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20];
         sscanf(line.c_str(), "%*s %s", str1);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1);
         cout<<Kind_of_HM_config_id<<endl;
         CBKRStrip_config_id.push_back(Kind_of_HM_config_id);
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
         CBKRStrip_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         CBKRStrip_waiting_time.push_back(Waiting_time_s);
         }
        if(nlines == 12) {
        cout << line << '\n';
        double temp_degC;
        sscanf(line.c_str(), "%*s %lf", &temp_degC);
        cout<<temp_degC<<endl;
        CBKRStrip_temp.push_back(temp_degC);
       }
        if(nlines == 13) {
        cout << line << '\n';
        double Av_temp_degC;
        sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
        cout<<Av_temp_degC<<endl;
        CBKRStrip_av_temp.push_back(Av_temp_degC);
       }
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double current_namp;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf ", date, time, &current_namp, &voltage, &temp_degC, &air_temp_degC, &rh_prcnt);
      //cout<<date<<","<<time<<","<<voltage<<","<<(float) current_namp<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
   //        in >>  CBKRStrip_voltage >> CBKRStrip_capacitance >> CBKRStrip_conductance;
      //  in >>  CBKRStrip_voltage >> CBKRStrip_capacitance;
      CBKRStrip_voltage_d[id][N_meas] = voltage*1E-9;
      CBKRStrip_current_d[id][N_meas] = current_namp; // this is not current in namp
      CBKRStrip_temp_meas.push_back(std::vector<double>());
      CBKRStrip_temp_meas[id].push_back(temp_degC);
      CBKRStrip_air_temp_meas.push_back(std::vector<double>());
      CBKRStrip_air_temp_meas[id].push_back(air_temp_degC);
      CBKRStrip_rh_prcnt_meas.push_back(std::vector<double>());
      CBKRStrip_rh_prcnt_meas[id].push_back(rh_prcnt);
      CBKRStrip_timestamp_meas.push_back(std::vector<string>());
      CBKRStrip_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      CBKRStrip_voltage[id][N_meas] = float(CBKRStrip_voltage_d[id][N_meas]);
      CBKRStrip_current[id][N_meas] = float(CBKRStrip_current_d[id][N_meas]);

      CBKRStrip_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));


    //        cout << N_meas << ")  "  <<  CBKRStrip_vv[N_meas] << "   " << CBKRStrip_cap[N_meas] << endl;
      cout << N_meas << ")  "  <<  CBKRStrip_voltage[id][N_meas] << "   " << CBKRStrip_current[id][N_meas]  << endl;
      N_meas++;

      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   CBKRStrip_Nmeas.push_back(N_meas);
   TGaxis::SetMaxDigits(3);

//   TCanvas *cc_spline = new TCanvas();
   CBKRStrip_iv[id] =  new TGraph(N_meas,CBKRStrip_current_d[id], CBKRStrip_voltage_d[id]);
   if((id % 4) == 0) CBKRStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC4_CBKRStrip_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 1) CBKRStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC4_CBKRStrip_s",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 2) CBKRStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_E_Left_PQC4_CBKRStrip_r",HM_Id,Structure_Id[file_id].c_str()));
   if((id % 4) == 3) CBKRStrip_iv[id]->SetTitle(Form("VPX%d_0%s_2-S_HM_W_Left_PQC4_CBKRStrip_r",HM_Id,Structure_Id[file_id].c_str()));
   CBKRStrip_iv[id]->GetXaxis()->SetTitle("CBKRStrip Current [A]");
   CBKRStrip_iv[id]->GetYaxis()->SetTitle("CBKRStrip Voltage [V]");
   CBKRStrip_iv[id]->SetDrawOption("AP");
   CBKRStrip_iv[id]->SetMarkerStyle(20+file_id);
   int colorid = 1+file_id;
   if(colorid == 5) colorid = 46;
   CBKRStrip_iv[id]->SetMarkerColor(colorid);
//   CBKRStrip_iv[id]->Draw();


   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   CBKRStrip_cv_low[id] = CBKRStrip_current_d[id][3];
   CBKRStrip_cv_high[id] = CBKRStrip_current_d[id][38];

// ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol1", CBKRStrip_cv_low[id], CBKRStrip_cv_high[id]);
   //   CBKRStrip_iv[id]->Fit("f_cv","0R+");
   CBKRStrip_iv[id]->Fit("f_cv","0R+");


   // Accumulation Region
   // ---------------------
   CBKRStrip_A_constant[id] = CBKRStrip_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   CBKRStrip_A_slope[id]  = CBKRStrip_iv[id]->GetFunction("f_cv")->GetParameter(1); // Accumulation line slope

   double R_sh = CBKRStrip_A_slope[id]; // CBKRStrip R_sh
   double R_sh_error = (CBKRStrip_iv[id]->GetFunction("f_cv")->GetParError(1));

   cout << "----->>>   R_sh Voltage = " << R_sh << " +/- " <<  R_sh_error<< " [Ohm/sq]" << endl;

   return R_sh;

}
void CBKRStrip_PQC_Flute4_Analysis_Final(int nFiles)
{
   gStyle->SetOptStat(0);

   VdP_PQC_Flute1_Analysis_Final(nFiles); // get the R_Geometry

   //TString CBKRStrip_DataFileName[600];
   for(int i=0;i<nFiles;i++) {
    CBKRStrip_DataFileName[4*i+0] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC4_CBKRStrip_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << CBKRStrip_DataFileName[4*i+0] << endl;
    CBKRStrip_DataFileName[4*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC4_CBKRStrip_s.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << CBKRStrip_DataFileName[4*i+1] << endl;

    CBKRStrip_DataFileName[4*i+2] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/VdP/VPX%d_0%s_2-S_HM_E_Left_PQC4_CBKRStrip_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << CBKRStrip_DataFileName[4*i+2] << endl;
    CBKRStrip_DataFileName[4*i+3] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/VdP/VPX%d_0%s_2-S_HM_W_Left_PQC4_CBKRStrip_r.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
    cout << "Analyze file : " << CBKRStrip_DataFileName[4*i+3] << endl;

   }

   float d = 13; // in um
   float W = 33; // in um

   for (int i=0;i<nFiles;i++) {

     CBKRStrip_s_Rsh_E[i] = CBKRStrip_Read_Info_from_File(CBKRStrip_DataFileName[4*i+0],(4*i+0),i);
     CBKRStrip_s_Rsh_Rgeom_Corrected_E[i] = CBKRStrip_s_Rsh_E[i]  - 4*(VDPStrip_s_Rsh_E[i])*d*d/(3*W*W)*(1+d/(2*(W-d)));

     CBKRStrip_s_Rsh_W[i] = CBKRStrip_Read_Info_from_File(CBKRStrip_DataFileName[4*i+1],(4*i+1),i);
     CBKRStrip_s_Rsh_Rgeom_Corrected_W[i] = CBKRStrip_s_Rsh_W[i]  - 4*(VDPStrip_s_Rsh_W[i])*d*d/(3*W*W)*(1+d/(2*(W-d)));

     CBKRStrip_r_Rsh_E[i] = CBKRStrip_Read_Info_from_File(CBKRStrip_DataFileName[4*i+2],(4*i+2),i);
     CBKRStrip_r_Rsh_Rgeom_Corrected_E[i] = CBKRStrip_r_Rsh_E[i]  - 4*VDPStrip_r_Rsh_E[i]*d*d/(3*W*W)*(1+d/(2*(W-d)));

     CBKRStrip_r_Rsh_W[i] = CBKRStrip_Read_Info_from_File(CBKRStrip_DataFileName[4*i+3],(4*i+3),i);
     CBKRStrip_r_Rsh_Rgeom_Corrected_W[i] = CBKRStrip_r_Rsh_W[i]  - 4*VDPStrip_r_Rsh_W[i]*d*d/(3*W*W)*(1+d/(2*(W-d)));

     CBKRStrip_RC_all[4*i+0] =  CBKRStrip_s_Rsh_Rgeom_Corrected_E[i];
     CBKRStrip_RC_all[4*i+1] =  CBKRStrip_s_Rsh_Rgeom_Corrected_W[i];
     CBKRStrip_RC_all[4*i+2] =  CBKRStrip_r_Rsh_Rgeom_Corrected_E[i];
     CBKRStrip_RC_all[4*i+3] =  CBKRStrip_r_Rsh_Rgeom_Corrected_W[i];

     CBKRStrip_R_all[4*i+0] =  CBKRStrip_s_Rsh_E[i];
     CBKRStrip_R_all[4*i+1] =  CBKRStrip_s_Rsh_W[i];
     CBKRStrip_R_all[4*i+2] =  CBKRStrip_r_Rsh_E[i];
     CBKRStrip_R_all[4*i+3] =  CBKRStrip_r_Rsh_W[i];

   }

   TCanvas *CBKRStrip_cc_final[4];
   TMultiGraph *CBKRStrip_mg[4];
   for (int i=0;i<4;i++) {
     CBKRStrip_cc_final[i] = new TCanvas(Form("CBKRStrip_cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     CBKRStrip_mg[i] = new TMultiGraph(Form("CBKRStrip_mg_%d",i),"CBKRStrip R_sh");
     if(i == 0) CBKRStrip_mg[i]->SetTitle("CBKRStrip_s_E R_sh ");
     if(i == 1) CBKRStrip_mg[i]->SetTitle("CBKRStrip_s_W R_sh ");
     if(i == 2) CBKRStrip_mg[i]->SetTitle("CBKRStrip_r_E R_sh ");
     if(i == 3) CBKRStrip_mg[i]->SetTitle("CBKRStrip_r_W R_sh ");

    // CBKRStrip_iv[4*i+]->Draw("ALP");
     for(int id=0;id<nFiles;id++) {

        CBKRStrip_mg[i]->Add(CBKRStrip_iv[4*id+i]);

     }
     CBKRStrip_mg[i]->Draw("LPsame");
     CBKRStrip_mg[i]->GetXaxis()->SetTitle("Current [A]");
     CBKRStrip_mg[i]->GetYaxis()->SetTitle("Voltage [V]");
     CBKRStrip_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     CBKRStrip_cc_final[i]->BuildLegend(0.15,0.5,0.5,0.85);


     gPad->Modified();
     gPad->Update();
   }


   cout <<endl;
   cout << "====================================================================================================================================================================================================" << endl;
   cout << "---------------------------------------------------------------------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  -------------------------------------------------------------------------------" << endl;
   cout << "====================================================================================================================================================================================================" << endl;

   cout << "      CBKRStrip   Id        |  CBKRStrip_s_Rsh [Ohm/sq]   | CBKRStrip_r_Rsh [Ohm/sq]|  CBKRStrip_s_Rsh_Corr [Ohm/sq]   | CBKRStrip_r_Rsh_Corr [Ohm/sq]|" <<endl;
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) << setw(18) <<  setprecision(4) << CBKRStrip_s_Rsh_E[id]  <<  setw(30)
	       << CBKRStrip_r_Rsh_E[id] << setw(30) <<  setprecision(4) << CBKRStrip_s_Rsh_Rgeom_Corrected_E[id]  <<  setw(30)
	       << CBKRStrip_r_Rsh_Rgeom_Corrected_E[id]   << endl;

         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << setw(18) <<  setprecision(4) << CBKRStrip_s_Rsh_W[id]  <<  setw(30)
	       << CBKRStrip_r_Rsh_W[id]  << setw(30) <<  setprecision(4) << CBKRStrip_s_Rsh_Rgeom_Corrected_W[id]  <<  setw(30)
	       << CBKRStrip_r_Rsh_Rgeom_Corrected_W[id]  << endl;
   }
   cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
}
//--------------------------------------------------------------------------------------------------//
// CBKRStrip xml production
// --------------------------------------------------------------------------------------------------//
void CBKRStrip_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  cout << "CBKRStrip Xml production is on" << endl;
  if(xml_on == 1){
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<4; j++){
          //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = CBKRStrip_DataFileName[4*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = CBKRStrip_DataFileName[4*i+j].Last(delimeter2);
         xml_filename = CBKRStrip_DataFileName[4*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[4*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute4/VPX%d/",HM_Id) +xml_filename);

        // First create engine
        TXMLEngine CBKRStrip_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = CBKRStrip_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = CBKRStrip_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = CBKRStrip_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = CBKRStrip_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = CBKRStrip_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = CBKRStrip_xml.NewChild(HEADER, 0, "RUN");
              CBKRStrip_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              CBKRStrip_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              CBKRStrip_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              CBKRStrip_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              CBKRStrip_xml.NewChild(RUN, 0, "INITIATED_BY_USER", CBKRStrip_Operators[4*i+j].c_str());
              CBKRStrip_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", CBKRStrip_Begin_Timestamp[4*i+j].c_str());
              CBKRStrip_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = CBKRStrip_xml.NewChild(ROOT, 0, "DATA_SET");
            CBKRStrip_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            CBKRStrip_xml.NewChild(DATA_SET, 0, "VERSION", "IV_measurement");
            XMLNodePointer_t PART = CBKRStrip_xml.NewChild(DATA_SET, 0, "PART");
              CBKRStrip_xml.NewChild(PART, 0, "NAME_LABEL", CBKRStrip_name_labels[4*i+j].c_str());
              CBKRStrip_xml.NewChild(PART, 0, "KIND_OF_PART", CBKRStrip_kind_of_parts[4*i+j].c_str());
            XMLNodePointer_t DATA = CBKRStrip_xml.NewChild(DATA_SET, 0, "DATA");
              CBKRStrip_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC4");
              CBKRStrip_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", CBKRStrip_struct_id[4*i+j].c_str());
              CBKRStrip_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", CBKRStrip_set_id[4*i+j].c_str());
              CBKRStrip_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", CBKRStrip_config_id[4*i+j].c_str());
              CBKRStrip_xml.NewChild(DATA, 0, "EQUIPMENT", CBKRStrip_equipment[4*i+j].c_str());
              CBKRStrip_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(CBKRStrip_waiting_time[4*i+j], 3).c_str());
              CBKRStrip_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(CBKRStrip_temp[4*i+j], 3).c_str());
              CBKRStrip_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(CBKRStrip_av_temp[4*i+j], 3).c_str());
              CBKRStrip_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", CBKRStrip_struct_id[4*i+j].c_str());
              CBKRStrip_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = CBKRStrip_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = CBKRStrip_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = CBKRStrip_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  CBKRStrip_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  CBKRStrip_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = CBKRStrip_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                CBKRStrip_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                CBKRStrip_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "IV_measurement");
                XMLNodePointer_t PART_CDS1 = CBKRStrip_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  CBKRStrip_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", CBKRStrip_name_labels[4*i+j].c_str());
                  CBKRStrip_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", CBKRStrip_kind_of_parts[4*i+j].c_str());
                  for (int k=0; k<CBKRStrip_Nmeas[4*i+j]-1; k++){
                      XMLNodePointer_t DATA_CDS1 = CBKRStrip_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                      cout<< Conv_float_to_string(CBKRStrip_current[4*i+j][k], 4).c_str() << endl;
                        CBKRStrip_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string_scientific_format(CBKRStrip_voltage[4*i+j][k], 3).c_str());
                        CBKRStrip_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(CBKRStrip_current[4*i+j][k]*1E+9, 3).c_str());
                        CBKRStrip_xml.NewChild(DATA_CDS1, 0, "TIME", CBKRStrip_timestamp_meas[4*i+j][k].c_str());
                        CBKRStrip_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(CBKRStrip_temp_meas[4*i+j][k], 3).c_str());
                        CBKRStrip_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(CBKRStrip_air_temp_meas[4*i+j][k], 3).c_str());
                        CBKRStrip_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(CBKRStrip_rh_prcnt_meas[4*i+j][k], 3).c_str());
                  }
                XMLNodePointer_t CHILD_DATA_SET2 = CBKRStrip_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = CBKRStrip_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = CBKRStrip_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      CBKRStrip_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      CBKRStrip_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = CBKRStrip_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    CBKRStrip_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    CBKRStrip_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement");
                    XMLNodePointer_t PART_CDS2 = CBKRStrip_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      CBKRStrip_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", CBKRStrip_name_labels[4*i+j].c_str());
                      CBKRStrip_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", CBKRStrip_kind_of_parts[4*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = CBKRStrip_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      CBKRStrip_xml.NewChild(DATA_CDS2, 0, "RC_OHM", Conv_float_to_string(CBKRStrip_RC_all[4*i+j], 3).c_str());
                      CBKRStrip_xml.NewChild(DATA_CDS2, 0, "R_OHM", Conv_float_to_string(CBKRStrip_R_all[4*i+j], 3).c_str());

        XMLDocPointer_t CBKRStrip_xmldoc = CBKRStrip_xml.NewDoc();
        CBKRStrip_xml.DocSetRootElement(CBKRStrip_xmldoc, ROOT);
        // Save document to file
        CBKRStrip_xml.SaveDoc(CBKRStrip_xmldoc, Form("./XML_Info_Flute4/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        CBKRStrip_xml.FreeDoc(CBKRStrip_xmldoc);
      }
    }
  }
}
//==========================================================================================

void GCD05_Read_IV_from_File(TString GCD05_DataFileName, int id, int file_id)
{
   ifstream in;
   in.open(GCD05_DataFileName);

   epsilon_0 = 0.0885418782; // [pF/cm]
   epsilon_Si = 11.68*epsilon_0; // [pF/cm]
   epsilon_SiO2 = 3.9*epsilon_0; // [pF/cm]
   GCD05_area = 0.723e-2; // [cm^2]
   //GCD05_area = 0.723e-2; // [cm^2]
   q = 1.60219E-7; // [pF.V]
   kB = 8.6173324E-5; // [eV/K]
   //   TKelvin = 273 + Temp_measurement[id]; // [K]
   //   Nintrisic = 1.45e10; // [cm^-3]
   cout << endl << "=============>>>>>> Analyze File : " << GCD05_DataFileName << endl;

   int nlines = 0;
   int N_meas = 0;

   string line;
   while (in.good()) {
     getline (in,line);
      if(nlines < 18) {
         if(nlines == 1) {
         cout << line << '\n';
         char First_Name[20], Last_Name[20];
         sscanf(line.c_str(), "%*s %s %s", First_Name, Last_Name);
         cout<<First_Name<<Last_Name<<endl;
         GCD05_Operators.push_back(Conv_Char_To_String(First_Name)+" "+Conv_Char_To_String(Last_Name));
         }
         if(nlines == 2) {
         cout << line << '\n';
         char Date[20], Time[20];
         sscanf(line.c_str(), "%*s %s %s", Date, Time);
         cout<<Date<<Time<<endl;
         GCD05_Begin_Timestamp.push_back(Conv_Char_To_String(Date)+" "+Conv_Char_To_String(Time));
         }
         if(nlines == 3) {
         cout << line << '\n';
         char Name_Label[25];
         sscanf(line.c_str(), "%*s %s", Name_Label);
         cout<<Name_Label<<endl;
         GCD05_name_labels.push_back(Conv_Char_To_String(Name_Label));
         }
         if(nlines == 4) {
         cout << line << '\n';
         char batch_type[20], loc[20];
         sscanf(line.c_str(), "%*s %s %*s %s",batch_type, loc);
         string Kind_of_part = std::string()+Conv_Char_To_String(batch_type)+" Halfmoon "+loc;
         cout<<Kind_of_part<<endl;
         GCD05_kind_of_parts.push_back(Kind_of_part);
         }
         if(nlines == 5) {
         cout << line << '\n';
         char Kind_of_HM_flute_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_flute_id);
         cout<<Kind_of_HM_flute_id<<endl;
         GCD05_kind_of_HM_flute_id.push_back(Kind_of_HM_flute_id);
         }
         if(nlines == 6) {
         cout << line << '\n';
         char Kind_of_HM_struct_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_struct_id);
         cout<<Kind_of_HM_struct_id<<endl;
         //std::string data = Conv_Char_To_String(Kind_of_HM_struct_id);
         //boost::to_upper(data);
         // convert string to upper case
         char uppercase[100];
         int i=0;
         for (int i=0; i<=20; i++){
           uppercase[i]=toupper(Kind_of_HM_struct_id[i]);
         }
         GCD05_struct_id.push_back(Conv_Char_To_String(uppercase));
         }
         if(nlines == 7) {
         cout << line << '\n';
         char Kind_of_HM_set_id[20];
         sscanf(line.c_str(), "%*s %s", Kind_of_HM_set_id);
         cout<<Kind_of_HM_set_id<<endl;
         GCD05_set_id.push_back(Conv_Char_To_String(Kind_of_HM_set_id));
         }
         if(nlines == 8) {
         cout << line << '\n';
         char str1[20], str2[20];
         sscanf(line.c_str(), "%*s %s %s", str1, str2);
         string Kind_of_HM_config_id = Conv_Char_To_String(str1)+" "+Conv_Char_To_String(str2);
         cout<<Kind_of_HM_config_id<<endl;
         GCD05_config_id.push_back(Kind_of_HM_config_id);
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
         GCD05_equipment.push_back(Equipment);
         }
         if(nlines == 11) {
         cout << line << '\n';
         double Waiting_time_s;
         sscanf(line.c_str(), "%*s %lf", &Waiting_time_s);
         cout<<Waiting_time_s<<endl;
         GCD05_waiting_time.push_back(Waiting_time_s);
         }
        if(nlines == 12) {
        cout << line << '\n';
        double bias;
        sscanf(line.c_str(), "%*s %lf", &bias);
        cout<<bias<<endl;
        GCD05_Bias.push_back(bias);
       }
         if(nlines == 13) {
         cout << line << '\n';
         double temp_degC;
         sscanf(line.c_str(), "%*s %lf", &temp_degC);
         cout<<temp_degC<<endl;
         GCD05_temp.push_back(temp_degC);
         TKelvin = 273 + temp_degC;
       }
         if(nlines == 14) {
         cout << line << '\n';
         double Av_temp_degC;
         sscanf(line.c_str(), "%*s %lf", &Av_temp_degC);
         cout<<Av_temp_degC<<endl;
         GCD05_av_temp.push_back(Av_temp_degC);
        }
        if(nlines == 15) {
        cout << line << '\n';
        double compliance;
        sscanf(line.c_str(), "%*s %lf", &compliance);
        cout<<compliance<<endl;
        GCD05_compliance.push_back(compliance);
       }
      }
      else {
      //  cout<< line.c_str() << endl;
      char date[20], time[20];
      double voltage, temp_degC, air_temp_degC, rh_prcnt;
      double current_namp;
      sscanf(line.c_str(), "%s %s %lf %lf %lf %lf %lf ", date, time, &voltage, &current_namp, &temp_degC, &air_temp_degC, &rh_prcnt);
      //cout<<date<<","<<time<<","<<voltage<<","<<(float) current_namp<<","<<temp_degC<<","<<air_temp_degC<<","<<rh_prcnt<< endl;
   //        in >>  GCD05_vv >> GCD05_capacitance >> GCD05_conductance;
      //  in >>  GCD05_vv >> GCD05_capacitance;
      GCD05_vv_d[id][N_meas] = voltage;
      GCD05_current_d[id][N_meas] = current_namp*1E+3; // this is not current in namp
      GCD05_temp_meas.push_back(std::vector<double>());
      GCD05_temp_meas[id].push_back(temp_degC);
      GCD05_air_temp_meas.push_back(std::vector<double>());
      GCD05_air_temp_meas[id].push_back(air_temp_degC);
      GCD05_rh_prcnt_meas.push_back(std::vector<double>());
      GCD05_rh_prcnt_meas[id].push_back(rh_prcnt);
      GCD05_timestamp_meas.push_back(std::vector<string>());
      GCD05_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));
      GCD05_vv[id][N_meas] = float(GCD05_vv_d[id][N_meas]);
      GCD05_current[id][N_meas] = float(GCD05_current_d[id][N_meas]);

      GCD05_timestamp_meas[id].push_back(Conv_Char_To_String(date)+" "+Conv_Char_To_String(time));


    //        cout << N_meas << ")  "  <<  GCD05_vv[N_meas] << "   " << GCD05_cap[N_meas] << endl;
      cout << N_meas << ")  "  <<  GCD05_vv[id][N_meas] << "   " << GCD05_current[id][N_meas]  << endl;
      N_meas++;

      }
      if (!in.good()) break;
      nlines++;
   }
   printf(" Found %d lines\n",nlines);
   in.close();
   GCD05_Nmeas.push_back(N_meas);
   cout << endl << "=====================>>>>>>>>>>>>>>>>>>>>>>>   Number of Measurements = " << N_meas << endl << endl;

   //TKelvin = 273 + GCD05_Temp_measurement[id]; // [K]
   Nintrisic = 5.29E+19*pow((TKelvin/300),2.54)*exp(-6726/TKelvin);

   cout << endl << ">>>>>>>>>>> Nintrisic  = " << Nintrisic << " cm^{-3} and area = " << GCD05_area <<  " cm^{-2} " << endl;
   int  N_meas_reduced = N_meas-1;

   if(N_meas_reduced > 300 ) {
      cout << "Serious Problem : Measurements > 300 that TSpline3 can handle...  " << endl;
      cout << " Solution: Reduce number of measurements for TSpline3 rootine.... " << endl;
      cout << " Now Abort..." << endl;
      abort();
   }
   double GCD05_vv_d_reduced[300];
   double GCD05_current_d_reduced[300];
   for (Int_t i=0; i<N_meas_reduced; i++) {
      GCD05_vv_d_reduced[i] = GCD05_vv_d[id][i];
      GCD05_current_d_reduced[i] = GCD05_current_d[id][i];
      //      cout << i << ")  "  <<  GCD05_vv_d[id][i] << "   " << GCD05_current_d[id][i] << endl;
   }
    //   CSpline_Interpolator[id] = new ROOT::Math::Interpolator();
    //   CSpline_Interpolator[id]->SetData(N_meas_reduced,vv_d_reduced,current_d_inversed);
    GCD05_CuSpl_iv[id] = new TSpline3("Cubic Spline", GCD05_vv_d_reduced, GCD05_current_d_reduced, N_meas_reduced, "cb2e2", 0, 0);  // Seems that works well for less than 300 points
    //GCD05_cc_spline[id] = new TCanvas(Form("GCD05_cc_spline_%d",id),"", 700, 500, 700,700);
    GCD05_iv[id] =  new TGraph(N_meas-2,GCD05_vv_d_reduced,GCD05_current_d_reduced);
    if((id % 2) == 0) GCD05_iv[id]->SetTitle(Form("GCD05%d_0%s E IV",HM_Id,Structure_Id[file_id].c_str()));
    if((id % 2) == 1) GCD05_iv[id]->SetTitle(Form("GCD05%d_0%s W IV",HM_Id,Structure_Id[file_id].c_str()));
    GCD05_iv[id]->GetXaxis()->SetTitle("Gate Voltage [V]");
    GCD05_iv[id]->GetYaxis()->SetTitle("GCD05 Current [pA]");
    GCD05_iv[id]->GetYaxis()->SetRangeUser(-1,15);
    GCD05_iv[id]->SetDrawOption("AP");
    GCD05_iv[id]->SetMarkerStyle(20+file_id);
    int colorid = 1+file_id;
    if(colorid == 5) colorid = 46;
    GCD05_iv[id]->SetMarkerColor(colorid);
    //GCD05_iv[id]->Draw();
    GCD05_CuSpl_iv[id]->SetLineColor(kBlue);
    //GCD05_CuSpl_iv[id]->Draw("same");
    // Get V_Fb grosso - modo
    // --------------------------
    for(int i = 0;i<N_meas_reduced; i++) {
       GCD05_dev_spline_simple[i] = GCD05_CuSpl_iv[id]->Derivative(GCD05_vv_d_reduced[i]);
    }

    //GCD05_cc_spline_der[id] = new TCanvas(Form("GCD05_cc_spline_der_%d",id),"", 50, 200, 500,500);
    TGraph *GCD05_iv_deriv =  new TGraph(N_meas_reduced,GCD05_vv_d_reduced,GCD05_dev_spline_simple);
    if((id % 2) == 0) GCD05_iv_deriv->SetTitle(Form("GCD05%d_0%s E IV derivarive",HM_Id,Structure_Id[file_id].c_str()));
    if((id % 2) == 1) GCD05_iv_deriv->SetTitle(Form("GCD05%d_0%s W IV derivarive",HM_Id,Structure_Id[file_id].c_str()));
    //GCD05_iv_deriv->Draw();

    /*
        for (Int_t i=0; i<N_meas_reduced; i++) {

          cout << "Data points : " << i << ")  "  <<  GCD05_vv_d_reduced[i] << "   " << "  "
               << GCD05_current_d_reduced[i] <<  "  " << GCD05_CuSpl_iv[id]->Eval(GCD05_vv_d_reduced[i]) << "   " << GCD05_dev_spline_simple[i] << endl;

        }
    */
    // get the minimum of the derivative
    int maxdev_index_simple = 0;
    float max_GCD05_dev_spline_simple = 0;
    for (int i = 1;i<N_meas_reduced; i++) {
       if(GCD05_dev_spline_simple[i] >  max_GCD05_dev_spline_simple && GCD05_vv_d_reduced[i] > -9 &&  GCD05_vv_d_reduced[i] < -7) {
         max_GCD05_dev_spline_simple = GCD05_dev_spline_simple[i];
	 maxdev_index_simple = i;
       }
    }
    cout  << "---------------- >>>>> Max derivative Simple at i = " << maxdev_index_simple
          << " with value = " << GCD05_dev_spline_simple[maxdev_index_simple] << endl;

    // Get V_Fb fine tuned
    // -------------------------
    int Nub_steps = 1000;

    float vv_initial = GCD05_vv_d_reduced[0];
    float vv_final = GCD05_vv_d_reduced[N_meas-1];

    Double_t dev_spline[2000];
    Float_t dev_spline_f[2000];
    Float_t voltage_extented[2000];
    for(int i = 0;i<Nub_steps; i++) {
      voltage_extented[i] =  vv_initial + i*(vv_final - vv_initial)/(Nub_steps-1);
      dev_spline[i] = GCD05_CuSpl_iv[id]->Derivative(voltage_extented[i]);
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
    double V_FB_derivative = double(voltage_extented[max_index]) - GCD05_Bias[id];
    double I_FB_derivative =GCD05_CuSpl_iv[id]->Eval(V_FB_derivative);
    cout  << "---------------- >>>>> Max derivative accurate at i = " << max_index
          << "  with value = " << dev_spline_f[max_index]
	        << "  V_flatband = " <<  V_FB_derivative
	        << " and I_flatband = " << I_FB_derivative << endl;

    double V_FB_plus_V_Bias_derivative = double(voltage_extented[min_index]);
    double I_FB_plus_V_Bias_derivative =GCD05_CuSpl_iv[id]->Eval(V_FB_plus_V_Bias_derivative);
    cout  << "---------------- >>>>> Min derivative accurate at i = " << min_index
          << "  with value = " << dev_spline_f[min_index]
	        << " V_flatband = " <<  V_FB_plus_V_Bias_derivative
	        << " and I_flatband_plus_B_Bias = " << I_FB_plus_V_Bias_derivative << endl << endl;

   GCD05_VFB_Derivative[id] = V_FB_derivative; // in V
   cout << endl<<  "======>>> V_flatband = " <<  V_FB_derivative << endl;

   // find the maximum of the IV curve in range maxdev_index_simple + 16
   int IV_max_index = 0;
   double IV_Max = 0.0;
   for (int i = maxdev_index_simple;i<maxdev_index_simple+16; i++) {
       if(GCD05_current_d[id][i] > IV_Max ) {
        IV_Max  = GCD05_current_d[id][i];
	      IV_max_index  = i;
       }
    }
    cout  << "---------------- >>>>> Max IV at i = " <<  IV_max_index
          << " with value = " << IV_Max << endl;
   // Set linear fit ranges
   // -------------------------

   // Accumulation Region Linear Fit with slop set to 0
   // ---------------------------------------------------

   GCD05_iv_low[id] = GCD05_vv_d_reduced[IV_max_index-20];
   GCD05_iv_high[id] = GCD05_vv_d_reduced[IV_max_index-30];

   // Depletion Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   int Dep_Central_Point = IV_max_index;

   GCD05_iv1_low[id] = GCD05_vv_d_reduced[Dep_Central_Point];
   GCD05_iv1_high[id] = GCD05_vv_d_reduced[Dep_Central_Point+5];


   // Inversion Region Linear Fit with slop set to 0
   // ---------------------------------------------------
   GCD05_iv2_low[id] = GCD05_vv_d_reduced[N_meas-60];;
   GCD05_iv2_high[id] = GCD05_vv_d_reduced[N_meas-45];

    // Get tree branches to find the fit points
    // -------------------------------------------

    //   Int_t nentries = (Int_t)T->GetEntries();
    //   cout << "Number of Measurements = " << nentries << endl;

    // ---------------------------------------------------------------------------------------------------------------

   TF1 *f_cv = new TF1("f_cv","pol0", GCD05_iv_low[id], GCD05_iv_high[id]);
//   f_cv->SetParameter(1,0.0);

   GCD05_iv[id]->Fit("f_cv","R+");

   TF1 *f_cv1 = new TF1("f_cv1","pol0",GCD05_iv1_low[id], GCD05_iv1_high[id]);

   GCD05_iv[id]->Fit("f_cv1","R+");

   TF1 *f_cv2 = new TF1("f_cv2","pol0", GCD05_iv2_low[id], GCD05_iv2_high[id]);

   GCD05_iv[id]->Fit("f_cv2","R+");

   // Accumulation Region
   // ---------------------
   GCD05_A_constant[id] = GCD05_iv[id]->GetFunction("f_cv")->GetParameter(0); // Accumulation constant term
   GCD05_A_slope[id]  = 0.0; // Accumulation line slope

   GCD05_I_accumulation[id] = GCD05_A_constant[id]; // GCD05 Accumulation Current
   float I_accumulation_Error = GCD05_iv[id]->GetFunction("f_cv")->GetParError(0);

   cout << endl << "======>>>   Accumulation Current = " <<  GCD05_I_accumulation[id] << " +/- " << I_accumulation_Error << " [pA]" << endl;

   // Depletion Region
   // ---------------------
   GCD05_D_constant[id] = GCD05_iv[id]->GetFunction("f_cv1")->GetParameter(0); // Depletion Region constant term
   GCD05_D_slope[id] = 0.0; // Depletion Region slope

   GCD05_I_depletion[id] = GCD05_D_constant[id]; // GCD05 Depletion Current
   float GCD05_I_depletion_Error = GCD05_iv[id]->GetFunction("f_cv1")->GetParError(0);

   cout << "======>>>   Depletion Current = " << GCD05_I_depletion[id] << " +/- " << GCD05_I_depletion_Error << " [pA]" << endl;

   // Inversion Region
   // ---------------------
   GCD05_I_constant[id] = GCD05_iv[id]->GetFunction("f_cv2")->GetParameter(0); // Inversion constant term
   GCD05_I_slope[id]  = 0.0; // Inversion line slope

   GCD05_I_inversion[id] = GCD05_I_constant[id]; // GCD05 High Frequency Current - Inversion Current
   float GCD05_I_inversion_Error = GCD05_iv[id]->GetFunction("f_cv2")->GetParError(0);

   cout << "======>>>   Inversion Current = " << GCD05_I_inversion[id] << " +/- " << GCD05_I_inversion_Error << " [pA]" << endl;

   GCD05_I_surface[id] = GCD05_I_depletion[id] - GCD05_I_inversion[id];

   cout << "======>>>   I_surface = " << GCD05_I_surface[id] << " +/- "
        << sqrt(GCD05_I_depletion_Error*GCD05_I_depletion_Error + GCD05_I_inversion_Error*GCD05_I_inversion_Error) << " [pA]" << endl;

   GCD05_S_interface_recombination_Velocity[id] = GCD05_I_surface[id]/(q*Nintrisic*GCD05_area);

   cout << "======>>>   Interface Recombination Velocity (S_{0}) = " << GCD05_S_interface_recombination_Velocity[id] << " [cm/sec]" <<endl;
   cout << "======>>>   V_flatband = " <<  V_FB_derivative << endl;

   float sigma = 2.5e-16; // [cm^2]
   float u_thermal = 2.0e7; // [cm/sec]
   GCD05_D_it[id] = GCD05_S_interface_recombination_Velocity[id]/(3.14159 * sigma * u_thermal *  TKelvin * kB);
   GCD05_N_it[id] = 1.12*GCD05_D_it[id]/2;;


}
void GCD05_PQC_Analysis_Final(int nFiles)
{
   gStyle->SetOptStat(0);

   //TString GCD05_DataFileName[100];
   for(int i=0;i<nFiles;i++) {
      GCD05_DataFileName[2*i] = Data_Dir + Form("VPX%d_2-S_HM_E/txt_files/IV/VPX%d_0%s_2-S_HM_E_Left_PQC4_GCD05_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
      GCD05_DataFileName[2*i+1] = Data_Dir + Form("VPX%d_2-S_HM_W/txt_files/IV/VPX%d_0%s_2-S_HM_W_Left_PQC4_GCD05_IV.txt",HM_Id,HM_Id,Structure_Id[i].c_str());
   }


   for (int i=0;i<nFiles;i++) {

     GCD05_Read_IV_from_File(GCD05_DataFileName[2*i],2*i,i);
     GCD05_Read_IV_from_File(GCD05_DataFileName[2*i+1],2*i+1,i);
   }

   TCanvas *GCD05_cc_final[2];
   TMultiGraph *GCD05_mg[2];
   for (int i=0;i<2;i++) {
     GCD05_cc_final[i] = new TCanvas(Form("cc_final_%d",i),"",50*i,10,1000,1000);
     gPad->SetGrid();

     GCD05_mg[i] = new TMultiGraph(Form("GCD05_mg_%d",i),"GCD05 IV");
     if(i == 0) GCD05_mg[i]->SetTitle("GCD05_E IV ");
     if(i == 1) GCD05_mg[i]->SetTitle("GCD05_W IV ");

     for(int id=0;id<nFiles;id++) {

        GCD05_mg[i]->Add(GCD05_iv[2*id+i]);

     }
     GCD05_mg[i]->Draw("LPsame");
     GCD05_mg[i]->GetXaxis()->SetTitle("Gate GCD05_vv [V]");
     GCD05_mg[i]->GetYaxis()->SetTitle("GCD05_current [pA]");
     GCD05_mg[i]->GetYaxis()->SetRangeUser(-1,15);
     GCD05_mg[i]->GetYaxis()->SetTitleOffset(1.35);
     gPad->Modified();
     gPad->Update();

     GCD05_cc_final[i]->BuildLegend(0.6,0.6,0.85,0.85);
     gPad->Modified();
     gPad->Update();
   }

   cout <<endl;
   cout << "=====================================================================================" << endl;
   cout << "------------ VPX28442: HMSet_Left GCD05 Rectangular Irradiated  -------------------------" << endl;
   cout << "=====================================================================================" << endl;

   cout << "         GCD05 Id                 | I_acc [pA]   | I_dep [pA]   |  I_inv [pA]  | I_Sur [pA]   | S_0 [cm/sec]    | V_FB [V]  |     N_it " << endl;
   cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

   for(int id=0;id<nFiles;id++) {

         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       <<  setw(15) <<  setprecision(4) << GCD05_I_depletion[2*id]
	       <<  setw(15) <<  setprecision(4) << GCD05_I_depletion[2*id]
	       <<  setw(15) <<  setprecision(2) << GCD05_I_inversion[2*id]
	       <<  setw(15) <<  setprecision(3)<< GCD05_I_surface[2*id]
	       <<  setw(15) << GCD05_S_interface_recombination_Velocity[2*id]
	       <<  setw(15) << GCD05_VFB_Derivative[2*id]
	       <<  setw(15) << GCD05_N_it[2*id]    << endl;


         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       <<  setw(15) <<  setprecision(4) << GCD05_I_depletion[2*id+1]
	       <<  setw(15) <<  setprecision(4) << GCD05_I_depletion[2*id+1]
	       <<  setw(15) <<  setprecision(2) << GCD05_I_inversion[2*id+1]
	       <<  setw(15) <<  setprecision(3)<< GCD05_I_surface[2*id+1]
	       <<  setw(15) << GCD05_S_interface_recombination_Velocity[2*id+1]
	       <<  setw(15) << GCD05_VFB_Derivative[2*id+1]
	       <<  setw(15) << GCD05_N_it[2*id+1]    << endl;


   }
   cout << "-------------------------------------------------------------------------------------" << endl;

}
void GCD05_PQC_xml_production(int xml_on, int nFiles){
  TString xml_filename;
  if(xml_on == 1){
    cout << "GCD05 Xml production is on" << endl;
    for(int i=0;i<nFiles;i++) {
      for(int j=0; j<=1; j++){
        //       TString xml_filename[nFiles];
         char delimeter1 = '/';
         int pos1 = GCD05_DataFileName[2*i+j].Last(delimeter1);
         char delimeter2 = '.';
         int pos2 = GCD05_DataFileName[2*i+j].Last(delimeter2);
         xml_filename = GCD05_DataFileName[2*i+j](pos1+1, ((pos2-pos1)-1)) + ".xml";
         //        cout << MOS_DataFileName[2*i+j](0,pos1) << endl;

         cout << xml_filename << endl;
         ofstream xml_file;

         xml_file.open(Form("./XML_Info_Flute4/VPX%d/",HM_Id) +xml_filename);

         // First create engine
        TXMLEngine GCD05_xml;

        // Create main node of document tree
        XMLNodePointer_t ROOT = GCD05_xml.NewChild(0, 0, "ROOT");
          XMLNodePointer_t HEADER = GCD05_xml.NewChild(ROOT, 0, "HEADER");
            XMLNodePointer_t TYPE = GCD05_xml.NewChild(HEADER, 0, "TYPE");
              XMLNodePointer_t EXTENSION_TABLE_NAME = GCD05_xml.NewChild(TYPE, 0, "EXTENSION_TABLE_NAME", "HALFMOON_METADATA");
              XMLNodePointer_t NAME = GCD05_xml.NewChild(TYPE, 0, "NAME", "Tracker Halfmoon Metadata");
            XMLNodePointer_t RUN = GCD05_xml.NewChild(HEADER, 0, "RUN");
              GCD05_xml.NewAttr(RUN, 0, "mode", "SEQUENCE_NUMBER");
              GCD05_xml.NewAttr(RUN, 0, "sequence", "TRK_OT_RUN_SEQ");
              GCD05_xml.NewChild(RUN, 0, "RUN_TYPE", "PQC");
              GCD05_xml.NewChild(RUN, 0, "LOCATION", "Demokritos");
              GCD05_xml.NewChild(RUN, 0, "INITIATED_BY_USER", GCD05_Operators[2*i+j].c_str());
              GCD05_xml.NewChild(RUN, 0, "RUN_BEGIN_TIMESTAMP", GCD05_Begin_Timestamp[2*i+j].c_str());
              GCD05_xml.NewChild(RUN, 0, "COMMENT_DESCRIPTION", "Test");
          XMLNodePointer_t DATA_SET = GCD05_xml.NewChild(ROOT, 0, "DATA_SET");
            GCD05_xml.NewChild(DATA_SET, 0, "COMMENT_DESCRIPTION", "Metadata with flute and structure");
            GCD05_xml.NewChild(DATA_SET, 0, "VERSION", "v2");
            XMLNodePointer_t PART = GCD05_xml.NewChild(DATA_SET, 0, "PART");
              GCD05_xml.NewChild(PART, 0, "NAME_LABEL", GCD05_name_labels[2*i+j].c_str());
              GCD05_xml.NewChild(PART, 0, "KIND_OF_PART", GCD05_kind_of_parts[2*i+j].c_str());
            XMLNodePointer_t DATA = GCD05_xml.NewChild(DATA_SET, 0, "DATA");
              GCD05_xml.NewChild(DATA, 0, "KIND_OF_HM_FLUTE_ID", "PQC4");
              GCD05_xml.NewChild(DATA, 0, "KIND_OF_HM_STRUCT_ID", "GCD05");
              GCD05_xml.NewChild(DATA, 0, "KIND_OF_HM_SET_ID", GCD05_set_id[2*i+j]);
              GCD05_xml.NewChild(DATA, 0, "KIND_OF_HM_CONFIG_ID", GCD05_config_id[2*i+j].c_str());
              GCD05_xml.NewChild(DATA, 0, "EQUIPMENT", GCD05_equipment[2*i+j].c_str());
              GCD05_xml.NewChild(DATA, 0, "WAITING_TIME_S", Conv_float_to_string(GCD05_waiting_time[2*i+j],3).c_str());
              GCD05_xml.NewChild(DATA, 0, "TEMP_SET_DEGC", Conv_float_to_string(GCD05_temp[2*i+j],3).c_str());
                cout << "Xml production is on" << endl;
              GCD05_xml.NewChild(DATA, 0, "AV_TEMP_DEGC", Conv_float_to_string(GCD05_av_temp[2*i+j],3).c_str());
              GCD05_xml.NewChild(DATA, 0, "PROCEDURE_TYPE", "GCD051");
              GCD05_xml.NewChild(DATA, 0, "FILE_NAME", "Demokritos_"+xml_filename);
            XMLNodePointer_t CHILD_DATA_SET1 = GCD05_xml.NewChild(DATA_SET, 0, "CHILD_DATA_SET");
              XMLNodePointer_t HEADER_CDS1 = GCD05_xml.NewChild(CHILD_DATA_SET1, 0, "HEADER");
                XMLNodePointer_t TYPE_CDS1 = GCD05_xml.NewChild(HEADER_CDS1, 0, "TYPE");
                  GCD05_xml.NewChild(TYPE_CDS1, 0, "EXTENSION_TABLE_NAME", "TEST_SENSOR_IV");
                  GCD05_xml.NewChild(TYPE_CDS1, 0, "NAME", "Tracker Halfmoon IV Test");
              XMLNodePointer_t DATA_SET_CDS1 = GCD05_xml.NewChild(CHILD_DATA_SET1, 0, "DATA_SET");
                GCD05_xml.NewChild(DATA_SET_CDS1, 0, "COMMENT_DESCRIPTION", "Test");
                GCD05_xml.NewChild(DATA_SET_CDS1, 0, "VERSION", "V2");
                XMLNodePointer_t PART_CDS1 = GCD05_xml.NewChild(DATA_SET_CDS1, 0, "PART");
                  GCD05_xml.NewChild(PART_CDS1, 0, "NAME_LABEL", GCD05_name_labels[2*i+j].c_str());
                  GCD05_xml.NewChild(PART_CDS1, 0, "KIND_OF_PART", GCD05_kind_of_parts[2*i+j].c_str());
                for (int k=0; k<GCD05_Nmeas[2*i+j]-1; k++){
                XMLNodePointer_t DATA_CDS1 = GCD05_xml.NewChild(DATA_SET_CDS1, 0, "DATA");
                  GCD05_xml.NewChild(DATA_CDS1, 0, "VOLTS", Conv_float_to_string(GCD05_vv_d[2*i+j][k],3).c_str());
                  //cout << Conv_float_to_string(GCD05_vv_d[id]k],4).c_str() << endl;
                  //GCD05_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", to_string(GCD05_current[k]).c_str());
                  GCD05_xml.NewChild(DATA_CDS1, 0, "CURRNT_NAMP", Conv_float_to_string_scientific_format(GCD05_current[2*i+j][k]*1E-3, 3).c_str()); // The current for gthe analysis is in pA
                  GCD05_xml.NewChild(DATA_CDS1, 0, "TIME", GCD05_timestamp_meas[2*i+j][k].c_str());
                  GCD05_xml.NewChild(DATA_CDS1, 0, "TEMP_DEGC", Conv_float_to_string(GCD05_temp_meas[2*i+j][k], 3).c_str());
                  GCD05_xml.NewChild(DATA_CDS1, 0, "AIR_TEMP_DEGC", Conv_float_to_string(GCD05_air_temp_meas[2*i+j][k], 3).c_str());
                  GCD05_xml.NewChild(DATA_CDS1, 0, "RH_PRCNT", Conv_float_to_string(GCD05_rh_prcnt_meas[2*i+j][k], 3).c_str());
                }
                XMLNodePointer_t CHILD_DATA_SET2 = GCD05_xml.NewChild(DATA_SET_CDS1, 0, "CHILD_DATA_SET");
                  XMLNodePointer_t HEADER_CDS2 = GCD05_xml.NewChild(CHILD_DATA_SET2, 0, "HEADER");
                    XMLNodePointer_t TYPE_CDS2 = GCD05_xml.NewChild(HEADER_CDS2, 0, "TYPE");
                      GCD05_xml.NewChild(TYPE_CDS2, 0, "EXTENSION_TABLE_NAME", "HALFMOON_IV_PAR");
                      GCD05_xml.NewChild(TYPE_CDS2, 0, "NAME", "Tracker Halfmoon IV Parameters");
                  XMLNodePointer_t DATA_SET_CDS2 = GCD05_xml.NewChild(CHILD_DATA_SET2, 0, "DATA_SET");
                    GCD05_xml.NewChild(DATA_SET_CDS2, 0, "COMMENT_DESCRIPTION", "Test");
                    GCD05_xml.NewChild(DATA_SET_CDS2, 0, "VERSION", "IV_measurement-004");
                    XMLNodePointer_t PART_CDS2 = GCD05_xml.NewChild(DATA_SET_CDS2, 0, "PART");
                      GCD05_xml.NewChild(PART_CDS2, 0, "NAME_LABEL", GCD05_name_labels[2*i+j].c_str());
                      GCD05_xml.NewChild(PART_CDS2, 0, "KIND_OF_PART", GCD05_kind_of_parts[2*i+j].c_str());
                    XMLNodePointer_t DATA_CDS2 = GCD05_xml.NewChild(DATA_SET_CDS2, 0, "DATA");
                      GCD05_xml.NewChild(DATA_CDS2, 0, "ISURF_PAMPR", Conv_float_to_string(GCD05_I_surface[2*i+j], 3).c_str());
                      GCD05_xml.NewChild(DATA_CDS2, 0, "S0_CMSEC", Conv_float_to_string(GCD05_S_interface_recombination_Velocity[2*i+j], 3).c_str());
                      GCD05_xml.NewChild(DATA_CDS2, 0, "VFB_ACC_V", Conv_float_to_string(GCD05_VFB_Derivative[2*i+j], 3).c_str());
        XMLDocPointer_t GCD05_xmldoc = GCD05_xml.NewDoc();
        GCD05_xml.DocSetRootElement(GCD05_xmldoc, ROOT);
        // Save document to file
        GCD05_xml.SaveDoc(GCD05_xmldoc, Form("./XML_Info_Flute4/VPX%d/",HM_Id) +xml_filename);
        // Release memory before exit
        GCD05_xml.FreeDoc(GCD05_xmldoc);
      }
    }
  }
}


void Flute4_2_S_Extended_Characterization_with_xml(int nFiles)
{
   VDP_CCStructures_PQC_Flute4_Analysis_Final(nFiles);
   VDP_CCStructures_PQC_xml_production(xml_on, nFiles);

   CBKRPoly_PQC_Flute4_Analysis_Final(nFiles);
   CBKRPoly_PQC_xml_production(xml_on, nFiles);

   CBKRStrip_PQC_Flute4_Analysis_Final(nFiles);
   CBKRStrip_PQC_xml_production(xml_on, nFiles);

   GCD05_PQC_Analysis_Final(nFiles);
   GCD05_PQC_xml_production(xml_on, nFiles);

   CVS_Output_File.open(Form("CSV_Info/VPX%d/Flute4_VPX%d_0xx_2-S_HM_E_W_Left_Extended_Info.csv",HM_Id,HM_Id));

   cout << "==========================================================================================================================================" << endl;
   cout << "--------------------------  " <<  Form("VPX%d_0xx_2-S_HM_E/W_Left",HM_Id) << "  ----------------------------------------------------------" << endl;
   cout << "==========================================================================================================================================" << endl;
   cout << "     Flute4  Id            | CCEdge Rc | CCPoly Rc | CCStrip Rc |   CBKRPoly   |  CBKRStrip   |  GCD05 I_surf | GCD05 S_0  | GCD05 Vfb_Acc_Dep |" <<endl;
   cout << "                           |  [kOhm]   |  [MOhm]   |   [kOhm]   |   [kOhm/sq]  |   [Ohm/sq]  |     [pA]    | [cm/sec] |       [V]       |" <<endl;
   cout << "                           |           |           |            |  stad | rot  |  stad | rot  |             |          |                 |" <<endl ;
   cout << "     Spec Limit            |           |           |            |              |              |             |   0.60   |       <5        |" <<endl ;
   cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
   cout << showpoint;
   for(int id=0;id<nFiles;id++) {
         cout  <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(9) <<  setprecision(3) << VDP_CCEdge_Rc_E[id]*1E-3
	       << setw(12) << VDP_CCPoly_Rc_E[id]*1E-6
	       << setw(13) << VDP_CCStrip_Rc_E[id]*1E-3
	       << setw(11) << setprecision(4) << CBKRPoly_s_Rsh_Rgeom_Corrected_E[id]*1E-3
	       << setw(7) <<  CBKRPoly_r_Rsh_Rgeom_Corrected_E[id]*1E-3
	       << setw(8) <<  setprecision(4) << CBKRStrip_s_Rsh_Rgeom_Corrected_E[id]
	       << setw(7) <<  CBKRStrip_r_Rsh_Rgeom_Corrected_E[id]
	       << setw(11) <<  setprecision(3)<< GCD05_I_surface[2*id]
	       << setw(12) << GCD05_S_interface_recombination_Velocity[2*id]
	       << setw(14) << GCD05_VFB_Derivative[2*id]
	       << endl;


         cout  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str())
	       << setw(9) <<  setprecision(3) << VDP_CCEdge_Rc_W[id]*1E-3
	       << setw(12) << VDP_CCPoly_Rc_W[id]*1E-6
	       << setw(13) << VDP_CCStrip_Rc_W[id]*1E-3
	       << setw(11) << setprecision(4) << CBKRPoly_s_Rsh_Rgeom_Corrected_W[id]*1E-3
	       << setw(7) <<  CBKRPoly_r_Rsh_Rgeom_Corrected_W[id]*1E-3
	       << setw(8) <<  setprecision(4) << CBKRStrip_s_Rsh_Rgeom_Corrected_W[id]
	       << setw(7) <<  CBKRStrip_r_Rsh_Rgeom_Corrected_W[id]
	       << setw(11) <<  setprecision(3)<< GCD05_I_surface[2*id+1]
	       << setw(12) << GCD05_S_interface_recombination_Velocity[2*id+1]
	       << setw(14) << GCD05_VFB_Derivative[2*id+1]
	       << endl;
   }

   cout << "----------------------------------------------------------------------------------------------------------------------------------------------" << endl;

   CVS_Output_File << "     Flute4  Id            ;Mode ; CCEdge Rc ; CCPoly Rc ; CCStrip Rc  ; CBKRPoly;     ; CBKRStrip ;  ;  GCD05 I_surf ;  GCD05 S_0   ; GCD05 Vfb_Acc_Dep" << endl;
   CVS_Output_File << "                           ;     ;  [kOhm]   ; [MOhm]    ; [kOhm]      ; [kOhm/sq]  ;  ; [Ohm/sq] ;  ;   [pA]      ;  [cm/sec]  ;    [V]      "<< endl;
   CVS_Output_File << "                           ;     ;           ;           ;             ; stad    ; rot ; stad  ;  rot ;             ;            ;              " << endl;
   CVS_Output_File << "     Spec Limit            ;     ;           ;           ;             ;         ;     ;       ;      ;             ;            ;  " << endl;

   for(int id=0;id<nFiles;id++) {
           CVS_Output_File <<  Form("VPX%d_0%s_2-S_HM_E_Left",HM_Id,Structure_Id[id].c_str()) << ";" << "Extended" << ";"
	                   << VDP_CCEdge_Rc_E[id]*1E-3 << ";" << VDP_CCPoly_Rc_E[id]*1E-6  << ";" << VDP_CCStrip_Rc_E[id]*1E-3 << ";"
	                   << CBKRPoly_s_Rsh_Rgeom_Corrected_E[id]*1E-3  << ";" <<  CBKRPoly_r_Rsh_Rgeom_Corrected_E[id]*1E-3  << ";"
			               << CBKRStrip_s_Rsh_Rgeom_Corrected_E[id] << ";" <<  CBKRStrip_r_Rsh_Rgeom_Corrected_E[id]  << ";"
	                   << GCD05_I_surface[2*id] << ";" << GCD05_S_interface_recombination_Velocity[2*id] << ";"  << GCD05_VFB_Derivative[2*id]
	                   << endl;


          CVS_Output_File  <<  Form("VPX%d_0%s_2-S_HM_W_Left",HM_Id,Structure_Id[id].c_str()) << ";" << "Extended" << ";"
	                   << VDP_CCEdge_Rc_W[id]*1E-3 << ";" << VDP_CCPoly_Rc_W[id]*1E-6 << ";" << VDP_CCStrip_Rc_W[id]*1E-3 << ";"
	                   << CBKRPoly_s_Rsh_Rgeom_Corrected_W[id]*1E-3  << ";" <<  CBKRPoly_r_Rsh_Rgeom_Corrected_W[id]*1E-3  << ";"
			               << CBKRStrip_s_Rsh_Rgeom_Corrected_W[id] << ";" <<  CBKRStrip_r_Rsh_Rgeom_Corrected_W[id]  << ";"
	                   << GCD05_I_surface[2*id+1] << ";" << GCD05_S_interface_recombination_Velocity[2*id+1] << ";"  << GCD05_VFB_Derivative[2*id+1]
	                   << endl;
   }

}
