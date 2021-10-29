#!/bin/bash
if [ $# -lt 2 ]
then
  echo 'The script needs two arguments:'
  echo 'The First argument is the Flute Id folder eg: XML_Info_Flute1'
  echo 'The Second argument is the HM Batch ID folder eg: VPX37400'
  echo 'Thus the correct syndax is : > source Upload_files_to_DB.sh  XML_Info_Flute1 VPX37400'
  echo 'Thus try again -- Exit'  
else 

  export current_dir=`pwd`

  echo "I am in folder : ${current_dir}"

  echo "I will upload all files of folder : $1/S2"

  cd $1/$2/

  for file in *.xml
  do
    echo "I will upload in DB file :  ${file}"
    #python3 ${current_dir}/cmsdbldr_client.py --login --url=https://cmsdca.cern.ch/trk_loader/trker/cmsr  ${file}
  done 
  cd ${current_dir}
  ll

fi
