#!/bin/bash
#AssemblyForEveryone.sh developed by Yiwen Zhou (Univeristy of Adelaide)
#+with help from Dr Jimmy Breen, Dr Stephen Pederson, Dr Rick Tearle
#the purpose of this script is for my own project.
#some bash bullitins may have compatability problem between Mac and Linux! be cautcious

usage () {           #####needs to add information
echo "1. this scripts only works for Paired-end reads(FR). ./AssemblyForEveryone.sh [forward_reads.fastq.gz] [reverse_reads.fastq.gz]
2. when type more than one Args, plz put white space in between"
}

echo "#command for the first step              "1. Fastqc check"
fastqc -t 2 -q \$Read1 \$Read2
#command for the second step             "2. Preproccess raw data"
java -jar ~/Downloads/trimmomatic-0.36/trimmomatic-0.36.jar PE \$Read1 \$Read2 \${Read1/%fastq.gz/Fw_Paired.fastq.gz} \${Read1/%fastq.gz/Fw_UnPaired.fastq.gz} \${Read2/%fastq.gz/Rv_Paired.fastq.gz} \${Read2/%fastq.gz/Rv_UnPaired.fastq.gz} MINLEN:35
#command for the third step              "3. Assembly"
spades.py -o \$location2 --sc --pe1-1 \${Preproccess1-\$global_var1} --pe1-2 \${Preproccess2-\$global_var2} --cov-cutoff auto
#command for the forth step              "4. Extend contigs"

#command for the fifth step              "5. Re-assembly with raw reads."

#command for the sixth step              "6. Evaluation"
reapr smaltmap -k 13 -s 2 -n 4 -m 14 \$facheck \$Read_fq1 \$Read_fq2 \$location/\$Dir6/\${ANSWER1}_ReapR.bam
reapr pipeline \$Final_contigs \$mapped_file \$(pwd)/ReapR_pipeline_results  " > $(pwd)/AFEconfig.txt
#exit 0

x1=$1
x2=$2
func () {
  global_var1=$x1
  global_var2=$x2
}
func

#Read1=$global_var1
#Read2=$global_var2
Read1=$(pwd)/$global_var1
Read2=$(pwd)/$global_var2

if [ "$1" == --help ]
   then
      usage
      exit
fi

echo "      ############################################################
      ###In order to achieve disirable results, please finish#####
      ###saving config file in the same folder as scripts file####
      ############################################################"

      #echo "Checking Files $L1 $L2"
      #if [ ! -f "$1" ]
      #    then
      #      echo "plz check the file $1"
      #fi

      #if [ ! -f "$2" ];
      #    then
      #      echo "plz check the file $2"
      #      exit
      #fi
##############################################
#####   starting with quality check of   #####
#####   the raw data                     #####
##############################################
echo
# ${var/%position/pattern}
echo "Give your directory a name OR Enter <00> to choose provious created directories"
read ANSWER1
Date=$(date '+%y-%m-%d')

if [ "$ANSWER1" == 00 ]

then
        find ~/ -type d -name "10-*-*" | xargs -n 1 echo "" > $(pwd)/list.txt
        nl -b a $(pwd)/list.txt # show the lane number
        echo "Choose the lane number for the wanted Directory to start "
        read ANSWERR
        Choose=($(sed -n "${ANSWERR}p" "$(pwd)/list.txt"))
        echo "=====> Moving to Directory <$Choose>..."
        cd $Choose

else
       Dir="10-$ANSWER1-$Date"
       mkdir -p $(pwd)/$Dir
       echo "=====> Creating Directory <$Dir>..."
       mkdir $Dir
       cd $(pwd)/$Dir
echo
fi

#if [ ! $? -eq 0 ]
#          then
#              echo "=====> Creating Directory <$Dir> failed"
#              exit 1
#          else
#              echo "=====> Creating Directory <$Dir> succeeded"
#fi

#echo "=====> making <$ANSWER1-$Date> completed"

#this step is for jumping each section of scripts
echo "===========Specify the number to start with==========="
echo "0. this allow you run the whole pipeline"
echo "1. Fastqc check (QC) -- Fastqc"
echo "2. Preproccess raw data (PR) -- Trimmomatic"
echo "3. Assembly (AS) -- Spades"
echo "4. Extend contigs (EC) -- (Not available)"
echo "5. Re-assembly with raw reads (RA)-- (Not available)."
echo "6. Evaluation (EA) -- ReapR"
echo "7. Report (Not available)"
echo
read ANSWER ANSWER2 ANSWER3 ANSWER4 ANSWER5 ANSWER6 #$ANSWER is for navigation

function QC () {
echo "#=================================   1 start  =================================#"
echo
   echo "we are now at" ; pwd
      location=$(pwd)
        echo
        #  cd $location/$Dir
            mkdir Fastqc_Check
              Dir1="Fastqc_Check"
                echo ">>>>===Start changing Directory to <$location/$Dir1>"
                  cd $location/$Dir1
                    echo
              echo ">>>>===Changing Directory to <$location/$Dir1> succeeded!!!"
echo
echo
echo "##############################################
Starting with quality check for the raw data
=== Paired-end Forward $global_var1
=== Paired-end Reverse $global_var2
##############################################"
echo
                    echo "searching for file AFEconfig.txt"
                       find "$location" -name "AFEconfig.txt"
                       sample_info=($(find "$location/.." -name "AFEconfig.txt"))
                          if [ "$?" -eq 0 ] ##double quates are really important
                            then
                              echo "AFEconfig.txt found"
                            else
                              echo "cannot file AFEconfig.txt"
                              exit
                          fi

                    echo "====>Reading AFEconfig.txt file ..."
                 sample_info1=($(cut -f 1 "$sample_info" | sed -n '2p'))
                 # echo $sample_info1 (for debugging step)
                    which ${sample_info1}
                        if [ "$?" -eq 0 ]
                          then
                            echo "====>Calling ${sample_info1} succeeded"
                          else
                            echo "####Unable to excute, cannot find ${sample_info1} Plz Check wether ${sample_info1}
                            is in your variable path otherwise state the path to the tool in the config file####"
                            exit
                        fi
                #Tool_Matrix=($(cut -f 1 "$sample_info"))
                       echo "$location"
                       eval $(sed -n '2p' "$sample_info")
                                         #      2> QCerror.txt
                #        if [ -s QCerror.txt ] ##negate true if the file size greater than 0
                #            then
                #              echo "######${sample_info1} works failed, plz check the QCerror.txt file######"
                #              exit
                #            else
                #              ls -lh
                find "$location/.." -maxdepth 1 -name "*fastqc*" | xargs -n 1 -I {} mv -f {} $location/$Dir1

                              echo "#=================================   1 end    =================================#"
                    #    fi
}

#                                                                              #
#                                                                              #
#echo "Checking Files $L1 $L2"                                                 #
#if [ ! -f "$1" ]                                                              #
#    then                                                                      #
#      echo "plz check the file $1"                                            #
#fi                                                                            #
                                                                              #
#if [ ! -f "$2" ];                                                             #
#    then                                                                      #
#      echo "plz check the file $2"                                            #
#      exit                                                                    #
#fi                                                                            #
#                                                                              #
function PR () {
echo "#=================================   2 start   ================================#"
echo
   #location=$(pwd)
   #cd $location/$Dir
   echo "we are now at" ; pwd
  # cd $(pwd)/..
   #echo $location
   #exit 0
   if [ ! -d Fastqc_Check ]
   then
     location=$(pwd)
  #   cd $location/$Dir
         #  exit 0
        echo
          mkdir -p $location/Preprocess_Raw_Data
            Dir2="Preprocess_Raw_Data"
        echo ">>>>===Start changing Directory to <$location/$Dir2>"
                echo
                cd $location/$Dir2
                location1=$(pwd)
              #  echo $location1
              #  exit 0
              echo ">>>>===Changing Directory to <$location1> succeeded 111 !!!"
      #     exit 0
    else
       location=$(pwd)
       echo
       mkdir $location/Preprocess_Raw_Data
         Dir2="Preprocess_Raw_Data"
            cd $location/$Dir2
               location1=$(pwd)
              echo ">>>>===Start changing Directory to <$location1>"
                cd $location1
                   echo
              echo ">>>>===Changing Directory to <$location1> succeeded!!!"
      #exit 0
    fi
echo
echo
echo "##############################################
Starting with Preprocess_Raw_Data for the raw data
=== Paired-end Forward $global_var1
=== Paired-end Reverse $global_var2
##############################################"
echo
    echo "searching for file AFEconfig.txt"
      find "$location/.." -name "AFEconfig.txt"
        sample_info=($(find "$location/.." -name "AFEconfig.txt"))
          echo "Reading AFEconfig.txt file ..."
            sample_info2=($(cut -f 1 "$sample_info" | sed -n '4p'))
            # echo $sample_info1 (for debugging step)
            which ${sample_info2}
                   if [ "$?" -eq 0 ]  ##double quates are really important
                    then
                     echo
                    else
                     echo "####Unable to excute, cannot find ${sample_info2} Plz make sure that ${sample_info2}
                     is in your variable path otherwise state the path to the tool in the config file####"
                     exit
                   fi
           #Tool_Matrix=($(cut -f 1 "$sample_info"))
           eval "$(sed -n '4p' "$sample_info")" > PRerror.txt

           if [ -s PRerror.txt ] ##negate true if the file size greater than 0
               then
                 echo "######${sample_info2} works failed, plz check the PRerror.txt file######"
                 exit 0
               else
                # contents=($(find "$location" -name "*_Paired*"))
                    #     until [ -w "$contents"]
                    #          do
			      echo $location
                              find "$location/.." -maxdepth 1 -name "*_Paired*" | xargs -n 1 -I {} mv {} $location/$Dir2

                              find "$location/$Dir2" -name "*_Paired*" | xargs -n 1 -I {} echo {} > $location/$Dir2/reads_contents.txt

                              find "$location/.." -name "*_UnPaired.fastq.gz" | xargs -n 1 -I {} rm -f {}

                              ls -lh

                              echo "#=================================   2 end   ================================#"
                    #          done
           fi
}


function AS () {
echo "#=================================   3 start   ================================#"
echo
echo "we are now at" ; pwd
#exit 0
if [ ! -d Preprocess_Raw_Data ]
then
  location=$(pwd)
#   cd $location/$Dir
      #  exit 0
     echo
       mkdir -p $location/Assembly
         Dir3="Assembly"
         echo ">>>>===Start changing Directory to <$location/$Dir3>"
          echo
             cd $location/$Dir3
             location2=$(pwd)
           #  echo $location1
           #  exit 0
         echo ">>>>===Changing Directory to <$location2> succeeded 111 !!!"
   #     exit 0
 else
    location=$(pwd) # $(pwd)/..
    echo
    mkdir -p $location/Assembly
      Dir3="Preprocess_Raw_Data"
      Dir3="Assembly"
         cd $location/$Dir3
            location2=$(pwd)
            echo ">>>>===Start changing Directory to <$location2>"
             cd $location2
                echo
                     echo ">>>>===Changing Directory to <$location2> succeeded!!!"
           read1=($(sed -n '1p' $location/$Dir2/reads_contents.txt))
           read2=($(sed -n '2p' $location/$Dir2/reads_contents.txt))
 #  exit 0
 fi
#if [ -s "${ANSWER1}_Paired_reads_contents.txt"]
#then
  echo ${read1-$global_var1}
  echo ${read2-$global_var2}
  #exit 0
  Preproccess1=${read1-$global_var1}
  Preproccess2=${read2-$global_var2}
  echo "##############################################
  Starting with Assemlby for the raw data
  === Paired-end Forward $Preproccess1
  === Paired-end Reverse $Preproccess2
  ##############################################"
#else
#  echo "Using supplied input files $global_var1 $global_var2"
#fi
echo
echo
echo
    echo "searching for file AFEconfig.txt"
      find "$location/.." -name "AFEconfig.txt"
        sample_info=($(find "$location/.." -name "AFEconfig.txt"))
          echo "Reading AFEconfig.txt file ..."
            sample_info3=($(cut -f 1 "$sample_info" | sed -n '6p'))
            # echo $sample_info1 (for debugging step)
            which ${sample_info3}
                   if [ "$?" -eq 0 ]  ##double quates are really important
                    then
                     echo
                    else
                     echo "####Unable to excute, cannot find ${sample_info3} Plz make sure that ${sample_info3}
                     is in your variable path otherwise state the path to the tool in the config file####"
                     exit 0
                   fi
            eval "$(sed -n '6p' "$sample_info")" 2> ASerror.txt
            if [ -s ASerror.txt ] ##negate true if the file size greater than 0
                then
                  echo "######${sample_info3} works failed, plz check the ASerror.txt file######"
                  exit
                else
                  ls -lh
                  echo "#=================================   3 end    =================================#"
            fi
          }

function EC () {
echo "#=================================   4 start   ================================#"
echo
echo "we are now at" ; pwd
#exit 0
SSPACE=($(find "~/" -iname "SSPACE_Standard_v*.pl"))

 if [ ! -d Assembly ]
 then
contigs=$global_var1
     location=$(pwd)
       echo
         mkdir Extend_contigs
           Dir4="Extend_contigs"
             echo ">>>>===Start changing Directory to <$location/$Dir4>"
               cd $location/$Dir4
                 echo
             echo ">>>>===Changing Directory to <$location/$Dir4> succeeded!!!"
echo
#exit 0
else
  Dir3="Assembly"
  contigs=($(find "$location/$Dir3" -maxdepth 1 -name "*scaffold*fa"))
        mkdir Extend_contigs
            Dir4="Extend_contigs"
                echo ">>>>===Start changing Directory to <$location/$Dir4>"
                      cd $location/$Dir4
                        echo
                            echo ">>>>===Changing Directory to <$location/$Dir4> succeeded!!!"
                              contigs=($(find "$location/$Dir3" -maxdepth 1 -name "*scaffold*fa"))
                                  echo $contigs

fi

echo
echo "##############################################
Starting with Extend contigs for the data
=== contigs file       $contigs
##############################################"
echo
    echo "searching for file AFEconfig.txt"
      find "$location" -name "AFEconfig.txt"
        sample_info=($(find "$location/.." -name "AFEconfig.txt"))
          echo "Reading AFEconfig.txt file ..."
            sample_info4=($(cut -f 1 "$sample_info" | sed -n '8p'))
        which ${sample_info4}
            if [ "$?" -eq 0 ]  ##double quates are really important
              then
                echo
              else
                echo "####Unable to excute, cannot find ${sample_info4} Plz make sure that ${sample_info4}
                is in your variable path otherwise state the path to the tool in the config file####"
                exit 0
            fi

            eval "$(sed -n '8p' "$sample_info")" 2> ECerror.txt
              if [ -s ECerror.txt ] ##negate true if the file size greater than 0
                then
                  echo "######${sample_info4} works failed, plz check the ECerror.txt file######"
                  exit
                else
                  ls -lh
                  echo "#=================================   4 end (skip 5)   =================================#"
                fi
}

#function RA () {
#echo "#=================================   5 start   ================================#"
#if [ ! -d Extend_contigs ]
#then
#echo "we are now at" ; pwd
#contigs=$global_var1
#    location=$(pwd)
#      echo
#else
#  contigs=($(find "$location/$Dir3" -maxdepth 1 -name "*scaffold*fa"))
#  echo $contigs
#fi
#      mkdir Re_Assembly_withRaw
#        Dir5="Re_Assembly_withRaw"
#            echo ">>>>===Start changing Directory to <$location/$Dir5>"
#              cd $location/$Dir5
#                echo
#            echo ">>>>===Changing Directory to <$location/$Dir5> succeeded!!!"
#echo
#exit 0
#echo
#echo
#echo "##############################################
#Starting with RE-Assembly back to contigs for the raw data
#=== Paired-end Forward ${Preproccess1-$global_var1}
#=== Paired-end Reverse ${Preproccess2-$global_var2}
#=== Congtigs file      $contigs
##############################################"
#echo
#    echo "searching for file AFEconfig.txt"
#      find "$location/.." -name "AFEconfig.txt"
#        sample_info=($(find "$location/.." -name "AFEconfig.txt"))
#          echo "Reading AFEconfig.txt file ..."
#            sample_info5=($(cut -f 1 "$sample_info" | sed -n '10p'))
#            which ${sample_info5}
#              if [ "$?" -eq 0 ]  ##double quotes are really important
#                then
#                  echo
##                  echo "####Unable to excute, cannot find ${sample_info5} Plz make sure that ${sample_info5}
  #                is in your variable path otherwise state the path to the tool in the config file####"
  #              fi
#
#                eval "$(sed -n '10p' "$sample_info")" 2> RAerror.txt
#                  if [ -s RAerror.txt ] ##negate true if the file size greater than 0
#                    then
#                      echo "######${sample_info5} works failed, plz check the RAerror.txt file######"
#                      exit
#                    else
#                      echo "#=================================   5 end    =================================#"
#                      ls -lh
#                  fi
#}


function EA () {
echo "#=================================   6 start   ================================#"
echo
#fastq1=($(find $location -samefile "$Rea
#fastq2=($(find $location -samefile "$Read2" | xargs -I {} gunzip {}))
#echo $fastq1
#echo $fastq2
#exit 0
echo "we are now at" ; pwd
   location=$(pwd)
     echo
       mkdir Evaluation
	 Dir3="Assembly"
	 Dir4="Extend_contigs"
         Dir6="Evaluation"
       zcat $(find ~/ -name $global_var1) > $(pwd)/../${global_var1/%.gz}
       zcat $(find ~/ -name $global_var2) > $(pwd)/../${global_var2/%.gz}

       Read_fq1=($(find "$(pwd)/.." -name "${global_var1/%.gz}"))
       Read_fq2=($(find "$(pwd)/.." -name "${global_var2/%.gz}"))
#echo ${Read_fq1}
#echo ${Read_fq2}


if [ ! -d Extend_contigs ]

then
          contigs=($(find "$location/$Dir3" -maxdepth 1 -iname "*scaffold*fasta"))
       echo $contigs
#	exit 0
 #       sed -n '1~4s/^@/>/p;2~4p' $contigs > $location/$Dir6/Spades_scaffold.fasta

else
          contigs1=($(find "$location/$Dir4" -maxdepth 1 -iname "*SSPACE*scaffold*.fasta"))
#          sed -n '1~4s/^@/>/p;2~4p' $contigs1 > $location/$Dir6/SSPACE_extended_scaffold.fasta

fi
#exit 0

#R1_fq=($(gunzip ${Preproccess1-$global_var1}))
#R2_fq=($(gunzip ${Preproccess2-$global_var2}))
             echo ">>>>===Start changing Directory to <$location/$Dir6>"
               cd $location/$Dir6
                 echo
             echo ">>>>===Changing Directory to <$location/$Dir6> succeeded!!!"
echo
Final_contigs=${contigs-$contigs1}
echo "Final contigs $Final_contigs"
echo "##############################################
Starting with Evaluation Preprocess
=== Contigs files      $Final_contigs
=== Forward Raw Reads  $Read1
=== Reversed Raw Reads $Read2
##############################################"
echo
echo "searching for file AFEconfig.txt"
  find "$location" -name "AFEconfig.txt"
    sample_info=($(find "$location/.." -name "AFEconfig.txt"))
      echo "Reading AFEconfig.txt file ..."
        sample_info6=($(cut -f 1 "$sample_info" | sed -n '12p'))
        which ${sample_info6}
     #   exit 0
              if [ "$?" -eq 0 ]  ##double quates are really important
                 then
                 echo
                 else
                 echo "####Unable to excute, cannot find ${sample_info6} Plz make sure that ${sample_info6}
                 is in your variable path otherwise state the path to the tool in the config file####"
              fi

           # eval "$(sed -n '12p' "$sample_info")" 2> EAerror.txt
           # waiit
           reapr facheck $Final_contigs scaffold.facheck
	   wait
	   echo "facheck done"
	#	mv "$location/../"SPades.facheck.fasta $location/$Dir6/
            facheck=$(find $location/$Dir6/ -name "*facheck.fa")
#exit 0
	    eval "$(sed -n '12p' "$sample_info")"
            mapped_file=($(find $location/$Dir6/ -name "*.bam"))

      #      eval "$(sed -n '11p' "$sample_info")" 2>> EAerror.txt
        #    reapr smaltmap -k 13 -s 2 -n 4 -m 14 $facheck $Read_fq1 $Read_fq2 ${Final_contigs/%fasta/bam}
            wait
           echo "ReapR pipeline start"
         reapr pipeline $Final_contigs $mapped_file $(pwd) 2>> EAerror.txt
	   echo "ReapR pipeline finished"
              if [ -s EAerror.txt ] ##negate true if the file size greater than 0
                 then
                    echo "######${sample_info6} works failed, plz check the EAerror.txt file######"
                      exit
                         else
                  ls -lh
              fi

              if [ "$?" -gt 0 ]
                 then
                    echo "######${sample_info6} works failed, plz check the EAerror.txt file######"
                     exit 0
 #                      echo "#=================================   6 end    =================================#"
	        else
		ls -lh
echo "#=================================   6 end    =================================#"
 fi
          }

function_arrary=(exit QC PR AS EC RA EA)
#                 0   1  2  3  4  5  6
G=0
###this is for selecting one of the step
#for i in $function_arrary[@]
#do

if [ $ANSWER -eq 0 ]
          then
              echo "=====> starting pipeline>>"
              eval ${function_arrary[1]} ; wait
              cd $(pwd)/..
              eval ${function_arrary[2]} ; wait
              cd $(pwd)/..
              eval ${function_arrary[3]} ; wait
              cd $(pwd)/..
              eval ${function_arrary[4]} ; wait
              cd $(pwd)/..
            #  eval ${function_arrary[5]} ; wait
              eval ${function_arrary[6]}
          else
              eval ${function_arrary[${ANSWER-$G}]} ; wait
              cd $(pwd)/..
              eval ${function_arrary[${ANSWER2-$G}]} ; wait
              cd $(pwd)/..
              eval ${function_arrary[${ANSWER3-$G}]} ; wait
              cd $(pwd)/..
              eval ${function_arrary[${ANSWER4-$G}]} ; wait
              cd $(pwd)/..
              eval ${function_arrary[${ANSWER5-$G}]} ; wait
              cd $(pwd)/..
              eval ${function_arrary[${ANSWER6-$G}]} ; wait
              cd $(pwd)/..
fi
echo " have a nice day ;)"
