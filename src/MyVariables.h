#ifndef METAMINIMAC_MYVARIABLES_H
#define METAMINIMAC_MYVARIABLES_H
#include <iostream>
#include <stdio.h>
#include "StringBasics.h"
using namespace std;

class myUserVariables
{
    public:
        String inputFiles;
        String info;
        String outfile;
        String method;
        int window,overlap;
        double switchLimit;
        String FileDelimiter;
        double chunkLength;
        int PrintBuffer;
        bool infoDetails;
        bool  gzip, nobgzip;
        String formatString;
        bool GT,DS,GP,HDS,SD;
        char* MyCommandLine;
        string CommandLine;
        String formatStringForVCF;

    myUserVariables()
    {
        FileDelimiter=":";
        inputFiles = "";
        info = "";
        outfile = "MetaMinimac.Output";
        method="B";
        formatString = "GT,DS";
        GT=false;
        DS=false;
        GP=false;
        HDS=false;
        SD=false;
        infoDetails=true;
        formatStringForVCF="";
        window=5000000;
        overlap=100000;
        switchLimit=0.1;
        gzip = true;
        nobgzip = false;
        chunkLength=10;
        PrintBuffer = 100000000;
    };
    bool CheckValidity()
    {

        string formatPiece,formatTemp=formatString.c_str();
        char *end_str1;

        for(char * pch = strtok_r ((char*)formatTemp.c_str(),",", &end_str1);
            pch!=NULL;
            pch = strtok_r (NULL, ",", &end_str1))
        {

            formatPiece=(string)pch;

            if(formatPiece.compare("GT")==0)
            {
                GT=true;
            }
            else if(formatPiece.compare("DS")==0)
            {
                DS=true;
            }
            else if(formatPiece.compare("GP")==0)
            {
                GP=true;
            }
            else if(formatPiece.compare("HDS")==0)
            {
                HDS=true;
            }
            else if(formatPiece.compare("SD")==0)
            {
                SD=true;
            }
            else
            {
                cout << " ERROR !!! \n Cannot identify handle for -f [--format] parameter : "<<formatPiece<<endl;
                cout << " Available handles GT, DS, HDS and GP (for genotype, dosage, haplotype dosage and posterior probability). \n\n";
                cout<<" Program Exiting ..."<<endl<<endl;
                return false;
            }
        }

        bool colonIndex=false;
        if(GT)
        {
            formatStringForVCF+="GT";
            colonIndex=true;
        }
        if(DS)
        {
            formatStringForVCF+= (colonIndex?":DS":"DS");
            colonIndex=true;
        }
        if(HDS)
        {
            formatStringForVCF+= (colonIndex?":HDS":"HDS");
            colonIndex=true;
        }
        if(GP)
        {
            formatStringForVCF+= (colonIndex?":GP":"GP");
            colonIndex=true;
        }
        if(SD)
        {
            formatStringForVCF+= (colonIndex?":SD":"SD");
            colonIndex=true;
        }


        if(nobgzip)
            gzip=false;


        if (inputFiles == "")
        {
            cout<< " Missing -i [--input], a required parameter.\n\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        if(PrintBuffer<=100)
        {
            cout << " ERROR !!! \n Invalid input for -b [--buffer] = "<<PrintBuffer<<"\n";;
            cout << " Buffer for writing output files should be at least 1,000 characters long !!! \n\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<<  " Program Exiting ..."<<endl<<endl;
            return false;
        }

        if(window<=0)
        {
            cout<< " Invalid input for parameter -w [--window] : "<<window<<endl;
            cout<< " --window  can only take Positive Integers ...\n";
            cout<< " Try -h [--help] for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        if(overlap<=0)
        {
            cout<< " Invalid input for parameter -o [--overlap] : "<<overlap<<endl;
            cout<< " --overlap  can only take Positive Integers ...\n";
            cout<< " Try -h [-h [--help]] for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        return true;
    };


    void CreateCommandLine(int argc, char ** argv)
    {
        int len = 0;

        for (int i=0; i<argc; i++)
            len += strlen(argv[i]) + 1;



        char MyCommandLine[len];
        strcpy(MyCommandLine,argv[0]);

        for (int i=1; i<argc; i++)
        {
            strcat(MyCommandLine, " ");
            strcat(MyCommandLine, argv[i]);
        }
        CommandLine=MyCommandLine;
    }

};




#endif //METAMINIMAC_MYVARIABLES_H
