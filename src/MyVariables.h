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
        String format ;
        int window,overlap;
        double switchLimit;
        String FileDelimiter;

        bool  gzip, nobgzip;
        bool   help, params, log;


    myUserVariables()
    {
        FileDelimiter=":";
        inputFiles = "";
        info = "";
        outfile = "MetaMinimac.Output";
        method="B";
        format = "DS";
        window=5000000;
        overlap=100000;
        switchLimit=0.1;
        gzip = true;
        nobgzip = false;
        help = false;
        params = false;
        log=false;

    };
    bool CheckValidity()
    {

        if(nobgzip)
            gzip=false;


        if (inputFiles == "")
        {
            cout<< " Missing \"--inputFiles\", a required parameter.\n\n";
            cout<< " Try \"--help\" for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        if(format!="GP" && format!="DS")
        {
            cout<< " Invalid input for \"--format\" parameter : "<<format<<endl;
            cout<< " Parameter must be equal to \"GP\" or \"DS\" ..."<<endl<<endl;
            cout<< " Try \"--help\" for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        if(window<=0)
        {
            cout<< " Invalid input for parameter \"--window\" : "<<window<<endl;
            cout<< " \"--window\"  can only take Positive Integers ...\n";
            cout<< " Try \"--help\" for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        if(overlap<=0)
        {
            cout<< " Invalid input for parameter \"--overlap\" : "<<overlap<<endl;
            cout<< " \"--overlap\"  can only take Positive Integers ...\n";
            cout<< " Try \"--help\" for usage ...\n\n";
            cout<< " Program Exiting ...\n\n";
            return false;
        }

        return true;
    };


};




#endif //METAMINIMAC_MYVARIABLES_H
