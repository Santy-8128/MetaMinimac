
#include "HaplotypeSet.h"
#include "assert.h"






void HaplotypeSet::PrintFlankLength()
{

    cout<<"  Writing Flanking Length summary to file : "<<InfilePrefix + ".flankLength" + (gzip ? ".gz" : "")<<endl<<endl;
    cout<<"  Writing Flanking Fraction summary to file : "<<InfilePrefix + ".flankFraction" + (gzip ? ".gz" : "")<<endl<<endl;
    IFILE temp=ifopen(InfilePrefix + ".flankLength" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);
    IFILE temp2=ifopen(InfilePrefix + ".flankFraction" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

    for(int i=0;i<numHaplotypes;i++)
    {
        ifprintf(temp,"%s\tHAPLO%d",individualName[i/2].c_str(),i%2+1);
        ifprintf(temp2,"%s\tHAPLO%d",individualName[i/2].c_str(),i%2+1);
        for(int j=0;j<numMarkers;j++)
        {
            ifprintf(temp,"\t%d",FlankLength[i][j]);
            ifprintf(temp2,"\t%0.5f",FlankFrac[i][j]);
        }
        ifprintf(temp,"\n");
        ifprintf(temp2,"\n");
    }
    ifclose(temp);
    ifclose(temp2);
}




void HaplotypeSet::PrintOutputVCF()
{

    HaplotypeSet DosageForVcfPartial;
    IFILE vcfdosepartial=NULL;


    vcfdosepartial = ifopen(outFile + ".dose.vcf" + (gzip ? ".gz" : ""), "wb", gzip ?InputFile::BGZF:InputFile::UNCOMPRESSED);

    ifprintf(vcfdosepartial,"##fileformat=VCFv4.1\n");
    time_t t = time(0);
    struct tm * now = localtime( & t );
    ifprintf(vcfdosepartial,"##filedate=%d.%d.%d\n",(now->tm_year + 1900),(now->tm_mon + 1) ,now->tm_mday);
    ifprintf(vcfdosepartial,"##source=MetaMinimac\n");
    ifprintf(vcfdosepartial,"##contig=<ID=%s>\n",finChromosome.c_str());


    ifprintf(vcfdosepartial,"##FORMAT=<ID=HDS,Number=2,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n");

ifprintf(vcfdosepartial,"##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">\n");


    ifprintf(vcfdosepartial,"##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Estimated Alternate Allele Frequency\">\n");
    ifprintf(vcfdosepartial,"##INFO=<ID=R2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy\">\n");
    ifprintf(vcfdosepartial,"##INFO=<ID=ER2,Number=1,Type=Float,Description=\"Empirical (Leave-One-Out) R-square (available only for genotyped variants)\">\n");
    ifprintf(vcfdosepartial,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");


    for(int hapId=0;hapId<(int)individualName.size();hapId++)
    {
        ifprintf(vcfdosepartial,"\t%s",individualName[hapId].c_str());
    }

    ifprintf(vcfdosepartial,"\n");

    double p1,p2,p3,x,y;
    for (int i = 0; i < numMarkers; i++)
    {

        ifprintf(vcfdosepartial,"%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\tHDS:GP",
        VariantList[i].chr.c_str(),VariantList[i].bp,
        VariantList[i].name.c_str(),VariantList[i].refAlleleString.c_str(),
        VariantList[i].altAlleleString.c_str());


//   cout<<" CHECK = "<<HapDosage[1][0]<<endl;

        for(int hapId=0;hapId<(int)HapDosage.size();hapId+=2)
            {


                    x=HapDosage[hapId][i];
                    y=HapDosage[hapId+1][i];
                    p3=x*y;
                    p2=x*(1-y)+y*(1-x);
                    p1=(1-x)*(1-y);

//                    ifprintf(vcfdosepartial,"\t%.3f|%.3f",x,y);
                    ifprintf(vcfdosepartial,"\t%.3f|%.3f:%.3f,%.3f,%.3f",x,y,p1,p2,p3);
            }

        ifprintf(vcfdosepartial,"\n");
    }


    ifclose(vcfdosepartial);
}



bool HaplotypeSet::printInfoError(int RowNo,string filename)
{
    cout<<"\n Info file does NOT have 14 columns at Row : "<<RowNo<<endl;
    cout<<" Please verify the following file : "<<filename<<endl;
    return false;

}




bool HaplotypeSet::InitializeDosageChunkData(VcfFileReader &ThisVcfFileStream,
                                             VcfRecord &ThisRecord)
{

    if(doesExistFile((InfilePrefix+".dose.vcf").c_str()))
        InfileName=(InfilePrefix+".dose.vcf").c_str();
    else if(doesExistFile((InfilePrefix+".dose.vcf.gz").c_str()))
        InfileName=(InfilePrefix+".dose.vcf.gz").c_str();
    else
    {
        cout<<"\n No VCF file found ("<<InfilePrefix<<".dose.vcf or "<<InfilePrefix+".dose.vcf.gz) "<<endl;
        cout<<" Please check input file prefix ["<< InfilePrefix <<"] properly ... "<<endl;
        return false;
    }

    VcfHeader TempHeader;
    ThisVcfFileStream.open(InfileName, TempHeader);
    ThisVcfFileStream.setSiteOnly(false);
    ThisVcfFileStream.readRecord(ThisRecord);
    TotalNoReadRecods=0;
    LooDosage.resize(numHaplotypes);
    TypedGT.resize(numHaplotypes);


    return true;
}



bool HaplotypeSet::LoadHapDoseVariant(VcfRecordGenotype &ThisGenotype,int &numReadRecords)
{
    int Counter=0;
    for (int i = 0; i<(numSamples); i++)
    {
        string temp=*ThisGenotype.getString("LDS",i);
        char *end_str;

        if(SampleNoHaplotypes[i]==2) {
            char *pch = strtok_r((char *) temp.c_str(), "|", &end_str);
            HapDosage[numReadRecords][Counter++] = atof(pch);

            pch = strtok_r(NULL, "\t", &end_str);
            HapDosage[numReadRecords][Counter++] = atof(pch);
        }
        else
        {
            HapDosage[numReadRecords][Counter++] = atof(temp.c_str());
        }

    }
}


bool HaplotypeSet::LoadHapDoseVariant(VcfRecordGenotype &ThisGenotype)
{
    for (int i = 0; i<(numSamples); i++)
    {
        string temp=*ThisGenotype.getString("HDS",i);
        char *end_str;

        if(SampleNoHaplotypes[i]==2) {
            char *pch = strtok_r((char *) temp.c_str(), ",", &end_str);
            CurrentHapDosage[2*i] = atof(pch);

            pch = strtok_r(NULL, "\t", &end_str);
            CurrentHapDosage[2*i+1] = atof(pch);
        }
        else
        {
            CurrentHapDosage[2*i] = atof(temp.c_str());
        }

    }
}





bool HaplotypeSet::LoadLooVariant(VcfRecordGenotype &ThisGenotype,int loonumReadRecords)
{
    int looCounter=0;
    for (int i = 0; i<(numSamples); i++)
    {
        string temp=*ThisGenotype.getString("LDS",i);
        if(SampleNoHaplotypes[i]==2)
        {
            char *end_str;
            char *pch = strtok_r((char *) temp.c_str(), "|", &end_str);
            LooDosage[2*i][loonumReadRecords] = atof(pch);
            pch = strtok_r(NULL, "\t", &end_str);
            LooDosage[2*i + 1][loonumReadRecords] = atof(pch);


            temp = *ThisGenotype.getString("GT", i);
            char *end_str1;
            char *pch1 = strtok_r((char *) temp.c_str(), "|", &end_str1);
            TypedGT[2*i][loonumReadRecords] = atof(pch1);
            pch1 = strtok_r(NULL, "\t", &end_str1);
            TypedGT[2*i][loonumReadRecords] = atof(pch1);
        }
        else
        {
            LooDosage[2*i][loonumReadRecords] = atof(temp.c_str());
            temp = *ThisGenotype.getString("GT", i);
            TypedGT[2*i][loonumReadRecords] = atof(temp.c_str());

        }
    }
}






bool HaplotypeSet::LoadDosageChunkData(ThisChunk &MyChunk,VcfFileReader &ThisVcfFileStream,
                                             VcfRecord &ThisRecord)
{

    int EndPos=MyChunk.End;
    int numReadRecords=0;
    int loonumReadRecords=0;
    TotalNoReadRecodsBegin=TotalNoReadRecods;
    HapDosage.resize(MyChunk.NoVariants);
    for(int i=0;i<MyChunk.NoVariants;i++)
    {
        HapDosage[i].resize(numHaplotypes);
    }

    for(int i=0;i<numHaplotypes;i++)
    {
        LooDosage[i].resize(MyChunk.NoGenoAllStudies);
        TypedGT[i].resize(MyChunk.NoGenoAllStudies);
    }

    if(!ThisVcfFileStream.isEOF())
    {
        do
        {

            if(ThisRecord.get1BasedPosition()>EndPos)
                break;

            VcfRecordGenotype &ThisGenotype=ThisRecord.getGenotypeInfo();


            LoadHapDoseVariant(ThisGenotype,numReadRecords);


            if(InterAllTypedSitesIndicator[TotalNoReadRecods])
            {
                LoadLooVariant(ThisGenotype,loonumReadRecords);
            }


            numReadRecords++;
            TotalNoReadRecods++;

        }while(ThisVcfFileStream.readRecord(ThisRecord));
    }


    cout<<" File = "<<InfileName<<endl;
    cout<<"   Number of Markers in Chunk = "<<numReadRecords<<endl;
    return true;

}




bool HaplotypeSet::LoadDosageData(string prefix)
{

	VcfFileReader inFile;
	VcfHeader header;
	VcfRecord record;

    string filename;

    if(doesExistFile(prefix+".dose.vcf"))
        filename=prefix+".dose.vcf";
    else if(doesExistFile(prefix+".dose.vcf.gz"))
        filename=prefix+".dose.vcf.gz";
    else
    {
        cout<<"\n No VCF file found ("<<prefix<<".dose.vcf or "<<prefix+".dose.vcf.gz) "<<endl;
        cout<<" Please check input file prefix ["<< prefix <<"] properly ... "<<endl;
//        cout<<" Program aborting ... "<<endl;
        return false;
    }


    int numSamples = 0,markerCount=0,TypedMarkerCount=0, totmarkerCount=0,NumGP=0,NumDS=0,NumGT=0;
    if (!inFile.open(filename.c_str(), header))
	{
		cout << "\n Program could NOT open file : " << filename << endl<<endl;
		return false;
	}

    inFile.setSiteOnly(false);
	numSamples = header.getNumSamples();
    char * pch_split,* pch_split_split,* pch_split3;
    char * pch;
    char *end_str1,*end_str3;

    bool GT,DS,GP,LDS,GT0;
    int doseIndex,gpIndex,gtIndex,ldsIndex,gt0Index;


    HapDosage.resize(2*numSamples);
    LooDosage.resize(2*numSamples);
    GTDosage.resize(2*numSamples);
    for(int i=0;i<(2*numSamples);i++)
    {
        HapDosage[i].resize(numMarkers);
        LooDosage[i].resize(NoTypedSites);
        GTDosage[i].resize(NoTypedSites);
    }




    int part=1;

    string name,refAlleleString,altAlleleString;
    int bp;
    string chr;
    int factor=10000;

    IFILE DoseRead = ifopen(filename.c_str(), "r");
    string line;

    std::cout << "\n Reading Dosage Data from VCF File        = "<<filename<<endl<<endl ;

    if(DoseRead)
    {
        bool Header=true;

        while(Header)
        {
            line.clear();
            DoseRead->readLine(line);
            if(line.substr(0,1).compare("#")==0)
                Header=true;
            else
                Header=false;
        }




        markerCount=0;
        TypedMarkerCount=0;

//            cout<<line<<endl;
//            cout<<markerCount<<" "<<BufferSize<<endl;
        while(markerCount<numMarkers && line.compare("")!=0)
        {

//                abort();

            if(totmarkerCount>=numMarkers)
            {
                std::cout << "\n ERROR !!! Number of markers found in info file = "<<numMarkers ;
                std::cout << "\n           More than "<<numMarkers <<" markers are found in dosage VCF file ... " ;
                std::cout << "\n           Please use info file from same imputation run !!! " << endl;
                return false;
            }


//                cout<<totmarkerCount<<endl;

            char *end_str9;
            pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);
            if(pch==NULL)
            {
                cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
                cout<<" Please verify the following file : "<<filename<<endl;
                return false;
            }
            else
                chr=pch;

            if(finChromosome=="")
                finChromosome=chr;
            else
            {

                if(chr!=finChromosome)
                {
                    cout << "\n Error !!! Input VCF File ["<<filename<<"] contains multiple chromosomes : "<<chr<<", "<<finChromosome<<", ... "<<endl;
                    cout << " Please use VCF file with single chromosome !!! "<<endl;
                    return false;


                }

            }

            pch = strtok_r (NULL, "\t", &end_str1);

            if(pch==NULL)
            {
                cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
                cout<<" Please verify the following file : "<<filename<<endl;
                return false;
            }
            else
                bp=atoi(pch);

            pch = strtok_r (NULL, "\t", &end_str1);
            if(pch==NULL)
            {
                cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
                cout<<" Please verify the following file : "<<filename<<endl;
                return false;
            }
            else
                name=pch;

            if(VariantList[totmarkerCount].name!=name)
            {
                std::cout << "\n ERROR !!! Mismatching variant name between VCF and INFO file !!!" ;
                std::cout << "\n           Marker #"<<totmarkerCount+1 <<" ["<<name<<"] in VCF file does NOT match with"<<
                " Marker #"<<totmarkerCount+1 <<" ["<<VariantList[totmarkerCount].name <<"] in INFO file" ;
                std::cout << "\n           Please use info file from same imputation run !!! " << endl;
                return false;
            }


            variant tempVariant(name,chr,bp);

            pch = strtok_r (NULL, "\t", &end_str1);
            pch = strtok_r (NULL, "\t", &end_str1);

            VariantList[totmarkerCount].assignValues(name,chr,bp);

            pch = strtok_r (NULL, "\t", &end_str1);
            pch = strtok_r (NULL, "\t", &end_str1);
            pch = strtok_r (NULL, "\t", &end_str1);


            pch = strtok_r (NULL, "\t", &end_str1);

//cout<<" FIRSTS  = "<<pch<<endl;

            pch_split = strtok_r (pch,":", &end_str9);


//cout<<" FIRST  = "<<pch<<endl;


            DS=false;
            GP=false;
            GT=false;
            LDS=false;
            GT0=false;

            int index=0;



            while (pch_split != NULL)
            {
                if(strcmp(pch_split,"DS")==0)
                {
                    DS=true;
                    doseIndex=index;
                }
                if(strcmp(pch_split,"GP")==0)
                {
                    GP=true;
                    gpIndex=index;
                }
                if(strcmp(pch_split,"GT")==0)
                {
                    GT=true;
                    gtIndex=index;
                }


                if(TypedSitesIndicator[markerCount])
                {
                     if(strcmp(pch_split,"LDS")==0)
                    {
                        LDS=true;
                        ldsIndex=index;
                    }
                     if(strcmp(pch_split,"GT0")==0)
                    {
                        GT0=true;
                        gt0Index=index;
                    }

                }



                index++;
                pch_split = strtok_r (NULL, ":", &end_str9);
            }



            if(!DS && !GP && !GT)
            {
                cout<<"\n Input VCF file has NONE of the following fields : GT, DS or GP "<<endl;
                cout<<" VCF file must have at least one of these fields to work !!! "<<endl;
                cout<<" Please check the following file : "<<filename<<endl;
                return false;
            }

            int indCount=0;
//
//              if(markerCount==0)
//                        cout<<" LINE = "<<line<<endl;


//cout<<" DERP = "<<pch<<endl;
            pch = strtok_r (NULL, "\t", &end_str1);
//            if(markerCount==20683)
//                        cout<<" UMMAAA = "<<pch<<endl;


            while (pch != NULL)
            {

                char *end_str2;
//                                  if(markerCount==20683)
//                        cout<<" UMMAAA = "<<pch<<endl;

                pch_split = strtok_r (pch,":", &end_str2);
                index=0;


                while (pch_split != NULL)
                {
//
//                    if(markerCount==20683)
//                        cout<<" UMM = "<<pch_split<<endl;

                    if(DS)
                    {
                        if(doseIndex==index)
                            {

                                char *end_str4;
                                pch_split_split = strtok_r (pch_split,"|", &end_str4);


//                cout<<" Ummm = "<<markerCount<<" "<<indCount<<endl;
                                assert(markerCount<numMarkers);
                                assert(indCount<(2*numSamples));

                                HapDosage[2*indCount][markerCount]=atof(pch_split_split);
//                                   if(HapDosage[0][0]!=0.002)
//                                        abort();

                                pch_split_split = strtok_r (NULL,"|", &end_str4);

                                if(pch_split_split==NULL)
                                {
                                    cout<<"\n ERROR: Input VCF file must have haplotype dosages for each sample !!! "<<endl;
                                    std::cout << " Marker #"<<totmarkerCount+1 ;
                                    cout<<" ["<<name<<"] for Indiv "<<indCount<<" in VCF file does NOT have haplotype dosages !!! "<<endl;
                                    cout<<" Please check the following file : "<<filename<<endl;
                                    abort();
                                    return false;
                                }

                                HapDosage[2*indCount+1][markerCount]=atof(pch_split_split);
                                NumDS++;

//                                cout<<HapDosage[2*indCount+1][markerCount]<<endl;
//                                   if(HapDosage[0][0]!=0.002)
//                                        abort();

                            }
                    }
                    else if(GT)
                    {
                        if(gtIndex==index)
                        {
                            NumGT++;
                             end_str3=NULL;
                             pch_split3 = strtok_r (pch_split,"|",
                                                    &end_str3);

                            HapDosage[2*indCount][markerCount]=atoi(pch_split3);
                            pch_split3 = strtok_r (NULL,"|", &end_str3);

                            if(pch_split3==NULL)
                            {
                                cout<<"\n ERROR: Input VCF file must have phased genotype (|) for each sample !!! "<<endl;
                                std::cout << " Marker #"<<totmarkerCount+1 ;
                                cout<<" ["<<name<<"] in VCF file does NOT have phased genotype !!! "<<endl;
                                cout<<" Please check the following file : "<<filename<<endl;
                                return false;
                            }

                            HapDosage[2*indCount+1][markerCount]=atoi(pch_split3);

                        }
                    }
                    else if(GP)
                    {

                        cout<<"\n ERROR: Input VCF file must have phased dosage (HDS) or genotype (GT) for each sample !!! "<<endl;
                        std::cout << " Marker #"<<totmarkerCount+1 ;
                        cout<<" ["<<name<<"] in VCF file only has Genotype Probability (GP) !!! "<<endl;
                        cout<<" Please check the following file : "<<filename<<endl;
                        return false;

                    }


                    if(TypedSitesIndicator[markerCount])
                    {
                        if(LDS)
                        {
                            if(ldsIndex==index)
                            {

//                                cout<<pch_split<<endl;
                                char *end_str4;
                                char *pch_split_split2 = strtok_r (pch_split,"|", &end_str4);

                                assert(TypedSites[TypedMarkerCount]<numMarkers);
                                assert(indCount<(2*numSamples));

                                LooDosage[2*indCount][TypedMarkerCount]=atof(pch_split_split2);

                                pch_split_split2 = strtok_r (NULL,"|", &end_str4);

                                if(pch_split_split2==NULL)
                                {
                                    cout<<"\n ERROR: Input VCF file must have LOO dosages for each sample !!! "<<endl;
                                    std::cout << " Marker #"<<totmarkerCount+1 ;
                                    cout<<" ["<<name<<"] in VCF file does NOT have haplotype dosages !!! "<<endl;
                                    cout<<" Please check the following file : "<<filename<<endl;
                                    return false;
                                }

                                LooDosage[2*indCount+1][TypedMarkerCount]=atof(pch_split_split2);
                            }

                        }



                        if(GT0)
                        {
                            if(gt0Index==index)
                            {
                                char *end_str5;
                                char *pch_split_split3 =strtok_r (pch_split,"|", &end_str5);

                                GTDosage[2*indCount][TypedMarkerCount]=atof(pch_split_split3);
                                pch_split_split3 = strtok_r (NULL,"|", &end_str5);
//
                                if(pch_split_split3==NULL)
                                {
                                    cout<<"\n ERROR: Input VCF file must have phased Typed genotype (|) for each sample !!! "<<endl;
                                    std::cout << " Marker #"<<totmarkerCount+1 ;
                                    cout<<" ["<<name<<"] in VCF file does NOT have phased genotype !!! "<<endl;
                                    cout<<" Please check the following file : "<<filename<<endl;
                                    return false;
                                }

                                GTDosage[2*indCount+1][TypedMarkerCount]=atof(pch_split_split3);
                            }
                        }




                    }




                    pch_split = strtok_r (NULL,":", &end_str2);
                    index++;
                }

                pch = strtok_r (NULL, "\t", &end_str1);
                indCount++;
            }



            if(TypedSitesIndicator[markerCount])
            {
                TypedMarkerCount++;

            }
            markerCount++;
            totmarkerCount++;
            line.clear();
            DoseRead->readLine(line);
        }


    }
    else
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }


    if(totmarkerCount<numMarkers)
    {
        std::cout << "\n\n ERROR !!! Number of markers found in info file = "<<numMarkers ;
        std::cout << "\n           Less than "<<numMarkers <<" markers were read from dosage VCF file ... " ;
        std::cout << "\n           Please use info file from same imputation run !!! " << endl;
        return false;
    }
//	std::cout << "\n Number of Variants imported from GP (Genotype Prob)      : " << NumGP/numSamples;
//	std::cout << "\n Number of Variants imported from DS (Dosage)             : " << NumDS/numSamples;
//	std::cout << "\n Number of Variants imported from GT (Genotype)           : " << NumGT/numSamples<<endl<<endl;


    if((NumGP+NumDS)==0)
    {

        std::cout << "\n WARNING !!! All Dosage values imported from GT values " ;
        std::cout << "\n             No GP or DS values were found in VCF file ! " ;
    }


//abort();

    cout<<"  Dosage Data successfully read ... "<<endl;
    return true;

}


bool HaplotypeSet::CheckSuffixFile(string prefix, const char* suffix, string &FinalName)
{
    if(doesExistFile(prefix+"."+suffix+".vcf"))
        FinalName=prefix+"."+suffix+".vcf";
    else if(doesExistFile(prefix+"."+suffix+".vcf.gz"))
        FinalName=prefix+"."+suffix+".vcf.gz";
    else
    {
        cout<<"\n No VCF file found ("<<prefix<<"."<<suffix<<"dose.vcf or "<<prefix+"."<<suffix<<".vcf.gz) "<<endl;
        cout<<" Please check input file prefix ["<< prefix <<"] properly ... "<<endl;
        return false;
    }
    return true;
}



bool HaplotypeSet::CheckSampleConsistency(int tempNoSamples,
                                          vector<string> &tempindividualName,
                                          vector<int> tempSampleNoHaplotypes,
                                          string File1, string File2)
{
    if(tempNoSamples!=numSamples)
    {
        cout<<"\n ERROR !!! "<<endl;
        cout<< " "<<File1<<" has "<<tempNoSamples<<" samples, but "<< File2<<" has "<<numSamples<<" samples !!!"<<endl;
        return false;
    }
    for(int i=0; i<numSamples; i++)
    {
        if(tempindividualName[i]!=individualName[i])
        {
            cout<<"\n  ERROR !!! "<<endl;
            cout<< " "<<File1<<" and "<< File2<<" have different samples orders !!! "<<endl;
            return false;
        }
        if(tempSampleNoHaplotypes[i]!=SampleNoHaplotypes[i])
        {
            cout<<"\n  ERROR !!! "<<endl;
            cout<< " "<<DoseFileName<<" and "<< EmpDoseFileName<<" have different ploidy for sample ID ["<< individualName[i]<<"]  !!! "<<endl;
            return false;
        }
    }
    return true;
}




void HaplotypeSet::LoadEmpVariantList()
{
    LoadVariantList(EmpDoseFileName);
    TypedVariantList=VariantList;
    noTypedMarkers=numMarkers;
}



void HaplotypeSet::LoadVariantList(string FileName)
{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    inFile.open(FileName.c_str(), header);
    inFile.setSiteOnly(true);
    VariantList.clear();

    int bp,numReadRecords=0;
    string cno,id,refAllele,altAllele,prevID="",currID;

    while (inFile.readRecord(record))
	{
		++numReadRecords;

        variant tempVariant;
        tempVariant.chr=record.getChromStr();
        tempVariant.bp=record.get1BasedPosition();
        tempVariant.name=record.getIDStr();
        tempVariant.altAlleleString = record.getAltStr();
        tempVariant.refAlleleString = record.getRefStr();
        VariantList.push_back(tempVariant);
	}

    numMarkers=VariantList.size();
    inFile.close();
}

void HaplotypeSet::SortCommonGenotypeList(std::unordered_set<string> &CommonGenotypeVariantNameList,
                                          vector<string> &SortedCommonGenoList,
                                          vector<variant> &CommonTypedVariantList)

{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    VcfRecordGenotype *recordGeno;
    inFile.open(EmpDoseFileName.c_str(), header);
    inFile.setSiteOnly(false);
    SortedCommonGenoList.clear();
    int bp,numReadRecords=0;
    string cno,id,refAllele,altAllele,prevID="",currID;

    CommonTypedVariantList.clear();
    CommonTypedVariantList.resize(CommonGenotypeVariantNameList.size());
    LooDosage.clear();
    TypedGT.clear();
    CurrentHapDosage.resize(numHaplotypes);
    LooDosage.resize(numHaplotypes);
    TypedGT.resize(numHaplotypes);
    for(int i=0; i<numHaplotypes; i++)
    {
        LooDosage[i].resize(CommonGenotypeVariantNameList.size());
        TypedGT[i].resize(CommonGenotypeVariantNameList.size());
    }

    int numComRecord = 0;
    while (inFile.readRecord(record))
    {
       ++numReadRecords;

        if(CommonGenotypeVariantNameList.find(record.getIDStr())!=CommonGenotypeVariantNameList.end())
        {
            variant &tempVariant=CommonTypedVariantList[numComRecord];
            tempVariant.chr=record.getChromStr();
            tempVariant.bp=record.get1BasedPosition();
            tempVariant.name=record.getIDStr();
            tempVariant.altAlleleString = record.getAltStr();
            tempVariant.refAlleleString = record.getRefStr();
            VariantList.push_back(tempVariant);

            recordGeno=&record.getGenotypeInfo();
            LoadLooVariant(*recordGeno, numComRecord);

            SortedCommonGenoList.push_back(record.getIDStr());

            numComRecord++;
        }
    }

    if(SortedCommonGenoList.size()!=CommonGenotypeVariantNameList.size())
    {
        cout<<" ERROR CODE 4219: Please contact author with this code to help with bug fixing ..."<<endl;
        abort();
    }
    inFile.close();
}

void HaplotypeSet::ReadBasedOnSortCommonGenotypeList(vector<string> &SortedCommonGenoList)

{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    VcfRecordGenotype *recordGeno;
    inFile.open(EmpDoseFileName.c_str(), header);
    inFile.setSiteOnly(false);
    int bp,numReadRecords=0;
    string cno,id,refAllele,altAllele,prevID="",currID;

    LooDosage.clear();
    TypedGT.clear();
    CurrentHapDosage.resize(numHaplotypes);
    LooDosage.resize(numHaplotypes);
    TypedGT.resize(numHaplotypes);
    for(int i=0; i<numHaplotypes; i++)
    {
        LooDosage[i].resize(SortedCommonGenoList.size());
        TypedGT[i].resize(SortedCommonGenoList.size());
    }
    int SortIndex = 0;
    int numComRecord = 0;
    while (inFile.readRecord(record))
    {
        ++numReadRecords;

        if(SortedCommonGenoList[SortIndex]==record.getIDStr())
        {
            recordGeno=&record.getGenotypeInfo();
            LoadLooVariant(*recordGeno, numComRecord);
            numComRecord++;
            SortIndex++;
        }
    }
    if(SortedCommonGenoList.size()!=SortIndex)
    {
        cout<<" ERROR CODE 2819: Please contact author with this code to help with bug fixing ..."<<endl;
        abort();
    }

    inFile.close();
}

bool HaplotypeSet::LoadSampleNames(string prefix)
{
    InfilePrefix.Copy(prefix.c_str());
    if(!CheckSuffixFile(prefix,"dose", DoseFileName)) return false;
    if(!CheckSuffixFile(prefix,"empiricalDose", EmpDoseFileName)) return false;

    GetSampleInformation(DoseFileName);
    int tempNoSamples=numSamples;
    vector<string> tempindividualName=individualName;
    vector<int> tempSampleNoHaplotypes=SampleNoHaplotypes;
    GetSampleInformation(EmpDoseFileName);
    if(!CheckSampleConsistency(tempNoSamples,tempindividualName,tempSampleNoHaplotypes,DoseFileName,EmpDoseFileName)) return false;

    return true;
//
//
//
//    inFile.setSiteOnly(true);
//    VariantList.clear();
//
//    int bp,numReadRecords=0;
//    string cno,id,refAllele,altAllele,prevID="",currID;
//
//
//    while (inFile.readRecord(record))
//	{
//
//		++numReadRecords;
//
//
//		cno=record.getChromStr();
//        bp=record.get1BasedPosition();
//        id=record.getIDStr();
//        altAllele = record.getAltStr();
//		refAllele = record.getRefStr();
//
//        variant tempVariant(id,cno,bp);
//        tempVariant.assignRefAlt(refAllele,altAllele);
//
//
//        VcfRecordFilter& ThisFilter=record.getFilter();
//        bool typed=false;
//
//        for(int i=0;i<ThisFilter.getNumFilters();i++)
//        {
//            string Filter=ThisFilter.getString(i);
//
//            if(Filter=="GENOTYPED")
//            {
//                TypedSites.push_back(numReadRecords-1);
//                TypedSitesIndicator.push_back(true);
//                typed=true;
//            }
//        }
//        if(!typed)
//            TypedSitesIndicator.push_back(false);
//
//
//        stringstream strs1,strs2;
//        strs1<<(cno);
//        strs2<<(bp);
//        currID=(string)strs1.str()+"_"+(string)strs2.str()+
//                refAllele+altAllele;
//
//        if(prevID==currID)
//        {
//            cout << "\n Error !!! Duplicate Variant found chr:"
//                <<cno<<":"<<bp<<" with identical REF = "<<refAllele
//                 <<" and ALT = "<<altAllele <<".";
//        altAllele=altAllele+altAllele;
//tempVariant.assignRefAlt(refAllele,altAllele);
//
//
////            cout << " Program Aborting ... "<<endl;
////            return false;
//        }
//        else
//            prevID=currID;
//
//
//        VariantList.push_back(tempVariant);
//
//	}
//    NoTypedSites=(int)TypedSites.size();
//    numMarkers=VariantList.size();
//    inFile.close();
//    InterAllTypedSitesIndicator.resize(numMarkers,false);
//    InterAllTypedSitesReverseMap.clear();
//
//    cout << " Number of Samples in VCF File      : " << numSamples << endl;
//    cout << " Number of Markers in VCF File      : " << numMarkers << endl;
//    cout << " Number of Genotyped Markers        : " << NoTypedSites << endl<<endl;
//
//    return true;

}


bool HaplotypeSet::GetSampleInformation(string filename)
{
    VcfFileReader inFile;
    VcfHeader header;
    VcfRecord record;
    individualName.clear();

    if (!inFile.open(filename.c_str(), header))
    {
        cout << "\n Program could NOT open file : " << filename << endl<<endl;
        return false;
    }
    numSamples = header.getNumSamples();
    if(numSamples==0)
    {
        std::cout << "\n Number of Samples read from VCF File    : " << numSamples << endl;
        std::cout << "\n ERROR !!! "<<endl;
        cout << "\n NO samples found in VCF File !! \n Please Check Input File !!!  "<< endl;
        return false;
    }
    individualName.resize(numSamples);
    CummulativeSampleNoHaplotypes.resize(numSamples);
    SampleNoHaplotypes.resize(numSamples);

    for (int i = 0; i < numSamples; i++)
    {
        individualName[i]=header.getSampleName(i);
    }

    inFile.setSiteOnly(false);
    inFile.readRecord(record);
    int tempHapCount=0;
    for (int i = 0; i<(numSamples); i++)
    {
        if(record.getNumGTs(i)==0)
        {
            std::cout << "\n ERROR !!! \n Empty Value for Individual : " << individualName[i] << " at First Marker  " << endl;
            std::cout << " Most probably a corrupted VCF file. Please check input VCF file !!! " << endl;
            cout << "\n Program Exiting ... \n\n";
            return false;
        }
        else
        {
            CummulativeSampleNoHaplotypes[i]=tempHapCount;
            SampleNoHaplotypes[i]=(record.getNumGTs(i));
            tempHapCount+=SampleNoHaplotypes[i];
        }
    }
    inFile.close();
    numHaplotypes=tempHapCount;

    return true;

}


bool HaplotypeSet::LoadInfoFile(string prefix)
{

    string filename;

    if(doesExistFile(prefix+".info"))
        filename=prefix+".info";
    else if(doesExistFile(prefix+".info.gz"))
        filename=prefix+".info.gz";
    else
    {
        cout<<"\n No Info file <"<<prefix<<".info> or <"<<prefix+".info.gz> found "<<endl;
        cout<<" Please check input file prefix ["<< prefix <<"] properly ... "<<endl;
        return false;
    }


    IFILE ifs = ifopen(filename.c_str(), "r");

    VariantList.clear();
    string name,refAlleleString,altAlleleString;
//    string MajAlleleString,MinAlleleString;
    int bp=100;
    string chr="Z";
    char *pch,*end_str1;

    int Counter=0;

    int RowNo=1;
    TypedSites.clear();
    TypedSitesIndicator.clear();

    string line;
    if(ifs)
    {
        ifs->readLine(line);
        line.clear();
        int Indic=ifs->readLine(line);

        while(Indic==-1 ||  Indic==0)
        {

            if(Indic==-1 && line.length()==0)
            {
                break;
            }


            RowNo++;
            pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);

            if(pch==NULL)
                return printInfoError(RowNo,filename);
            else
                name=pch;

            variant tempVariant(name,chr,bp);

            cout<<" WHY " <<tempVariant.bp<<endl;

            pch = strtok_r (NULL, "\t", &end_str1);

            if(pch==NULL)
                return printInfoError(RowNo,filename);
            else
                refAlleleString=pch;

            pch = strtok_r (NULL, "\t", &end_str1);

            if(pch==NULL)
                return printInfoError(RowNo,filename);
            else
                altAlleleString=pch;

            pch = strtok_r (NULL, "\t", &end_str1);
            pch = strtok_r (NULL, "\t", &end_str1);
            pch = strtok_r (NULL, "\t", &end_str1);
            pch = strtok_r (NULL, "\t", &end_str1);

            if(pch==NULL)
                return printInfoError(RowNo,filename);
            else
                EstRsq.push_back(atof(pch));



            pch = strtok_r (NULL, "\t", &end_str1);

            if(pch==NULL)
                return printInfoError(RowNo,filename);
            else
            {
                if(((string)pch).compare("Genotyped")==0)
                {
                    TypedSites.push_back(RowNo-2);
                    TypedSitesIndicator.push_back(true);
                }
                else
                {
                    TypedSitesIndicator.push_back(false);
                }
            }

            tempVariant.assignRefAlt(refAlleleString,altAlleleString);
            VariantList.push_back(tempVariant);
            line.clear();
            Indic=ifs->readLine(line);

        }
    }
    else
    {
        cout<<" Following File File Not Available : "<<filename<<endl;
        return false;
    }

    NoTypedSites=(int)TypedSites.size();
    numMarkers=VariantList.size();

	std::cout << " Number of Markers in INFO File     : " << numMarkers << endl<<endl;

	ifclose(ifs);
	return true;

}


bool HaplotypeSet::CheckMarkerCount(string prefix)
{

    string filename;

    if(doesExistFile(prefix+".info"))
        filename=prefix+".info";
    else if(doesExistFile(prefix+".info.gz"))
        filename=prefix+".info.gz";
    else
    {
        cout<<"\n No Info file <"<<prefix<<".info> or <"<<prefix+".info.gz> found "<<endl;
        cout<<" Please check input file prefix ["<< prefix <<"] properly ... "<<endl;
        return false;
    }


    IFILE ifs = ifopen(filename.c_str(), "r");
    int RowNo=0;
    string line;
    if(ifs)
    {
        ifs->readLine(line);
        line.clear();
        int Indic=ifs->readLine(line);

        while(Indic==-1 ||  Indic==0)
        {

            if(Indic==-1 && line.length()==0)
            {
                break;
            }

            RowNo++;
            line.clear();
            Indic=ifs->readLine(line);

        }
    }
    else
    {
        cout<<" Following File File Not Available : "<<filename<<endl;
        return false;
    }

    numMarkers=RowNo;

    cout << " File Name : "<<filename<<endl;
	cout << " Number of Markers  : " << numMarkers << endl;

	ifclose(ifs);
	return true;

}


bool HaplotypeSet::UpdateFlankLength(int index,vector<double> tempRecom)
{

    // Use recursion for calculating flank length

    vector<int> &tempFlankLength=FlankLength[index];


    int StartCounter=0,EndCounter=0;

    double CurrentLeftProduct=1.0;
    double CurrentRightProduct=1.0;

    int prevLeftEnd=0,prevRightEnd=0;


    bool Startflag=true;

     // Update first and last markers to be 0
    tempFlankLength[0]=0;
    tempFlankLength[numMarkers-1]=0;

    for(int i=1;i<(numMarkers-1);i++)
    {

        //If prevRightEnd occurs before this marker, re-create it
        if(prevRightEnd==(i-1))
        {
            CurrentRightProduct=1.0-tempRecom[i];
            prevRightEnd++;

        }
        else
        {
            CurrentRightProduct/=(1-tempRecom[i-1]);
        }

        // Update by removing/adding left of the marker

        CurrentLeftProduct*=(1-tempRecom[i-1]);


        // Keep moving right till value is less than cutoff


        assert(prevLeftEnd<tempRecom.size());
        double tempVal=CurrentLeftProduct/(1-tempRecom[prevLeftEnd]);
        while(tempVal<(1-SplitLimit))
        {
            CurrentLeftProduct=tempVal;
            Startflag=false;
            ++prevLeftEnd;
            tempVal=CurrentLeftProduct/(1-tempRecom[prevLeftEnd]);
        }
        int leftLength=(i-prevLeftEnd-1);

        // If beginning of chromosome, prevLeftEnd measures less than cut-off, so leftLength should be increased by one.
        if(Startflag)
        {
            leftLength++;
        }


        // Keep moving right till value is less than cutoff
        while(CurrentRightProduct>=(1-SplitLimit) && prevRightEnd<(numMarkers-1))
        {
            CurrentRightProduct*=(1-tempRecom[++prevRightEnd]);

        }
        int rightLength=(prevRightEnd-i);

//        cout<<"MARKER = "<<" "<<i<<"\t"<<leftLength<<"\t"<<rightLength<<"\t"<<CurrentLeftProduct<<"\t"<<CurrentRightProduct<<endl;

        tempFlankLength[i]=leftLength<rightLength?leftLength:rightLength;

    }


}


bool HaplotypeSet::ProcessRecombination(int index,vector<double> tempRecom)
{



    UpdateFlankLength(index,tempRecom);


}

bool HaplotypeSet::LoadRecomSpectrum(string prefix)
{

    RecomSpec.clear();
    string RecomFile;
    std::string delimiter("\t") ;
    size_t pos = 0;

    RecomSpec.resize(numHaplotypes);
    FlankLength.resize(numHaplotypes);
    FlankFrac.resize(numHaplotypes);
    for(int i=0;i<numHaplotypes;i++)
    {
        RecomSpec[i].resize(numMarkers-1);
        FlankLength[i].resize(numMarkers,0);
        FlankFrac[i].resize(numMarkers,0.0);
    }


    if(doesExistFile(prefix+".metaRecom"))
        RecomFile=prefix+".metaRecom";
    else if(doesExistFile(prefix+".metaRecom.gz"))
        RecomFile=prefix+".metaRecom.gz";
    else
    {
        cout<<"\n No Recombination file <"<<prefix<<".metaRecom> or <"<<prefix+".metaRecom.gz> found "<<endl;
        cout<<" Please check input file prefix ["<< prefix <<"] properly ... "<<endl;
        return false;
    }


    cout<<"\n Reading Recombination Spectrum from File = "<<RecomFile<<endl;

    IFILE Recomifs = ifopen(RecomFile.c_str(), "r");
    int HapCount=0;

    string line;
    if(Recomifs)
    {
        while(Recomifs->readLine(line)==0)
        {
            int ColNo=0;

            vector<double> &tempRecom=RecomSpec[HapCount];

            char *end_str1;
            char *pch;

            pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);

            while (pch != NULL)
            {

                ColNo++;
                if(ColNo>2)
                {
                    assert(HapCount<numHaplotypes);
                    assert((ColNo-3)<numMarkers-1);
                    tempRecom[ColNo-3] = atof(pch);
                }
                pch = strtok_r (NULL, "\t", &end_str1);
            }

            ProcessRecombination(HapCount,tempRecom);

            HapCount++;
            line.clear();
        }

    }
    else
    {
        cout<<"\n Following File File Not Available : "<<RecomFile<<endl;
        return false;
    }

    cout<<"\n  Recombination Spectrum successfully read ... "<<endl<<endl;





    return true;

}





bool HaplotypeSet::GetSummary(string prefix, myUserVariables &ThisVariable)
{

    cout<<  " File Prefix                        : "<<prefix<<endl;
    if(!LoadSampleNames(prefix))
        return false;
    return true;
}



bool HaplotypeSet::ReadInputFiles(string prefix)
{
    HapDosage.clear();

    if(!LoadInfoFile(prefix))
        return false;

    if(!LoadSampleNames(prefix))
        return false;



    if(Method=="C")
        if(!LoadRecomSpectrum(prefix))
            return false;

    if(!LoadDosageData(prefix))
        return false;


    return true;
}






bool HaplotypeSet::FastLoadHaplotypes()
{
    String filename=VcfDose;



    cout<<"\n Detecting Dosage File Type ... "<<endl;

    string FileType=DetectReferenceFileType(filename);

    if(FileType.compare("NA")==0)
    {
        cout<<"\n Following File File Not Available : "<<filename<<endl;
        return false;
    }
    else if(FileType.compare("Invalid")==0)
    {

        cout << "\n Dosage File provided by \"--vcfDose\" must be a VCF file !!! \n";
        cout << " Please check the following file : "<<filename<<endl;
        return false;
    }

    cout<<"\n Format = VCF (Variant Call Format) "<<endl;


    cout<<"\n Reading Info File                       : "<<Info<<endl;

//
//    if(!LoadInfoFile(Info))
//        return false;



    std::cout << "\n Loading Dosage Data from VCF File       : " << filename << endl<<endl;

//    return WriteMachFile();

}



string HaplotypeSet::DetectReferenceFileType(String filename)
{
    IFILE fileStream = ifopen(filename, "r");
    string line;
    if(fileStream)
    {
        fileStream->readLine(line);
        if(line.length()<1)
            return "Invalid";
         string tempString;
        tempString=(line.substr(0,17));

        char temp[tempString.length() + 1];
        std::strcpy(temp,tempString.c_str());
        for (char *iter = temp; *iter != '\0'; ++iter)
        {
           *iter = std::tolower(*iter);
        }

        if(((string)temp).compare("##fileformat=vcfv")==0)
        {
            return "vcf";
        }
        else
            return "Invalid";

    }
    else
    {
        return "NA";
    }
    ifclose(fileStream);

    return "NA";
}



bool HaplotypeSet::doesExistFile(string filename)
{
    IFILE ifs = ifopen(filename.c_str(), "r");
    if (ifs)
    {
        ifclose(ifs);
        return true;
    }
    else
    {
        ifclose(ifs);
        return false;
    }
}





//
//    bool GT,DS,GP,LDS,GT0;
//    int doseIndex,gpIndex,gtIndex,ldsIndex,gt0Index;
//
//
//    HapDosage.resize(2*numSamples);
//    LooDosage.resize(2*numSamples);
//    GTDosage.resize(2*numSamples);
//    for(int i=0;i<(2*numSamples);i++)
//    {
//        HapDosage[i].resize(numMarkers);
//        LooDosage[i].resize(NoTypedSites);
//        GTDosage[i].resize(NoTypedSites);
//    }
//
//
//
//
//    int part=1;
//
//    string name,refAlleleString,altAlleleString;
//    int bp;
//    string chr;
//    int factor=10000;
//
//    IFILE DoseRead = ifopen(filename.c_str(), "r");
//    string line;
//
//    std::cout << "\n Reading Dosage Data from VCF File        = "<<filename<<endl<<endl ;
//
//    if(DoseRead)
//    {
//        bool Header=true;
//
//        while(Header)
//        {
//            line.clear();
//            DoseRead->readLine(line);
//            if(line.substr(0,1).compare("#")==0)
//                Header=true;
//            else
//                Header=false;
//        }
//
//
//
//
//        markerCount=0;
//        TypedMarkerCount=0;
//
////            cout<<line<<endl;
////            cout<<markerCount<<" "<<BufferSize<<endl;
//        while(markerCount<numMarkers && line.compare("")!=0)
//        {
//
////                abort();
//
//            if(totmarkerCount>=numMarkers)
//            {
//                std::cout << "\n ERROR !!! Number of markers found in info file = "<<numMarkers ;
//                std::cout << "\n           More than "<<numMarkers <<" markers are found in dosage VCF file ... " ;
//                std::cout << "\n           Please use info file from same imputation run !!! " << endl;
//                return false;
//            }
//
//
////                cout<<totmarkerCount<<endl;
//
//            char *end_str9;
//            pch = strtok_r ((char*)line.c_str(),"\t", &end_str1);
//            if(pch==NULL)
//            {
//                cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
//                cout<<" Please verify the following file : "<<filename<<endl;
//                return false;
//            }
//            else
//                chr=pch;
//
//            if(finChromosome=="")
//                finChromosome=chr;
//            else
//            {
//
//                if(chr!=finChromosome)
//                {
//                    cout << "\n Error !!! Input VCF File ["<<filename<<"] contains multiple chromosomes : "<<chr<<", "<<finChromosome<<", ... "<<endl;
//                    cout << " Please use VCF file with single chromosome !!! "<<endl;
//                    return false;
//
//
//                }
//
//            }
//
//            pch = strtok_r (NULL, "\t", &end_str1);
//
//            if(pch==NULL)
//            {
//                cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
//                cout<<" Please verify the following file : "<<filename<<endl;
//                return false;
//            }
//            else
//                bp=atoi(pch);
//
//            pch = strtok_r (NULL, "\t", &end_str1);
//            if(pch==NULL)
//            {
//                cout<<"\n Input VCF file has invalid format at Marker No : "<<totmarkerCount+1<<endl;
//                cout<<" Please verify the following file : "<<filename<<endl;
//                return false;
//            }
//            else
//                name=pch;
//
//            if(VariantList[totmarkerCount].name!=name)
//            {
//                std::cout << "\n ERROR !!! Mismatching variant name between VCF and INFO file !!!" ;
//                std::cout << "\n           Marker #"<<totmarkerCount+1 <<" ["<<name<<"] in VCF file does NOT match with"<<
//                " Marker #"<<totmarkerCount+1 <<" ["<<VariantList[totmarkerCount].name <<"] in INFO file" ;
//                std::cout << "\n           Please use info file from same imputation run !!! " << endl;
//                return false;
//            }
//
//
//            variant tempVariant(name,chr,bp);
//
//            pch = strtok_r (NULL, "\t", &end_str1);
//            pch = strtok_r (NULL, "\t", &end_str1);
//
//            VariantList[totmarkerCount].assignValues(name,chr,bp);
//
//            pch = strtok_r (NULL, "\t", &end_str1);
//            pch = strtok_r (NULL, "\t", &end_str1);
//            pch = strtok_r (NULL, "\t", &end_str1);
//
//
//            pch = strtok_r (NULL, "\t", &end_str1);
//
////cout<<" FIRSTS  = "<<pch<<endl;
//
//            pch_split = strtok_r (pch,":", &end_str9);
//
//
////cout<<" FIRST  = "<<pch<<endl;
//
//
//            DS=false;
//            GP=false;
//            GT=false;
//            LDS=false;
//            GT0=false;
//
//            int index=0;
//
//
//
//            while (pch_split != NULL)
//            {
//                if(strcmp(pch_split,"DS")==0)
//                {
//                    DS=true;
//                    doseIndex=index;
//                }
//                if(strcmp(pch_split,"GP")==0)
//                {
//                    GP=true;
//                    gpIndex=index;
//                }
//                if(strcmp(pch_split,"GT")==0)
//                {
//                    GT=true;
//                    gtIndex=index;
//                }
//
//
//                if(TypedSitesIndicator[markerCount])
//                {
//                     if(strcmp(pch_split,"LDS")==0)
//                    {
//                        LDS=true;
//                        ldsIndex=index;
//                    }
//                     if(strcmp(pch_split,"GT0")==0)
//                    {
//                        GT0=true;
//                        gt0Index=index;
//                    }
//
//                }
//
//
//
//                index++;
//                pch_split = strtok_r (NULL, ":", &end_str9);
//            }
//
//
//
//            if(!DS && !GP && !GT)
//            {
//                cout<<"\n Input VCF file has NONE of the following fields : GT, DS or GP "<<endl;
//                cout<<" VCF file must have at least one of these fields to work !!! "<<endl;
//                cout<<" Please check the following file : "<<filename<<endl;
//                return false;
//            }
//
//            int indCount=0;
////
////              if(markerCount==0)
////                        cout<<" LINE = "<<line<<endl;
//
//
////cout<<" DERP = "<<pch<<endl;
//            pch = strtok_r (NULL, "\t", &end_str1);
////            if(markerCount==20683)
////                        cout<<" UMMAAA = "<<pch<<endl;
//
//
//            while (pch != NULL)
//            {
//
//                char *end_str2;
////                                  if(markerCount==20683)
////                        cout<<" UMMAAA = "<<pch<<endl;
//
//                pch_split = strtok_r (pch,":", &end_str2);
//                index=0;
//
//
//                while (pch_split != NULL)
//                {
////
////                    if(markerCount==20683)
////                        cout<<" UMM = "<<pch_split<<endl;
//
//                    if(DS)
//                    {
//                        if(doseIndex==index)
//                            {
//
//                                char *end_str4;
//                                pch_split_split = strtok_r (pch_split,"|", &end_str4);
//
//
////                cout<<" Ummm = "<<markerCount<<" "<<indCount<<endl;
//                                assert(markerCount<numMarkers);
//                                assert(indCount<(2*numSamples));
//
//                                HapDosage[2*indCount][markerCount]=atof(pch_split_split);
////                                   if(HapDosage[0][0]!=0.002)
////                                        abort();
//
//                                pch_split_split = strtok_r (NULL,"|", &end_str4);
//
//                                if(pch_split_split==NULL)
//                                {
//                                    cout<<"\n ERROR: Input VCF file must have haplotype dosages for each sample !!! "<<endl;
//                                    std::cout << " Marker #"<<totmarkerCount+1 ;
//                                    cout<<" ["<<name<<"] for Indiv "<<indCount<<" in VCF file does NOT have haplotype dosages !!! "<<endl;
//                                    cout<<" Please check the following file : "<<filename<<endl;
//                                    abort();
//                                    return false;
//                                }
//
//                                HapDosage[2*indCount+1][markerCount]=atof(pch_split_split);
//                                NumDS++;
//
////                                cout<<HapDosage[2*indCount+1][markerCount]<<endl;
////                                   if(HapDosage[0][0]!=0.002)
////                                        abort();
//
//                            }
//                    }
//                    else if(GT)
//                    {
//                        if(gtIndex==index)
//                        {
//                            NumGT++;
//                             end_str3=NULL;
//                             pch_split3 = strtok_r (pch_split,"|",
//                                                    &end_str3);
//
//                            HapDosage[2*indCount][markerCount]=atoi(pch_split3);
//                            pch_split3 = strtok_r (NULL,"|", &end_str3);
//
//                            if(pch_split3==NULL)
//                            {
//                                cout<<"\n ERROR: Input VCF file must have phased genotype (|) for each sample !!! "<<endl;
//                                std::cout << " Marker #"<<totmarkerCount+1 ;
//                                cout<<" ["<<name<<"] in VCF file does NOT have phased genotype !!! "<<endl;
//                                cout<<" Please check the following file : "<<filename<<endl;
//                                return false;
//                            }
//
//                            HapDosage[2*indCount+1][markerCount]=atoi(pch_split3);
//
//                        }
//                    }
//                    else if(GP)
//                    {
//
//                        cout<<"\n ERROR: Input VCF file must have phased dosage (HDS) or genotype (GT) for each sample !!! "<<endl;
//                        std::cout << " Marker #"<<totmarkerCount+1 ;
//                        cout<<" ["<<name<<"] in VCF file only has Genotype Probability (GP) !!! "<<endl;
//                        cout<<" Please check the following file : "<<filename<<endl;
//                        return false;
//
//                    }
//
//
//                    if(TypedSitesIndicator[markerCount])
//                    {
//                        if(LDS)
//                        {
//                            if(ldsIndex==index)
//                            {
//
////                                cout<<pch_split<<endl;
//                                char *end_str4;
//                                char *pch_split_split2 = strtok_r (pch_split,"|", &end_str4);
//
//                                assert(TypedSites[TypedMarkerCount]<numMarkers);
//                                assert(indCount<(2*numSamples));
//
//                                LooDosage[2*indCount][TypedMarkerCount]=atof(pch_split_split2);
//
//                                pch_split_split2 = strtok_r (NULL,"|", &end_str4);
//
//                                if(pch_split_split2==NULL)
//                                {
//                                    cout<<"\n ERROR: Input VCF file must have LOO dosages for each sample !!! "<<endl;
//                                    std::cout << " Marker #"<<totmarkerCount+1 ;
//                                    cout<<" ["<<name<<"] in VCF file does NOT have haplotype dosages !!! "<<endl;
//                                    cout<<" Please check the following file : "<<filename<<endl;
//                                    return false;
//                                }
//
//                                LooDosage[2*indCount+1][TypedMarkerCount]=atof(pch_split_split2);
//                            }
//
//                        }
//
//
//
//                        if(GT0)
//                        {
//                            if(gt0Index==index)
//                            {
//                                char *end_str5;
//                                char *pch_split_split3 =strtok_r (pch_split,"|", &end_str5);
//
//                                GTDosage[2*indCount][TypedMarkerCount]=atof(pch_split_split3);
//                                pch_split_split3 = strtok_r (NULL,"|", &end_str5);
////
//                                if(pch_split_split3==NULL)
//                                {
//                                    cout<<"\n ERROR: Input VCF file must have phased Typed genotype (|) for each sample !!! "<<endl;
//                                    std::cout << " Marker #"<<totmarkerCount+1 ;
//                                    cout<<" ["<<name<<"] in VCF file does NOT have phased genotype !!! "<<endl;
//                                    cout<<" Please check the following file : "<<filename<<endl;
//                                    return false;
//                                }
//
//                                GTDosage[2*indCount+1][TypedMarkerCount]=atof(pch_split_split3);
//                            }
//                        }
//
//
//
//
//                    }
//
//
//
//
//                    pch_split = strtok_r (NULL,":", &end_str2);
//                    index++;
//                }
//
//                pch = strtok_r (NULL, "\t", &end_str1);
//                indCount++;
//            }
//
//
//
//            if(TypedSitesIndicator[markerCount])
//            {
//                TypedMarkerCount++;
//
//            }
//            markerCount++;
//            totmarkerCount++;
//            line.clear();
//            DoseRead->readLine(line);
//        }
//
//
//    }
//    else
//    {
//        cout<<"\n Following File File Not Available : "<<filename<<endl;
//        return false;
//    }
//
//
//    if(totmarkerCount<numMarkers)
//    {
//        std::cout << "\n\n ERROR !!! Number of markers found in info file = "<<numMarkers ;
//        std::cout << "\n           Less than "<<numMarkers <<" markers were read from dosage VCF file ... " ;
//        std::cout << "\n           Please use info file from same imputation run !!! " << endl;
//        return false;
//    }
////	std::cout << "\n Number of Variants imported from GP (Genotype Prob)      : " << NumGP/numSamples;
////	std::cout << "\n Number of Variants imported from DS (Dosage)             : " << NumDS/numSamples;
////	std::cout << "\n Number of Variants imported from GT (Genotype)           : " << NumGT/numSamples<<endl<<endl;
//
//
//    if((NumGP+NumDS)==0)
//    {
//
//        std::cout << "\n WARNING !!! All Dosage values imported from GT values " ;
//        std::cout << "\n             No GP or DS values were found in VCF file ! " ;
//    }
//
//
////abort();
//
//    cout<<"  Dosage Data successfully read ... "<<endl;
//    return true;
//
//}

