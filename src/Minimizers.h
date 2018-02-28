#ifndef MINIMIZERS_H_INCLUDED
#define MINIMIZERS_H_INCLUDED

#include "MetaMinimac.h"
#include <math.h>
#include <assert.h>

using namespace std;


class LeastSquareError
{
    private:

        String Method;
        int NoMarkers;
        int NoStudies;
        int NoDimensions;
        int SampleID;
        vector<vector<double>* > LooDosage;
        vector<vector<int*> > FlankLength;
        vector<vector<double>* > RecomSum;
        vector<double> *ChipGT;


    public:

        double initialize(int SampleId, vector<HaplotypeSet> &InputData,
                                    MetaMinimac *const ThisStudy,ThisChunk &MyChunk);
        double  operator()(vector<double> x);
};


class LeastSquareWithRecombination
{
    private:

        String Method;
        int NoMarkers;
        int NoStudies;
        int NoDimensions;
        int SampleID;
        vector<vector<double>* > LooDosage;
        vector<vector<int*> > FlankLength;
        vector<vector<double*> > FlankFrac;
        vector<vector<double>* > RecomSum;
        vector<double> *ChipGT;


    public:

        double initialize(int SampleId, vector<HaplotypeSet> &InputData,
                                    MetaMinimac *const ThisStudy,ThisChunk &MyChunk);
        double  operator()(vector<double> x);


double FixedSSQ();

};




class LogOddsModel
{
    private:

        String Method;
        int NoMarkers;
        int NoStudies;
        int NoDimensions;
        int NoParams;
        int SampleID;
        vector<vector<double>* > LooDosage;
        vector<vector<double> > LooDosageVal;
        vector<double> ChipGTVal;
        vector<vector<int*> > FlankLength;
        vector<vector<double*> > FlankFrac;
        vector<vector<double*> > FlankLogOddFrac;


        vector<vector<double>* > RecomSum;
        vector<double> *ChipGT;


    public:

    double metaInitialize(int SampleId, vector<HaplotypeSet> &InputData,
                      MetaMinimac *const ThisStudy,ThisChunk &MyChunk);
    double initialize(int SampleId, vector<HaplotypeSet> &InputData,
                      MetaMinimac *const ThisStudy,ThisChunk &MyChunk);
        double  operator()(vector<double> x);


double FixedSSQ();

};




void logitTransform(vector<double> &From,
                    vector<double> &To)
{

    double sum=1.0;
    int NoDimensions = (int)To.size();
    for(int i=0; i < (NoDimensions-1); i++) sum+=exp(From[i]);
    for(int i=0; i < (NoDimensions-1); i++)  To[i]=exp(From[i])/sum;
    To[NoDimensions-1]=1.0/sum;


    double checkSum=0.0;
    for(int i=0;i<To.size();i++)
        checkSum+=To[i];
    if(checkSum>1.0001)
        abort();


}



void logitTransform(vector<double> From,
                    vector<double> &To,
                    vector< vector<double> > &Covariates,
                    int NoDimensions,
                    String Method)
{

    double sum=0.0;

    if(Method=="A")
    {
        for(int i=0; i < (NoDimensions-1); i++)
        {
            sum+=exp(From[i]);
        }
    }
    else if(Method=="B")
    {
        for(int i=0; i < (NoDimensions-1); i++)
        {
            From[i]=From[i]+( Covariates[0][i] * From[NoDimensions-1]);

            sum+=exp(From[i]);
        }
    }
    else if(Method=="C")
    {
        for(int i=0; i < (NoDimensions-1); i++)
        {
            From[i]=From[i]+( Covariates[0][i] * From[i+NoDimensions-1]);
            sum+=exp(From[i]);
        }
    }
    else if(Method=="E")
    {
        for(int i=0; i < (NoDimensions-1); i++)
        {
            sum+=exp(From[i]+( Covariates[0][i] * From[i+NoDimensions-1])+( Covariates[1][i] * From[i+2*(NoDimensions-1)]));
        }
    }


    sum+=1.0;

    for(int i=0; i < (NoDimensions-1); i++)
    {
        To[i]=exp(From[i])/sum;
    }

    To[NoDimensions-1]=1.0/sum;



    double checkSum=0.0;

    for(int i=0;i<To.size();i++)
        checkSum+=To[i];

//    cout<<" CHECK = "<<checkSum<<endl;

    if(checkSum>1.0001)
        abort();
}



void logTransform(vector<double> &From,vector<double> &To,int NoDimensions)
{

    double p=1.0;
    for(int i=0; i < NoDimensions; i++)
    {
        double logit = 1.0/(1.0+exp(0.0-From[i]));
        To[i]=p*logit;
        p = p*(1.-logit);
    }
    To[NoDimensions]=p;

}


double LeastSquareError::initialize(int SampleId, vector<HaplotypeSet> &InputData,
                                    MetaMinimac *const ThisStudy,ThisChunk &MyChunk)
{

    Method=ThisStudy->Method;
    NoStudies=InputData.size();
    NoDimensions=NoStudies-1;
    NoMarkers=MyChunk.NoGenoAllStudies;
    SampleID=SampleId;



    LooDosage.resize(NoStudies);
    for(int i=0;i<NoStudies;i++)
    {
        LooDosage[i]=&(InputData[i].LooDosage[SampleID]);
    }
    ChipGT=&InputData[0].TypedGT[SampleID];

}



double LeastSquareError::operator()(vector<double> x)
{

    double sum=0.0;
    vector<double> tempVar(NoDimensions+1);

    logTransform(x,tempVar,NoDimensions);



//
//    for(int i=0; i < NoDimensions; ++i)
//    {
//        double logit = 1.0/(1.0+exp(0.0-x[i]));
//        tempVar[i]=p*logit;
//        cout<<tempVar[i]<<"\t";
//        p = p*(1.-logit);
//    }
//    tempVar[NoDimensions]=p;
//    cout<<tempVar[NoDimensions];

    for(int i=0;i<NoMarkers;i++)
    {
        double temp=(*ChipGT)[i];
        for(int j=0;j<=NoDimensions;j++)
        {

            temp-=(tempVar[j]*((*LooDosage[j])[i]));
//              cout<<" THERAPY2 = "<<j<<"\t"<<tempVar[j]<<"\t"<<(*LooDosage[j])[i]<<endl;


//        cout<<" D = "<<temp<<"\t";

        }
//        cout<<"s= "<<temp<<"\t";

        temp=temp*temp;

        sum+=temp;

    }
//  cout<<"\t= "<<sum<<"\n";
  return sum;
}



double LogOddsModel::initialize(int SampleId, vector<HaplotypeSet> &InputData,
                                    MetaMinimac *const ThisStudy,ThisChunk &MyChunk)
{

    Method=ThisStudy->Method;
    NoStudies=InputData.size();

    if(Method=="A")
        {
            NoParams=NoStudies-1;
            NoDimensions=NoStudies;
        }
    else if(Method=="B")
        {
            NoParams=NoStudies;
            NoDimensions=NoStudies;
        }
     else if(Method=="C")
        {
            NoParams=2*(NoStudies-1);
            NoDimensions=NoStudies;
        }

    // 1 extra for the flank-Length parameter



    NoMarkers=MyChunk.NoGenoAllStudies;
    SampleID=SampleId;


    LooDosage.resize(NoStudies);
    for(int i=0;i<NoStudies;i++)
    {
        LooDosage[i]=&(InputData[i].LooDosage[SampleID]);
    }
    ChipGT=&InputData[0].TypedGT[SampleID];


    if(Method!="A")
    {
        FlankLength.resize(NoStudies);
        FlankFrac.resize(NoStudies);
        FlankLogOddFrac.resize(NoStudies);

        RecomSum.resize(NoStudies);

        for(int i=0;i<NoStudies;i++)
        {
            vector<int> &tempMapper=MyChunk.ThisChunkInterAllTypedSitesReverseMap[i];
            vector<int> &FromFlank = InputData[i].FlankLength[SampleID];
            vector<double> &FromFlankFrac = InputData[i].FlankFrac[SampleID];
//            vector<double> &FromFlankLogOddFrac = InputData[i].FlankLogOddFrac[SampleID];


            vector<int*> &ToFlank = FlankLength[i];
            vector<double*> &ToFlankFrac = FlankFrac[i];
            vector<double*> &ToFlankLogOddFrac = FlankLogOddFrac[i];


            ToFlank.resize(NoMarkers);
            ToFlankFrac.resize(NoMarkers);
            ToFlankLogOddFrac.resize(NoMarkers);

            int j=0;
            while(j<NoMarkers)
            {


                assert(tempMapper[j]<FromFlank.size());

                ToFlank[j]=&FromFlank[tempMapper[j]];
                ToFlankFrac[j]=&FromFlankFrac[tempMapper[j]];
//                ToFlankLogOddFrac[j]=&FromFlankFrac[tempMapper[j]];
                j++;
            }

    }




    }


//    if(SampleID==5)
//        abort();


}



double LogOddsModel::metaInitialize(int SampleId, vector<HaplotypeSet> &InputData,
                                MetaMinimac *const ThisStudy,ThisChunk &MyChunk)
{

    NoStudies=InputData.size();


    NoParams=NoStudies;
    NoDimensions=NoStudies;

    NoMarkers=MyChunk.NoGenoAllStudies;
    SampleID=SampleId;

    LooDosageVal.resize(NoStudies);
    for(int i=0;i<NoStudies;i++)
    {
        LooDosageVal[i].resize(NoMarkers);

        for(int j=0; j<NoMarkers; j++)
        {
            assert(MyChunk.StartWithWindowIndex+j<ThisStudy->NoCommonGenoVariants);
            LooDosageVal[i][j]=InputData[i].LooDosage[SampleID][MyChunk.StartWithWindowIndex+j];
        }
    }

    ChipGTVal.resize(NoMarkers);
    for(int j=0; j<NoMarkers; j++)
    {
        assert(MyChunk.StartWithWindowIndex+j<ThisStudy->NoCommonGenoVariants);
        ChipGTVal[j]=InputData[0].TypedGT[SampleID][MyChunk.StartWithWindowIndex+j];
    }

}

double LogOddsModel::operator()(vector<double> x)
{


    vector<double> tempVar(NoDimensions);
    double sum=0.0;


    logitTransform(x,tempVar);

    for(int ThisMarker=0;ThisMarker<NoMarkers;ThisMarker++)
    {

        double temp=0.0;

        for(int j=0;j<NoDimensions;j++)
        {
            temp+=((tempVar[j])*(LooDosageVal[j][ThisMarker]));
        }

        temp-=(ChipGTVal)[ThisMarker];
        temp=temp*temp;
        sum+=temp;

    }


 //   cout<<x[0]<<"\t"<<x[1]<<"\t"<<sum<<endl;
  return sum;
}





double LeastSquareWithRecombination::initialize(int SampleId, vector<HaplotypeSet> &InputData,
                                    MetaMinimac *const ThisStudy,ThisChunk &MyChunk)
{

    Method=ThisStudy->Method;
    NoStudies=InputData.size();
    NoDimensions=NoStudies; // 1 extra for the flank-Length parameter
    NoMarkers=MyChunk.NoGenoAllStudies;
    SampleID=SampleId;


    LooDosage.resize(NoStudies);
    for(int i=0;i<NoStudies;i++)
    {
        LooDosage[i]=&(InputData[i].LooDosage[SampleID]);
    }
    ChipGT=&InputData[0].TypedGT[SampleID];


    FlankLength.resize(NoStudies);
    FlankFrac.resize(NoStudies);
    RecomSum.resize(NoStudies);

    for(int i=0;i<NoStudies;i++)
    {
        vector<int> &tempMapper=MyChunk.ThisChunkInterAllTypedSitesReverseMap[i];
        vector<int> &FromFlank = InputData[i].FlankLength[SampleID];
        vector<double> &FromFlankFrac = InputData[i].FlankFrac[SampleID];
        vector<int*> &ToFlank = FlankLength[i];
        vector<double*> &ToFlankFrac = FlankFrac[i];

        ToFlank.resize(NoMarkers);
        ToFlankFrac.resize(NoMarkers);

        int j=0;
        while(j<NoMarkers)
        {


            assert(tempMapper[j]<FromFlank.size());

            ToFlank[j]=&FromFlank[tempMapper[j]];
            ToFlankFrac[j]=&FromFlankFrac[tempMapper[j]];
            j++;
        }

    }

//    if(SampleID==5)
//        abort();


}








double LeastSquareWithRecombination::FixedSSQ()
{

    double sum=0.0;
//    vector<double> tempVar(NoDimensions+1);
//
//    logTransform(x,tempVar,NoDimensions);



//
//    for(int i=0; i < NoDimensions; ++i)
//    {
//        double logit = 1.0/(1.0+exp(0.0-x[i]));
//        tempVar[i]=p*logit;
//        cout<<tempVar[i]<<"\t";
//        p = p*(1.-logit);
//    }
//    tempVar[NoDimensions]=p;
//    cout<<tempVar[NoDimensions];

    for(int ThisMarker=0;ThisMarker<NoMarkers;ThisMarker++)
    {
        double temp=0.0;


//        double Sum=0.0;
//
//        for(int j=0;j<NoDimensions;j++)
//        {
//           Sum+=(1+*(FlankLength[j][ThisMarker]));
//        }

//if(ThisMarker==2)
//    abort();
//        cout<<Sum<<endl;

//        cout<<ThisMarker<<"\t";
        for(int j=0;j<NoDimensions;j++)
        {

//            double FlankEffect=tempVar[NoDimensions]*((double)(1+*(FlankLength[j][ThisMarker]))/Sum);


            double FlankEffect=(*(FlankFrac[j][ThisMarker]));





//            cout<<(*(FlankFrac[j][ThisMarker]))<<"\t-\t";



//            cout<<((double)(*(FlankLength[j][ThisMarker])+1)/Sum)<<"\t";
            temp+=((FlankEffect)*((*LooDosage[j])[ThisMarker]));


//        cout<<" THERAPY = "<<j<<"\t"<<*(FlankFrac[j][ThisMarker])<<"\t"<<(*LooDosage[j])[ThisMarker]<<endl;



//        cout<<" D = "<<temp<<"\t";

        }
//        cout<<endl;
//        cout<<"s= "<<temp<<"\t";


        temp-=(*ChipGT)[ThisMarker];

        temp=temp*temp;

        sum+=temp;


//        cout<<" FUNCTION = "<<ThisMarker<<"\t"<<temp<<"\t"<<sum<<endl;
    }
//  cout<<"\t= "<<sum<<"\n";
  return sum;
}



double LeastSquareWithRecombination::operator()(vector<double> x)
{

    double sum=0.0;
    vector<double> tempVar(NoDimensions+1);

    logTransform(x,tempVar,NoDimensions);



//
//    for(int i=0; i < NoDimensions; ++i)
//    {
//        double logit = 1.0/(1.0+exp(0.0-x[i]));
//        tempVar[i]=p*logit;
//        cout<<tempVar[i]<<"\t";
//        p = p*(1.-logit);
//    }
//    tempVar[NoDimensions]=p;
//    cout<<tempVar[NoDimensions];

    for(int ThisMarker=0;ThisMarker<NoMarkers;ThisMarker++)
    {
        double temp=0.0;


        double Sum=0.0;

        for(int j=0;j<NoDimensions;j++)
        {
           Sum+=(1+*(FlankLength[j][ThisMarker]));
        }

//if(ThisMarker==2)
//    abort();
//        cout<<Sum<<endl;

//        cout<<ThisMarker<<"\t";
        for(int j=0;j<NoDimensions;j++)
        {

//            double FlankEffect=tempVar[NoDimensions]*((double)(1+*(FlankLength[j][ThisMarker]))/Sum);


            double FlankEffect=tempVar[NoDimensions]*(*(FlankFrac[j][ThisMarker]));





//            cout<<(*(FlankFrac[j][ThisMarker]))<<"\t-\t";



//            cout<<((double)(*(FlankLength[j][ThisMarker])+1)/Sum)<<"\t";
            temp+=((FlankEffect+tempVar[j])*((*LooDosage[j])[ThisMarker]));


//        cout<<" THERAPY = "<<j<<"\t"<<*(FlankFrac[j][ThisMarker])<<"\t"<<(*LooDosage[j])[ThisMarker]<<endl;



//        cout<<" D = "<<temp<<"\t";

        }
//        cout<<endl;
//        cout<<"s= "<<temp<<"\t";


        temp-=(*ChipGT)[ThisMarker];

        temp=temp*temp;

        sum+=temp;


//        cout<<" FUNCTION = "<<ThisMarker<<"\t"<<temp<<"\t"<<sum<<endl;
    }
//  cout<<"\t= "<<sum<<"\n";
  return sum;
}






#endif // MINIMIZERS_H_INCLUDED
