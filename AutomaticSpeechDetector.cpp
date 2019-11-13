/*
  ==============================================================================

    AutomaticSpeechDetector.cpp
    Created: 7 Oct 2019 7:07:35pm
    Author:  Danjeli Schembri

  ==============================================================================
*/
#include <iostream>
#include <cmath>
#include<cstdlib>
#include<complex>
//Class files containing constructor and destructor for feature funtions
#include "AutomaticSpeechDetector.h"


using namespace std;
typedef complex<double> dcomp;

//-------------------------------------------------------------------------------------------//
//signum function
//-------------------------------------------------------------------------------------------//
int sig(double x)
{
    if (x>0)
    {
        return 1;
    }
    else if(x<0)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}
//-------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------//
// Median function
//-------------------------------------------------------------------------------------------//
double Median(double *arr, int n) // n is size of array/vect
{
   // check for even case
    if (n % 2 != 0)
    {
        return (double)arr[n/2];
    }
    else
    {
        return (double)(arr[(n-1)/2] + arr[n/2])/2.0;
    }
}
//-------------------------------------------------------------------------------------------//
// magnitude of the resulting transform coefficients (X) function (Equation 1 has been defined into a function here because it is being used in multiple features . Here power raised is either 1(X^1) or 2 (X^2)
//Also known as spectral power density
//-------------------------------------------------------------------------------------------//
double **magnitudeOfResultingTransformCoefficients(const float *arr,int startSampleNumber, int blocksInFrame, int SamplesPerBlock, int powerRaised)
{
    dcomp  I;
    dcomp a;
    I=-1;
    I=sqrt(I);
    int M=blocksInFrame;
    int N=SamplesPerBlock;
    int row=M;       // 0<=m<M
    int column=N/2;  //0<= k < N/2
    double **X;
    X= new double *[column];
    for (int m=0;m<row;m++)
    {
        X[m]=new double [column];
        for (int k=0;k<column; k++ )
        {
            double X_val=0;
            for (int n=0; n<N;n++)
            {
                int x=startSampleNumber + (m*N)+n;
                double win=(0.54 - 0.46 * cos(2 * juce::MathConstants<float>::pi * n / N ));
                double var= (2*juce::MathConstants<float>::pi*k*n/N);
                a= exp (var*(-1)*I);
                X_val+=(abs(arr[x]*win*a));
            }
            
            X[m][k]=pow(X_val,powerRaised);
        }
    }
    return X;
    
}
//-------------------------------------------------------------------------------------------//
//sorting function
//-------------------------------------------------------------------------------------------//
double *sortAscending(double *arr, int size){
    int i,j;
    for(i=0;i<size;i++)
        for(j=i+1;j<size;j++)
            if(*(arr+j)<*(arr+i))
              {
                double temp=*(arr+j);
                *(arr+j)=*(arr+i);
                *(arr+i)=temp;
                }
   return arr;
}
//-------------------------------------------------------------------------------------------//
//Zero Crossing Rate values.. function
//-------------------------------------------------------------------------------------------//
double *AutomaticSpeechDetector::ZeroCrossingRate(const float *arr, int startSample,int blocksInFrame, int samplesPerBlock)
{
    int M=blocksInFrame;
    int N=samplesPerBlock;
    double *zeroCrossingRate;
    zeroCrossingRate=new double[M];
    for (int m=0;m<M;m++)
    {
             int zeroCrossingRateCoefficient=0;
             for(int n=0;n<N;n++)
             {
                 int x =startSample+ ((m*N)+n);
                 if (x>startSample)
                 {
                     zeroCrossingRateCoefficient+=(abs(sig(arr[x])-sig(arr[x-1])));
                 }
             }
             zeroCrossingRate[m]=(1/(2*N))*zeroCrossingRateCoefficient;
    }
    return zeroCrossingRate;
}

//-------------------------------------------------------------------------------------------//
//process on audio data
// Features DSP Code definitions begin->


//-------------------------------------------------------------------------------------------//
// Feature 1: Average Squared L2 Norm Of Spectral Flux
//-------------------------------------------------------------------------------------------//

//here hopSamples will be 1024, as we will be sliding after every 1024 samples
// therefore (currentSample=hopSamples*numberOfBlocks+samplesinBlock)
double AutomaticSpeechDetector::AverageSquaredL2NormOfSpectralFlux(const  float *arr, int startSample,  int BlockNumbers, int SampleNumbers)
{
    // Constructor Declaration from file DeclareAverageSquaredL2NormOfSpectralFlux.cpp containing class DeclareAverageSquaredL2NormOfSpectralFlux
    
    
     // X is The magnitude of the DFT coefficients  (Equation (1))
    double **X=magnitudeOfResultingTransformCoefficients(arr, startSample, BlockNumbers, SampleNumbers,1);
    int M=BlockNumbers;
    int N=SampleNumbers;
    int row=M;       // 0<=m<M
    int column=N/2;  //0<= k < N/2
    
    // 'Weight' is  the weight for the current block m. (EQUATION 2)
    //'l2Norm'is l2-norm of the weighted spectral flux for block (EQUATION 3)
    
    double sumOfSquaredL2Norm=0;
    for (int m=0;m<row;m++)
    {
        double Weight_val=0;
        double l2_val=0;
        for (int k=0;k<column; k++ )
        {
            if (m==0)
            {
                Weight_val += pow(abs(X[m][k]),(2))/N;
                l2_val+= (pow(abs (0-X[m][k]),2)); //EQUATION 3 ONLY NUMERATOR VALUE
            }
        else
            {
                Weight_val += (pow((abs(X[m-1][k])),(2))+ pow(abs(X[m][k]),(2)))/N;
                l2_val+=(pow(abs ((X[m-1][k])-(X[m][k])),2)); //EQUATION 3 ONLY NUMERATOR VALUE
            }
        }
        normSpectralFluxF1.Weight[m]=Weight_val; // EQUATION (2)
        normSpectralFluxF1.l2Norm[m]=sqrt(l2_val/Weight_val); // Equation (3)
        //the feature for average squared l 2-norm of the weighted spectral flux for frame t.
        sumOfSquaredL2Norm+=normSpectralFluxF1.l2Norm[m];
    }
    return sumOfSquaredL2Norm;
}
//-------------------------------------------------------------------------------------------//



//-------------------------------------------------------------------------------------------//
// Feature 2: Skew of regressive line of best fit through estimated spectral power density
//-------------------------------------------------------------------------------------------//
//feature 2
double AutomaticSpeechDetector::SkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity(const  float *arr,int startSample,  int BlockNumbers, int SampleNumbers)
{
    //Constructor Declaration from file  DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity.cpp containing class  DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity
    DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity F2(BlockNumbers);
    
    //Number of blocks
    int M=BlockNumbers;
    //Number of samples in a block
    int N=SampleNumbers;
    int row= BlockNumbers;
    int column=N/2; // this is the value of k i.e. 0<k<N/2
    //Spectral Power density
    double **X_Square=magnitudeOfResultingTransformCoefficients(arr, startSample, BlockNumbers, SampleNumbers,2);
    
    //sum of k where 0<k<N/2
    double sumK=0;
    double sumK_Square=0;
    for (int k=0;k<column;k++)
    {
        sumK+=k;                  // will be used in Regressive coefficient calculation equation 7
        sumK_Square+=pow(k,2);
    }
    
    double regressiveCoefficientSum=0;
    for (int m=0;m<row;m++)
    {
        //insert class
        F2.Xdb[m]=new double [column];
        double  coefficient1=0;
        double coefficient2=0;
        for (int k=0;k<column; k++ )
        {
            F2.Xdb[m][k]=10*log10(X_Square[m][k]);
             coefficient1+=(N/2)*(k*F2.Xdb[m][k]);
            coefficient2+=sumK*F2.Xdb[m][k];
        }
        double denominator=(N/2)*sumK_Square-pow(sumK,2);
        F2.regressiveCoefficient[m]=(coefficient1-coefficient2)/denominator;
        regressiveCoefficientSum+=F2.regressiveCoefficient[m];
    }
    //equation 8
    double regressiveCoefficientValue=0;
    for (int m=0;m<M;m++)
    {
        double ValueCoefficient=F2.regressiveCoefficient[m]-(regressiveCoefficientSum/M);
        regressiveCoefficientValue+=abs(pow(ValueCoefficient,3));
    }
    //returns value of feature for the given frame
    return log(regressiveCoefficientValue);
        
}
//-------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------//
// Feature 3: Pause Count
//-------------------------------------------------------------------------------------------//
//feature 3
double AutomaticSpeechDetector::PauseCount(const  float *arr,int startSample,  int BlockNumbers, int SampleNumbers)
{
    DeclarePauseCount F3(BlockNumbers);
    int M=BlockNumbers;
    int N=SampleNumbers;
    
    double powerOfEachFrame=0;
    for (int m=0;m<M;m++)
    { double InputSum=0;
        for (int n=0;n<N;n++)
        {
            int x=startSample + (m*N)+n;
            InputSum+=pow(arr[x],2)/N;
        }
        F3.powerOfEachBlock[m]=InputSum;
        powerOfEachFrame+=F3.powerOfEachBlock[m]/M;
    }
    int count=0;
    // feature value is the number of times the value of power for each block is less than or equal to 1/4th value of power of that frame
    for (int m=0;m<M;m++)
    {
        if (F3.powerOfEachBlock[m]<= (1/4)*powerOfEachFrame)
            count++;
    }
    return count;
    
}
//-------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------//
// Feature 4: Skew coefficient of zero crossing rate
//-------------------------------------------------------------------------------------------//
//feature 4
double AutomaticSpeechDetector::SkewCoefficientOfZeroCrossingRate(const  float *arr,int startSample,  int BlockNumbers, int SampleNumbers)
{
    int M=BlockNumbers;
    int N=SampleNumbers;
    //calculation for zero crossing rate
    double *zeroCrossingRate=ZeroCrossingRate(arr, startSample, M, N);
    double zeroCrossingRateSum=0;
    for (int m=0;m<M;m++)
    {
         zeroCrossingRateSum+=zeroCrossingRate[m];
    }
    // the feature for skew coefficient of the zero crossing rate for frame i.e. implementation of equation 10
    double numerator=0;
    double denominator=0;
    for (int m=0;m<M;m++)
    {
        // equation 10 numerator calculation
        double numeratorCoefficient=zeroCrossingRate[m]-zeroCrossingRateSum/M;
        numerator+=pow(numeratorCoefficient,3);
        //equation 10 denominator calculation
        denominator+=pow(numeratorCoefficient,2);
    }
    denominator=pow(denominator,(3/2));
    return (numerator/denominator);
}
//-------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------//
// Feature 5: Mean-to-median ratio of zero crossing rate
//-------------------------------------------------------------------------------------------//
//feature 5
double AutomaticSpeechDetector::MeanToMedianRatioOfZeroCrossingRate(const  float *arr,int startSample, int BlockNumbers, int SampleNumbers)
{
    int M=BlockNumbers;
    //Zero Crossing Rate Values
    double *zeroCrossingRate=ZeroCrossingRate(arr, startSample, BlockNumbers, SampleNumbers);
    //Sorting of zero crossing rate values in ascending order
    double *zeroCrossingRateSorted=sortAscending(zeroCrossingRate, BlockNumbers);
    //Zero Crossing Rate Median Value
    double ZeroCrossingRateMedian= Median(zeroCrossingRateSorted, BlockNumbers);
    //Summation of zero crossing rate value
    double zeroCrossingRateSum=0;
    for (int m=0;m<M;m++)
    {
         zeroCrossingRateSum+=zeroCrossingRate[m];
    }
    //Feature value .....implementation of equation 11
    return ((ZeroCrossingRateMedian)/(zeroCrossingRateSum/M));
}
//-------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------//
// Feature 6: Rhythmic measure
//-------------------------------------------------------------------------------------------//
//feature 6
double AutomaticSpeechDetector::RhythmicMeasure(const  float *arr,int startSample, int BlockNumbers, int SampleNumbers)
{
    DeclareRhythmicMeasure F6(BlockNumbers);
    // Total number of blocks in a frame
    int M=BlockNumbers;
    // Total number of samples in a block
    int N=SampleNumbers;
  
    // mean of samples in a block M
    double meanSamplesinBlockM=0;
    
    //Mean of variance
    double VarianceMean=0;
    //Variance Calculation equation no.12
    for (int m=0;m<M;m++)
    {
        for (int n=0;n<N;n++)
        {
            int x=startSample + (m*N)+n;
            //mean of samples in a block
            meanSamplesinBlockM+=arr[x]/N;
        }
        double tempVarianceValue=0;
        for (int n=0;n<M;n++)
        {
            int x=startSample + (m*N)+n;
            //equation 12 Variance calculation
            double tempVarianceCoefficient=arr[x]-meanSamplesinBlockM;
            tempVarianceValue+=pow(tempVarianceCoefficient,2);
        }
        F6.variance[m]=tempVarianceValue/N;
        VarianceMean+=F6.variance[m];
    }
    // Zero Mean Sequence calculation equation no.13
    for (int m=0;m<M;m++)
    {
        double tempZeroMeanSequence=0;
        tempZeroMeanSequence=F6.variance[m]-VarianceMean;
        F6.zeroMeanSequence[m]=tempZeroMeanSequence;
    }
    //AutoCorrelation Value of zero mean sequence
    for (int l=0;l<M;l++)
    {
        double tempAutoCorrelationValue=0;
        for (int m=0;m<M-l;m++)
        {
            tempAutoCorrelationValue+=F6.zeroMeanSequence[m] * F6.zeroMeanSequence[m+l];
        }
        F6.autoCorrelationOfZeroMeanSequence[l]=(1/M)*tempAutoCorrelationValue;
    }
    int L=10; //Lag of 160 ms (minimum value) i.e. period of most rapid rhythm expected
    //calculation of the feature for short rhythmic measure for frame equation 15
    double LargestElement=F6.autoCorrelationOfZeroMeanSequence[0];
    for (int n=L;n<M;n++)
    {
       if (F6.autoCorrelationOfZeroMeanSequence[n]>LargestElement)
       {
           LargestElement=F6.autoCorrelationOfZeroMeanSequence[n];
       }
    }
    return (LargestElement/F6.autoCorrelationOfZeroMeanSequence[0]);
}
//-------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------//
// Feature 7: Long Rhythmic measure
//-------------------------------------------------------------------------------------------//
//In this feature spectral wieght has been calculated in a seperate function because in calculation of AutoCorrelation of score for Long Rhythmic measure feature, it requires values of spectral weight from the previous frame. and hence it would make it accessible to get the values of previous frame.
//feature 7
double *AutomaticSpeechDetector::spectralWeights(const  float *arr,int startSample,int BlockNumbers, int SampleNumbers)
{
     int M=BlockNumbers;
     int N=SampleNumbers;
     int row=BlockNumbers;
     int column=SampleNumbers/2;
     double **spectralPowerDensity_Squared=magnitudeOfResultingTransformCoefficients(arr, startSample, BlockNumbers, SampleNumbers, 2);
     // an empirically derived constant equal to 0.1.
     double alpha=0.1;
     for (int m=0;m<row;m++)
     {
         spectralWeightsF7.Xdb[m]=new double [column];
        // temporary maximum value for log domain power spectrum variable
         double tempMaximumValue=0;
         for (int k=0;k<column; k++ )
         {
             spectralWeightsF7.Xdb[m][k]=10*log10(spectralPowerDensity_Squared[m][k]);
              //Equation 16 calculation....Maximum log domain power spectrum value
             if (spectralWeightsF7.Xdb[m][k]>tempMaximumValue)
             {
                 tempMaximumValue=spectralWeightsF7.Xdb[m][k];
             }
         }
         spectralWeightsF7.maximumLogDomainPowerSpectrum[m]=tempMaximumValue;
         //Spectral Weights for block m ...equation 17
         double tempSpectralWeightValue=0;
         for (int k=0;k<N/2;k++)
         {
             tempSpectralWeightValue+=(sig(spectralWeightsF7.Xdb[m][k]- spectralWeightsF7.maximumLogDomainPowerSpectrum[m]*alpha) +1 )/2;
         }
         spectralWeightsF7.spectralWeight[m]=tempSpectralWeightValue;
     }
    return spectralWeightsF7.spectralWeight;
}
//-------------------------------------------------------------------------------------------//

// Implementation of Long Rhythmic Measure feature
//-------------------------------------------------------------------------------------------//

double AutomaticSpeechDetector::LongRhythmicMeasure(const  float *arr,int startSample,int BlockNumbers, int SampleNumbers, double *previousSpectralWeightValues)
{
    DeclareLongRhythmicMeasure F7(BlockNumbers);
    int M=BlockNumbers;
    int N=SampleNumbers;
    //Spectral weights for blocks
    double *spectralWeightForBlocks=spectralWeights(arr, startSample, BlockNumbers, SampleNumbers);
    
    for(int m=0;m<M;m++)
    {
        F7.totalSpectralWeights[m]=previousSpectralWeightValues[m];
        F7.totalSpectralWeights[M+m]=spectralWeightForBlocks[m];
    }
    
    //AutoCorrelation Score
     for (int l=0;l<2*M;l++)
        {
            double tempAutoCorrelationValue=0;
            for (int m=0;m<M-l;m++)
            {
                tempAutoCorrelationValue+=F7.totalSpectralWeights[m]* F7.totalSpectralWeights[m+l];
            }
            F7.autoCorrelationScoreValue[l]=(1/2*M)*tempAutoCorrelationValue;
        }
        int LL=10; //Lag of 160 ms (minimum value) i.e. period of most rapid rhythm expected
        //calculation of the feature for short rhythmic measure for frame equation 15
        double LargestElement=F7.autoCorrelationScoreValue[0];
        for (int n=LL;n<M;n++)
        {
           if (F7.autoCorrelationScoreValue[n]>LargestElement)
           {
               LargestElement=F7.autoCorrelationScoreValue[n];
           }
        }
        return (LargestElement/F7.autoCorrelationScoreValue[0]);
    }
    
//-------------------------------------------------------------------------------------------//


//-------------------------------------------------------------------------------------------//
// Boosting Algorithm
//-------------------------------------------------------------------------------------------//
    
//Boosting
int AutomaticSpeechDetector::featureClassificationBoostingAlgorithm(double F1, double F2, double F3, double F4, double F5, double F6, double F7)
{
    
    double featureVector[7]={F1,F2,F3,F4,F5,F6,F7};
    // Final Decision
    int C_finalDecision;
    for (int j=0;j<20;j++)
    {
        int index_feature=BoostingAlgo.featureVectorIndexi_j[j];
        double tempWeightedWeakClassifier=0;
        tempWeightedWeakClassifier=BoostingAlgo.weightingCoefficient[j]* sig(featureVector[index_feature]-BoostingAlgo.Threshold[j]);
        BoostingAlgo.weightedWeakClassifier[j]=tempWeightedWeakClassifier;
    }
    double tempC_final=0;
    for( int j=0;j<20;j++)
    {
        tempC_final+=BoostingAlgo.weightedWeakClassifier[j];
    }
    C_finalDecision=sig(tempC_final);
    return C_finalDecision;
}

AutomaticSpeechDetector::AutomaticSpeechDetector() : normSpectralFluxF1(32), spectralWeightsF7(256)
{
	for(int i =0; i < buffer.size(); i++)
	{
		buffer[i] = 0.f;
	}
}

void AutomaticSpeechDetector::process16K(const float * audioData, int numSamples)
{// therefore (currentSample=hopSamples*numberOfBlocks+samplesinBlock)
   
   
	if(bufferPointer + numSamples >= buffer.size())
	{
		int currentlyAvailable = buffer.size() - bufferPointer;
		int removeFromStart = numSamples - currentlyAvailable;
		memcpy(&buffer[removeFromStart], &buffer[0], bufferPointer - removeFromStart);
		bufferPointer = bufferPointer - removeFromStart;
	}
	for(int i = 0; i < numSamples; i++)
	{
		buffer[bufferPointer + i] = audioData[i];
	}
	bufferPointer += numSamples;
	
    // initialising previous spectral weights value for Long rhythmic measure feature, equal to zero only for first frame which has 128 blocks!
	int M = 128;
    double *previousSpectralWeights;
    previousSpectralWeights=new double[M];
    for (int m=0; m<M;m++)
    {
        previousSpectralWeights[m]=0;
    }
    
    for (int frameNumber=0;frameNumber< numberOfHops; frameNumber++)
    {
        int startSampleNumber= frameNumber*hopSize;
        //Every feature has output as single value for every frame
        // function has following properties:  feature1(array,startSampleNumber, numberOfBlocksPerFrame, NumberOfSamplesinBlock)
        //        feature (AUDIOdata,  start sample number, blocks, samples in a block)
        double F1=AverageSquaredL2NormOfSpectralFlux(&buffer[0], startSampleNumber,  32, 1024);
        double F2=SkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity(&buffer[0], startSampleNumber,64, 512);
        double F3=PauseCount(&buffer[0], startSampleNumber, 128,256);
        double F4=SkewCoefficientOfZeroCrossingRate(&buffer[0], startSampleNumber, 128, 256);
        double F5=MeanToMedianRatioOfZeroCrossingRate(&buffer[0], startSampleNumber, 128, 256);
        double F6=RhythmicMeasure(&buffer[0], startSampleNumber, 128,256);
        double F7= LongRhythmicMeasure(&buffer[0], startSampleNumber,128,256, previousSpectralWeights);
        
        previousSpectralWeights=spectralWeights(&buffer[0],startSampleNumber,128, 256);
        //output for boosting algorithm is either 1 or -1 for every frame
        int featureClassification=featureClassificationBoostingAlgorithm(F1,F2,F3,F4,F5,F6,F7);
		speechHistory.push_back(featureClassification == 1);
    }
   
}

bool AutomaticSpeechDetector::getLatestSpeech()
{
	if(speechHistory.size() > 0)
	{
		return speechHistory.back();
	}
	return false;
}
