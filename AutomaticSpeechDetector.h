/*
  ==============================================================================

    AutomaticSpeechDetector.h
    Created: 7 Oct 2019 7:07:35pm
    Author:  Danjeli Schembri

  ==============================================================================
*/

#pragma once
#include "../JuceLibraryCode/JuceHeader.h"
#include "DeclareAverageSquaredL2NormOfSpectralFlux.h"
#include "DeclareSkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity.h"
#include "DeclarePauseCount.h"
#include "DeclareRhythmicMeasure.h"
#include "DeclarespectralWeights.h"
#include "DeclareLongRhythmicMeasure.h"
#include "DeclarefeatureClassificationBoostingAlgorithm.h"

class AutomaticSpeechDetector
{
public:
	double * ZeroCrossingRate(const float * arr, int startSample, int blocksInFrame, int samplesPerBlock);
	double AverageSquaredL2NormOfSpectralFlux(const float * arr, int startSample, int BlockNumbers, int SampleNumbers);
	double SkewOfRegressiveLineOfBestFitThroughEstimatedSpectralPowerDensity(const float * arr, int startSample, int BlockNumbers, int SampleNumbers);
	double PauseCount(const float * arr, int startSample, int BlockNumbers, int SampleNumbers);
	double SkewCoefficientOfZeroCrossingRate(const float * arr, int startSample, int BlockNumbers, int SampleNumbers);
	double MeanToMedianRatioOfZeroCrossingRate(const float * arr, int startSample, int BlockNumbers, int SampleNumbers);
	double RhythmicMeasure(const float * arr, int startSample, int BlockNumbers, int SampleNumbers);
	double * spectralWeights(const float * arr, int startSample, int BlockNumbers, int SampleNumbers);
	double LongRhythmicMeasure(const float * arr, int startSample, int BlockNumbers, int SampleNumbers, double * previousSpectralWeightValues);
	int featureClassificationBoostingAlgorithm(double F1, double F2, double F3, double F4, double F5, double F6, double F7);
	AutomaticSpeechDetector();
    void process16K(const float * audioData, int numSamples);
    bool getLatestSpeech();
private:
    std::vector<bool> speechHistory;
    
    dsp::FFT fft1024{10};
    dsp::FFT fft512{9};
	static const int numberOfHops = 32;
	static const int hopSize=1024;//sliding frame after 1024 samples everytime
	static const int sizeOfFrame = numberOfHops * hopSize;
	std::array<float, sizeOfFrame> buffer;
	int bufferPointer = 0;
	DeclareAverageSquaredL2NormOfSpectralFlux normSpectralFluxF1;
	
     DeclarespectralWeights spectralWeightsF7;
	 DeclarefeatureClassificationBoostingAlgorithm BoostingAlgo;

};
