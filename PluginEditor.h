/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class AutomaticSpeechDetectionAudioProcessorEditor  : public AudioProcessorEditor, public AsyncUpdater
{
public:
    AutomaticSpeechDetectionAudioProcessorEditor (AutomaticSpeechDetectionAudioProcessor&);
    ~AutomaticSpeechDetectionAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;
    void handleAsyncUpdate() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    AutomaticSpeechDetectionAudioProcessor& processor;
    Colour backgroundColour{juce::Colours::black};
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (AutomaticSpeechDetectionAudioProcessorEditor)
};
