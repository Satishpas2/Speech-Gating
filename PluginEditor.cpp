/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
AutomaticSpeechDetectionAudioProcessorEditor::AutomaticSpeechDetectionAudioProcessorEditor (AutomaticSpeechDetectionAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (400, 300);
    processor.onSpeechChanged = [this]
    {
        this->backgroundColour = processor.getIsSpeech() ? Colours::green : Colours::red;
        triggerAsyncUpdate();
    };
    
}

AutomaticSpeechDetectionAudioProcessorEditor::~AutomaticSpeechDetectionAudioProcessorEditor()
{
    cancelPendingUpdate();
}

//==============================================================================
void AutomaticSpeechDetectionAudioProcessorEditor::paint (Graphics& g)
{
    g.fillAll (backgroundColour);

}

void AutomaticSpeechDetectionAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
}

void AutomaticSpeechDetectionAudioProcessorEditor::handleAsyncUpdate()
{
    repaint();
}
