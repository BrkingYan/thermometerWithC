package com.example.yy.thermometerwithc;

/**
 * Created by cc on 2017/11/29.
 */

public interface BorderVar {
    // message for the handler
    int MESSAGE_VELOCITY = 1;
    int MESSAGE_POI_FFT = 2;

    // sound speed
    double c = 331.4;

    //sampling rate
    int Fs = 48000;

   double endFrequnce = 22000;
    int searchEndFrequency = 1500;

    // FFT about
    int N = 48000; // N-point FFT
    int byteLength = 8192; // how many bytes in every call
    int FillCount = 1; // how many calls should we wait to fill the buffer for FFT

    double MinTemperature = -20;
    double MaxTemperature = 40;
}
