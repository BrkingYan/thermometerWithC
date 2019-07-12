package com.example.yy.thermometerwithc;

import org.jtransforms.fft.DoubleFFT_1D;

import java.util.Arrays;

/**
 * Created by cc on 2017/9/28.
 */

public class Spectrum {
    //private short[] pcm = null;

    private double signal[];
    private int N;
    private int Fs;
    private double df;
    private double f[];
    private double response[];
    DoubleFFT_1D fftClass = null;
    private int i = 0;

    /**
     *  initialization for fft computation
     * @param N : number of points to compute
     * @param Fs : sampling rate
     */

    public Spectrum(int N, int Fs){
        this.N = N;
        this.signal = new double[N];

        this.Fs = Fs;
        this.df = this.Fs/N;
        f = new double[N/2];
        for(i=0; i<N/2; i++){
            f[i] = df*i;
        }
        response = new double[N/2];
        fftClass = new DoubleFFT_1D(N);

    }

    /**
     * calculate FFT and the power spectrum
     * @param s
     */

    /**************************************/
    public void fft(double[] s) {

        //System.arraycopy( s, 0, this.signal,0, s.length );
        for (int i = 0; i < s.length; i++) {
            signal[i] = s[i];
        }
        fftClass.realForward(signal);
        for (i = 0; i < N / 2; i++) {
            response[i] = Math.sqrt(signal[2 * i] * signal[2 * i] +
                    signal[2 * i + 1] * signal[2 * i + 1]);
        }
    }
    public double[] getFreqResponse(){
        return response;
    }
}
