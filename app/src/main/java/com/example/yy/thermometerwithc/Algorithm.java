package com.example.yy.thermometerwithc;

public class Algorithm {
    /*cd
    *  将chirp与left进行相关，并找出startIndex
    * */
    public native static int correlationJni(int n,double[] chirp,double[] left);

    /*
    *  混频
    * */
    public native static double[] mixFrequenceJni(double[] s1, double[] s2);

    /*
    *  FFT
    * */
    public native static double[] fftJni(double[] signal);

    /*
    *  find max index
    * */
    public native static int maxIndexInRangeJni(double[] signal,int start,int end);

    public native static double[] firHelperJni(double[] in,double[] outData,short[] signal);

}
