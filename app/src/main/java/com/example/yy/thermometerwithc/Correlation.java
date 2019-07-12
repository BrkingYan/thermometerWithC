package com.example.yy.thermometerwithc;

import org.jtransforms.fft.DoubleFFT_1D;

import java.util.Arrays;


public class Correlation {

    private static DoubleFFT_1D fftClass = null;
    private static int N = 0;
    private static double[] x_fft = null;
    private static double[] y_fft = null;

    public Correlation(int N) {
        fftClass = new DoubleFFT_1D(N);

        this.N = N;
        // variables to store fft results
        x_fft = new double[N * 2];
        y_fft = new double[N * 2];
    }

    /*
    function: calculate the correlation using FFT
             correlation = IFFT (FFT(x) * conj(FFT(y)))
    input: double array value x and y
    output: correlation value

        PerformCorrelation(Chirp, left)
    */
    public static double[] PerformCorrelation(double[] x, double[] y) {

        // this is a very important initialization process
        Arrays.fill(x_fft, 0, N, 0);
        Arrays.fill(y_fft, 0, N, 0);

        // filling the array with data
        System.arraycopy(x, 0, x_fft, 0, x.length);
        System.arraycopy(y, 0, y_fft, 0, y.length);
        double[] ret = new double[N * 2];


        // compute the fft results
        fftClass.realForward(x_fft);
        fftClass.realForward(y_fft);

        // compute the conj of signal y_fft
        for (int i = 0; i < N / 2; i++) {
            y_fft[2 * i + N] = 0 - y_fft[2 * i + N];
        }

        // compute the multiplication
        for (int i = 0; i < N / 2; i++) {
            ret[i * 2] = x_fft[i * 2] * y_fft[i * 2] + x_fft[i * 2 + 1] * y_fft[i * 2 + 1];
            ret[i * 2 + 1] = x_fft[i * 2] * y_fft[i * 2 + 1] + x_fft[i * 2 + 1] * y_fft[i * 2];
        }

        // compute the ifft of the results
        fftClass.complexInverse(ret, false);
        double[] result = new double[N];
        System.arraycopy(ret, 0, result, 0, N);
        return result;
    }
}