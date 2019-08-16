package com.example.yy.thermometerwithc;

import android.util.Log;

/**
 * Created by Su Liu on 2018/2/2.
 */
public class ChirpSignal {
    private static double PI = Math.PI;
    private static int FEND = 22000;

    public static double[] upChirp(double fs,int fstart,double T){
        int len = (int)(T*fs);
        double[] t = new double[len];//len = 480
        double[] chirp = new double[len];
        for (int n = 0;n<len-1;n++){
            double k = (FEND-fstart)/T;
            t[n]=(double)n/fs;
            chirp[n] = Math.sin(2*PI*fstart*t[n]+ PI*k*t[n]*t[n]);
            System.out.println(t[n]);
        }
        Log.d("chirp","chirp len:" + len);
        System.out.println(len);

        return chirp;
    }

    public static double[] upChirpWithInterval(double fs,int fstart,double T){
        int chirpLen = (int)(T*fs);//chirp 长度 10ms
        int intervalLen = (int)(T*fs*2);//间隔长度 20ms
        //int intervalLen = 0;//间隔长度 20ms
        int totalLen = chirpLen + intervalLen;//总长度 30ms
        double[] t = new double[chirpLen];//len = 480 + 960 = 1440
        double[] signal = new double[totalLen];

        double k;
        for (int n = 0;n<chirpLen;n++){
            k = (FEND-fstart)/T;
            t[n]=(double)n/fs;
            signal[n] = Math.sin(2*PI*fstart*t[n]+ PI*k*t[n]*t[n]);
        }

        for (int m = chirpLen;m < totalLen-1;m++){
            signal[m] = 0;
        }
        for (int i = 0;i<signal.length;i++){
            Log.d("chirp","signal:" + signal[i]);
        }
        return signal;
    }

    /*public static short[] upChirp(double fs,int fstart,double T){
        int len = (int)(T*fs);
        double[] t = new double[len];
        short[] chirp = new short[len];
        for (int n = 0;n<len-1;n++){
            double k = (FEND-fstart)/T;
            t[n]=(double)n/fs;
            chirp[n] = (short)(Math.sin(2*PI*fstart*t[n]+ PI*k*t[n]*t[n]));
            System.out.println(t[n]);
        }
        System.out.println(len);
        return chirp;
    }*/
}
