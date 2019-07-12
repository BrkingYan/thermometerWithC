package com.example.yy.thermometerwithc;

/**
 * Created by Su Liu on 2018/2/2.
 */
public class ChirpSignal {
    private static double PI = Math.PI;
    private static int FEND = 22000;

    public static double[] upChirp(double fs,int fstart,double T){
        int len = (int)(T*fs);
        double[] t = new double[len];
        double[] chirp = new double[len];
        for (int n = 0;n<len-1;n++){
            double k = (FEND-fstart)/T;
            t[n]=(double)n/fs;
            chirp[n] = Math.sin(2*PI*fstart*t[n]+ PI*k*t[n]*t[n]);
            System.out.println(t[n]);
        }
        System.out.println(len);

        return chirp;
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
