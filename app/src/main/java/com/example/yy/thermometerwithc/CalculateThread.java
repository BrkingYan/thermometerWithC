package com.example.yy.thermometerwithc;

import android.os.Handler;
import android.os.Message;
import android.util.Log;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

public class CalculateThread extends Thread implements BorderVar,AudioRecorder.RecordingCallback {

    private int N = 48000;
    private short[] real; // audio data sample & real part of the result
    private boolean isRunning = true;
    private Spectrum mSpectrum;
    private short[] leftChannelAudioData;
    private short[] rightChannelAudioData;
    private double Duration;
    private double Band;
    private double MicDistance = 0.138;
    private double C = 331.3;
    private double[] Chirp;
    private double[] leftSegment;
    private double[] rightSegment;
    private double[] fpData;
    private int fp;
    private int temperature;
    private boolean isBufferFilled = false; // indicate wehter the input buffer is filled
    private boolean isDataReady = false;
    private short[] calData; // buffer to be divided into two channel


    private int curIdx = 0; // pointer of last stored byted in the buffer
    private int FillCount = 1; // how many calls should we wait to fill the buffer for FFT

    double[] mixSignalFFT;

    private static final int TIME_TO_UPDATE_UI = 1;

    private Handler mHandler;

    private GraphDataCallBack mGraphDataCallBack;
    private AudioRecorder.RecordingCallback mRecordingCallback;

    private List<Integer> fList;

    private static final String TAG = "cal";


    static {
        System.loadLibrary("JNITest");
    }

    public void setParameter(int N, double Duration, double Band,double micDistance,Handler handler) {
        this.N = N;
        this.Band = Band;
        this.Duration = Duration;
        this.MicDistance = micDistance;
        this.real = new short[this.N];
        mSpectrum = new Spectrum(N, Fs);
        leftChannelAudioData = new short[4096];//一次性收集的左声道的数据为4096个
        rightChannelAudioData = new short[4096];//一次性收集的右声道的数据为4096个
        Chirp = ChirpSignal.upChirp(Fs, (int) (endFrequnce - this.Band), this.Duration);
        leftSegment = new double[Chirp.length];
        rightSegment = new double[Chirp.length];
        fpData = new double[200];
        calData = new short[FillCount*byteLength];
        this.mHandler = handler;

        fList = new ArrayList<>();
    }

    @Override
    public void run() {
        while (isRunning) {
            if(isBufferFilled == false){
                //fillData()之后
                if (isDataReady) {
                    //mGraphDataCallBack.onGraphDataReady(graphData(),fpGraphData());
                    //文件下载完成后更新UI
                    Message msg = Message.obtain(mHandler);
                    msg.what = TIME_TO_UPDATE_UI;
                    msg.obj = estimate_temperatureJni();
                    mHandler.sendMessage(msg);
                    synchronized (this) {
                        isDataReady = false;
                    }
                }
            }

        }
    }

    /*
     *  从AudioRecorder回调录音数据
     * */
    @Override
    public void onDataReady(short[] data, int bytelen) {
        if(!isBufferFilled){
            System.arraycopy(data, 0, calData, curIdx * bytelen, bytelen);
            curIdx++;
            if(curIdx >= FillCount){
                curIdx = 0;
                isBufferFilled = false;
                synchronized (this) {
                    fillData(calData);
                }
            }
        }
    }

    /*
     * 读取buffer的音频数据，并将左右信道的信号分开
     * 输出到leftChannelAudioData，rightChannelAudioData中
     * */
    public void fillData(short[] x) {
        for (int i = 0; i < x.length / 2; i = i + 2) {
            leftChannelAudioData[i] = x[2 * i];
            leftChannelAudioData[i + 1] = x[2 * i + 2];
            rightChannelAudioData[i] = x[2 * i + 1];
            rightChannelAudioData[i + 1] = x[2 * i + 3];
        }
        synchronized (this) {
            isDataReady = true;
        }
    }


    public double estimate_temperature() {
        //对信号进行滤波
        Highpass a = new Highpass();

        double[] leftFilterData = a.FIRFilter2KStopBand(leftChannelAudioData);
        double[] rightFilterData = a.FIRFilter2KStopBand(rightChannelAudioData);

        //对信号进行同步，选取chirp长度的左右声道信号。
        synchronusSignal(leftFilterData, rightFilterData);
        //Log.e(TAG,"run here");

        double[] mixSignal = mixFrequence(leftSegment, rightSegment);//把左右信道的信号混合

        mSpectrum.fft(mixSignal);//对混合的信号做傅里叶变换
        double[] mixSignalFFT = mSpectrum.getFreqResponse();
        //double[] mixSignalFFT = Algorithm.fftJni(mixSignal);

        fp = maxIndex(mixSignalFFT, 0, searchEndFrequency);
        double distance = 17.216 * 20.075 * fp / 2000000;
        Log.e(TAG, "mic distance : " + distance);

//        String s1 = MixFrequnceActivity.mDistanceEditText.getText().toString();
//        if (!s1.equals(""))
//            MicDistance = Double.parseDouble(s1);
        Log.e(TAG, "" + fpThread(MaxTemperature) + "...fpThread......" + fpThread(MinTemperature));
        //求fp,混合信号频谱最大值所处的频率，只取低频部分的频谱，相当于通过一个低通滤波器。
        Log.e(TAG, fp + "  fp。。。。。。。。。。。。。");
        //Log.e("TAG","Mic distance = " + MicDistance);
        double temperature = ((MicDistance * Band) / (fp * Duration) - C) / 0.606;
        return temperature;
    }

    public double estimate_temperatureJni() {
        //对信号进行滤波
        Highpass a = new Highpass();
//        double[] leftFilterData = a.FIRFilter(leftChannelAudioData);
//        double[] rightFilterData = a.FIRFilter(rightChannelAudioData);
        Log.d(TAG,"leftChannelData:" + leftChannelAudioData[200]);
        double[] leftFilterData = a.FIRFilter2KStopBandWithJni(leftChannelAudioData);
        double[] rightFilterData = a.FIRFilter2KStopBandWithJni(rightChannelAudioData);

        Log.d(TAG,"leftData:" + leftFilterData[200]);
        //对信号进行同步，选取chirp长度的左右声道信号。
        synchronizeSignalWithJni(leftFilterData, rightFilterData);
        //Log.e(TAG,"run here");

        double[] mixSignal = Algorithm.mixFrequenceJni(leftSegment, rightSegment);//把左右信道的信号混合

        mixSignalFFT = Algorithm.fftJni(mixSignal);//对混合的信号做傅里叶变换
        //Log.d(TAG,"mixSignalFFT[200]: " + mixSignalFFT[200]);

        fp = Algorithm.maxIndexInRangeJni(mixSignalFFT, 200, 1200);
        Log.d(TAG,"mixSignalFFT[fp] = " + mixSignalFFT[fp]);
        Log.d(TAG,"mixSignalFFT[200] = " + mixSignalFFT[200]);
        fList.add(fp);
        if (fList.size() == 300){
            Message msg = Message.obtain();
            msg.what = 2;
            msg.obj = calAverage(fList);
            mHandler.sendMessage(msg);
        }

        double distance = 17.216 * 20.075 * fp / 2000000;
        Log.d(TAG, "mic distance : " + distance);

//        String s1 = MixFrequnceActivity.mDistanceEditText.getText().toString();
//        if (!s1.equals(""))
//            MicDistance = Double.parseDouble(s1);
        Log.e(TAG, "" + fpThread(MaxTemperature) + "...fpThread......" + fpThread(MinTemperature));
        //求fp,混合信号频谱最大值所处的频率，只取低频部分的频谱，相当于通过一个低通滤波器。
        Log.e(TAG, fp + "  fp。。。。。。。。。。。。。");
        //Log.e("TAG","Mic distance = " + MicDistance);
        double temperature = ((MicDistance * Band) / (fp * Duration) - C) / 0.606;
        return temperature;
    }

    private int calAverage(List<Integer> list){
        int av = 0;
        int sum = 0;
        for (int i = 0;i<list.size();i++){
            sum += list.get(i);
        }
        av = sum / 1000;
        return av;
    }


    public void synchronizeSignalWithJni(double[] s1,double[] s2){
        //对信号进行分段
        int segLen = Chirp.length * 3;
        double[] left = new double[segLen];
        double[] right = new double[segLen];
        System.arraycopy(s1,segLen,left,0,segLen);
        System.arraycopy(s2,segLen,right,0,segLen);
        int n = segLen + Chirp.length;
        int startIndex = Algorithm.correlationJni(n,Chirp,left);

        System.arraycopy(left, startIndex, leftSegment, 0, leftSegment.length);
        System.arraycopy(right, startIndex, rightSegment, 0, rightSegment.length);
    }
    public void synchronusSignal(double[] s1, double[] s2) {
        //对信号进行分段
        int segLen = Chirp.length * 3;

        double[] left = new double[segLen];
        double[] right = new double[segLen];

        System.arraycopy(s1, segLen, left, 0, segLen);
        System.arraycopy(s2, segLen, right, 0, segLen);

        long tstart = System.currentTimeMillis();
        Correlation t = new Correlation((segLen + Chirp.length) - 1);

        //double[]t1= DSP.xcorr(Chirp, left);
        double[] t1 = new Correlation(segLen + Chirp.length).PerformCorrelation(Chirp, left);
        long tstop = System.currentTimeMillis();
        Log.e("互相关所需的时间", "" + (tstop - tstart));

        int startindex = maxIndex(t1, 0, t1.length);
        if (startindex > segLen) {
            startindex = startindex - segLen;
            if (startindex > segLen - Chirp.length)
                startindex = startindex - Chirp.length;
        } else {
            startindex = segLen - startindex;
            if (startindex > segLen - Chirp.length)
                startindex = startindex - Chirp.length;
        }

        System.arraycopy(left, startindex, leftSegment, 0, leftSegment.length);
        System.arraycopy(right, startindex, rightSegment, 0, rightSegment.length);

    }


    /**
     * @param s1
     * @param s2
     * @return
     */
    public double[] mixFrequence(double[] s1, double[] s2) {
        double[] mix = new double[N];
        for (int i = 0; i < s1.length; i++) {
            mix[i] = s2[i] * s1[i];
        }
        return mix;
    }

    /**
     * @param in
     * @param startIndex
     * @param endIndex
     * @return
     */
    public int maxIndex(double[] in, int startIndex, int endIndex) {
        double max = Double.MIN_VALUE;
        int maxLoc = 0;
        for (int i = startIndex; i < endIndex; i++) {
            if (in[i] > max) {
                max = in[i];
                maxLoc = i;
            }
        }
        return maxLoc;
    }

    public int estimate_frequnce() {
        double[] ff = mSpectrum.getFreqResponse();
        return maxIndex(ff, 0, 20000);
    }



    public synchronized double[] graphData() {

        //double[] ff = mSpectrum.getFreqResponse();
        //double[] ff = mixSignalFFT;
        //return Arrays.copyOfRange(ff,fpThread(MaxTemperature)-100,fpThread(MinTemperature)+100);
        double[] data = new double[512];
        System.arraycopy(mixSignalFFT,0,data,0,data.length);
        return data;
        //return Arrays.copyOfRange(mixSignalFFT, 0, 1200);

    }

    public double[] fpGraphData() {
        System.arraycopy(fpData, 1, fpData, 0, fpData.length - 1);//新的fpdata取上一次的fpdata的第2个到最后一个数
        fpData[fpData.length - 1] = fp;//把最新获得的fp加到fpdata的最后一位
        return fpData;
    }


    public void finish() {
        isRunning = false;
        try {
            this.join();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public int fpThread(double temperature) {
        double c = C + 0.606 * temperature;
        double f = (MicDistance * Band) / (c * Duration);
        return (int) f;
    }



    public interface GraphDataCallBack{
        void onGraphDataReady(double[] graphData, double[] fpGraphData);
    }

    public void registerGraphDataCallback(GraphDataCallBack callBack){
        this.mGraphDataCallBack = callBack;
    }
}
