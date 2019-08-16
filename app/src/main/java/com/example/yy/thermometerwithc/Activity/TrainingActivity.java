package com.example.yy.thermometerwithc.Activity;

import android.content.Intent;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.support.annotation.Nullable;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.Button;

import com.example.yy.thermometerwithc.Algorithm;
import com.example.yy.thermometerwithc.AudioPlayer;
import com.example.yy.thermometerwithc.AudioRecorder;
import com.example.yy.thermometerwithc.BorderVar;
import com.example.yy.thermometerwithc.CalculateThread;
import com.example.yy.thermometerwithc.R;
import com.jjoe64.graphview.GraphView;
import com.jjoe64.graphview.series.DataPoint;
import com.jjoe64.graphview.series.LineGraphSeries;

public class TrainingActivity extends AppCompatActivity implements CalculateThread.TrainingGraphDataCallBack,BorderVar {

    private double bandWidth;
    private double duration;
    private double micDistance;

    private static final String bandKey = "band_key";
    private static final String durationKey = "duration_key";
    private static final String micKey = "mic_key";


    private LineGraphSeries<DataPoint> freSeries;
    private LineGraphSeries<DataPoint> timeFreSeries;
    private GraphView spectrumGraph;
    private GraphView timeFreGraph;
    private double[] spectrumData;
    private double[] timeFreData;

    private Button startTraning;
    private Button stopTraining;


    private CalculateThread mCalculateThread;
    private AudioRecorder mRecorder;
    private AudioPlayer mPlayer;

    private static final String TAG = "training";

    private Handler mHandler = new Handler(){
        @Override
        public void handleMessage(Message msg) {
            switch (msg.what){
                case 3:
                    updateSpecView(spectrumData);
                    updateTimeFreView(timeFreData);
                    break;

            }
        }
    };


    static {
        System.loadLibrary("JNITest");
    }

    @Override
    protected void onCreate(@Nullable Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_training);
        //从上个Activity获取Intent数据
        Intent intent = getIntent();
        if (intent!=null){
            bandWidth = endFrequnce - intent.getIntExtra(bandKey,0);
            duration = intent.getDoubleExtra(durationKey,0);
            micDistance = intent.getDoubleExtra(micKey,0);
        }else {
            Log.e(TAG,"intent is null");
        }
        Log.d(TAG,"bandWidth: "+bandWidth);
        Log.d(TAG,"duration: "+duration);
        Log.d(TAG,"mic distance: "+micDistance);

        spectrumGraph = findViewById(R.id.graph_spectrum_training);
        timeFreGraph = findViewById(R.id.graph_timeFre);
        startTraning = findViewById(R.id.startTraning_btn);
        stopTraining = findViewById(R.id.stopTraning_btn);

        init();

        startTraning.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {

                initCalculateThread();
                mCalculateThread.start();

                mRecorder = new AudioRecorder();
                mRecorder.registerCallback(mCalculateThread);

                if (!mRecorder.isRecording()) {
                    mRecorder.startRecord();
                }

                mPlayer = new AudioPlayer(TrainingActivity.this,mRecorder);
                mPlayer.initAudioTrack();
                //mPlayer.playWavFile();
                mPlayer.playChirpData();
            }
        });


        stopTraining.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                /*mRecorder.finishRecord();
                mPlayer.stopPlay();*/

                /*double[] arr = new double[]{3,4,5,1,2};
                int result = Algorithm.maxIndexInRangeJni(arr,0,3);
                Log.d(TAG,result + "");*/

                double[] arr = new double[]{0,1,2,0,0,0,0,1};
                double[] chirp = new double[]{1,2};
                int index = Algorithm.correlationJni(5,arr,chirp);
                Log.d(TAG,"index!!!!!!!!!!!!!!!!" + index);
            }
        });
    }

    private void init(){

        freSeries = new LineGraphSeries<>(generateData());
        spectrumGraph.addSeries(freSeries);

        timeFreSeries = new LineGraphSeries<>(generateData());
        timeFreSeries.setThickness(2);
        timeFreGraph.addSeries(timeFreSeries);
        //tempTaker.start();
    }

    private void initCalculateThread(){
        mCalculateThread = new CalculateThread();
        mCalculateThread.registerTrainingGraphDataCallback(this);
        mCalculateThread.setParameter(N,duration,bandWidth,micDistance,mHandler);
    }

    /*
    *  产生静态数据
    * */
    private DataPoint[] generateData() {
        int count = 200;
        DataPoint[] values = new DataPoint[count];
        for (int i=0; i<count; i++) {
            double x = i;
            double y = Math.random()+0.3;
            DataPoint v = new DataPoint(x, y);

            values[i] = v;
        }
        return values;
    }


    /**
     * 更新频谱图
     */
    private void updateSpecView(double[] s){
        DataPoint[] values = new DataPoint[s.length];
        for(int i = 0 ; i < s.length ; i++){
            double xx = i;
            double yy = s[i];
            DataPoint vv = new DataPoint(xx, yy);
            values[i] = vv;
        }
        Log.e(TAG,"graph update");
        Log.e(TAG, "sample length = " + s.length);
        freSeries.resetData(values);
    }

    /*
    * 更新时频图
    * */
    private void updateTimeFreView(double[] s){
        DataPoint[] values = new DataPoint[s.length];
        for(int i = 0 ; i < s.length ; i++){
            double xx = i;
            double yy = s[i];
            DataPoint vv = new DataPoint(xx, yy);
            values[i] = vv;
        }
        //Log.e(TAG, "sample length = " + s.length);
        timeFreSeries.resetData(values);
    }

    @Override
    protected void onDestroy() {
        super.onDestroy();
        if(mRecorder != null){
            mRecorder.finishRecord();
        }
        if(mCalculateThread != null){
            mCalculateThread.finish();
        }
        android.os.Process.killProcess(android.os.Process.myPid());
    }


    /*
     *  训练数据回调
     * */
    @Override
    public void onTrainingGraphDataReady(double[] spectrumData, double[] timeFreData) {
        this.spectrumData = spectrumData;
        this.timeFreData = timeFreData;
    }

}
