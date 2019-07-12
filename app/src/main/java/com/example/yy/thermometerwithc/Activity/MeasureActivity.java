package com.example.yy.thermometerwithc.Activity;

import android.content.Intent;
import android.os.Build;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.support.annotation.Nullable;
import android.support.annotation.RequiresApi;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.*;

import com.example.yy.thermometerwithc.Algorithm;
import com.example.yy.thermometerwithc.AudioPlayer;
import com.example.yy.thermometerwithc.AudioRecorder;
import com.example.yy.thermometerwithc.BorderVar;
import com.example.yy.thermometerwithc.CalculateThread;
import com.example.yy.thermometerwithc.R;
import com.jjoe64.graphview.GraphView;
import com.jjoe64.graphview.series.DataPoint;
import com.jjoe64.graphview.series.LineGraphSeries;

public class MeasureActivity extends AppCompatActivity implements CalculateThread.GraphDataCallBack,BorderVar {

    /*static {
        System.loadLibrary("JNITest");
    }*/

    private double bandWidth;
    private double duration;
    private double micDistance;

    private Button startBtn;
    private Button stopBtn;
    private Button calBtn;
    private TextView mTemperatureResult;
    private int mTemperature;
    private TextView mFResult;
    private LineGraphSeries<DataPoint> mSeries;
    private LineGraphSeries<DataPoint> mFpSeries;
    private GraphView spectrumGraph;
    private GraphView frePeakGraph;

    private CalculateThread mCalculateThread;
    private AudioRecorder mRecorder;
    private AudioPlayer mPlayer;

    private double[] graphData;
    private double[] fpGraphData;


    private static final String bandKey = "band_key";
    private static final String durationKey = "duration_key";
    private static final String micKey = "mic_key";
    private static final String TAG = "MeasureActivity";

    int globalCount = 1;
    private Handler mHandler = new Handler(){
        @Override
        public void handleMessage(Message msg) {
            switch (msg.what){
                case 1:
                    double temp = (double)msg.obj;
                    if (Double.isInfinite(temp)){
                        mTemperature = 0;
                    }else {
                        mTemperature = (int)((double)msg.obj * 10);
                    }
                    //mTemperatureResult.setText("" + mTemperature/10);
                    mTemperatureResult.setText("" + (globalCount++));
                    //updateView(mCalculateThread.graphData());
                    updateFpView(mCalculateThread.fpGraphData());
                    break;
                case 2:
                    int f = (int)msg.obj;
                    mFResult.setText("" + f);
            }
        }
    };

    @RequiresApi(api = Build.VERSION_CODES.JELLY_BEAN)
    @Override
    protected void onCreate(@Nullable Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_measurement);
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


        //初始化视图
        spectrumGraph = findViewById(R.id.graph_spectrum);
        frePeakGraph = findViewById(R.id.graph_frePeak);
        mTemperatureResult = findViewById(R.id.temperature_result);
        mFResult = findViewById(R.id.calibrate_result);
        startBtn = findViewById(R.id.btn_start);
        stopBtn = findViewById(R.id.btn_stop);
        //calBtn = findViewById(R.id.btn_cal);


        init();
        startBtn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {

                initCalculateThread();
                mCalculateThread.start();

                mRecorder = new AudioRecorder();
                mRecorder.registerCallback(mCalculateThread);

                if (!mRecorder.isRecording()) {
                    mRecorder.startRecord();
                }

                mPlayer = new AudioPlayer(MeasureActivity.this,mRecorder);
                mPlayer.initAudioTrack();
                mPlayer.playWavFile();
            }
        });

        stopBtn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                /*mRecorder.finishRecord();
                mPlayer.stopPlay();*/

                double[] arr = new double[]{3,4,5,1,2};
                int result = Algorithm.maxIndexInRangeJni(arr,0,3);
                Log.d(TAG,result + "");

            }
        });

        /*calBtn.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                //mCalculateThread.start();
            }
        });*/
    }

    private void init(){

        mSeries = new LineGraphSeries<>(generateData());
        spectrumGraph.addSeries(mSeries);

        mFpSeries = new LineGraphSeries<>(generateData());
        mFpSeries.setThickness(2);
        frePeakGraph.addSeries(mFpSeries);
        //tempTaker.start();
    }

    private void initCalculateThread(){
        mCalculateThread = new CalculateThread();
        mCalculateThread.registerGraphDataCallback(this);
        mCalculateThread.setParameter(N,duration,bandWidth,micDistance,mHandler);
    }


    /**
     * generate data for spectrum
     * @pram double[]
     * @return
     */
    private void updateView(double[] s){
        DataPoint[] values = new DataPoint[s.length];
        for(int i = 0 ; i < s.length ; i++){
            double xx = i;
            double yy = s[i];
            DataPoint vv = new DataPoint(xx, yy);
            values[i] = vv;
        }
        Log.e(TAG,"graph update");
        Log.e(TAG, "sample length = " + s.length);
        mSeries.resetData(values);
    }

    private void updateFpView(double[] s){
        DataPoint[] values = new DataPoint[s.length];
        for(int i = 0 ; i < s.length ; i++){
            double xx = i;
            double yy = s[i];
            DataPoint vv = new DataPoint(xx, yy);
            values[i] = vv;
        }
        //Log.e(TAG, "sample length = " + s.length);
        mFpSeries.resetData(values);
    }


    @Override
    public void onGraphDataReady(double[] graphData,double[] fpGraphData) {
        this.graphData = graphData;
        this.fpGraphData = fpGraphData;
    }

    /**
     * generate data for test
     * @return
     */
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
     * Kill all the process
     */
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

}
