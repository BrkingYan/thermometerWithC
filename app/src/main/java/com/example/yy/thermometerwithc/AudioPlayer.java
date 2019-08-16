package com.example.yy.thermometerwithc;

import android.content.Context;
import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioTrack;
import android.os.Environment;
import android.util.Log;
import android.widget.Toast;

import com.example.yy.thermometerwithc.Activity.MeasureActivity;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class AudioPlayer {
    private AudioTrack mPlayer;
    private AudioTrack dataPlayer;

    private static final int SAMPLE_RATE = 48000;
    private static final int CHANNEL = AudioFormat.CHANNEL_OUT_MONO;
    private static final int AUDIO_FORMAT = AudioFormat.ENCODING_PCM_16BIT;
    private static final int STREAM_TYPE = AudioManager.STREAM_MUSIC;
    private int mBufferSizeInBytes;
    private volatile Status mStatus = Status.STATUS_NO_READY;

    private ExecutorService exec;

    private AudioRecorder mAudioRecorder;


    private String filePath = Environment.getExternalStorageDirectory() + "/main.wav";

    private Context mMainActivity;
    private static final String TAG = "player";


    public AudioPlayer(Context context,AudioRecorder mAudioRecorder){
        exec = Executors.newCachedThreadPool();
        mMainActivity = context;
        this.mAudioRecorder = mAudioRecorder;
    }

    public void playWavFile(){

        if (mStatus == Status.STATUS_NO_READY || mPlayer == null){
            /*Toast.makeText(mMainActivity,"player is not inited yet",Toast.LENGTH_SHORT).show();
            return;*/
            throw new IllegalStateException("player is not inited yet");
        }
        if (mStatus == Status.STATUS_START){
            /*Toast.makeText(mMainActivity,"play is already start",Toast.LENGTH_SHORT).show();
            return;*/
            throw new IllegalStateException("player is busy now");
        }

        mStatus = Status.STATUS_START;

        exec.execute(new Runnable() {
            @Override
            public void run() {
                DataInputStream dis = null;
                byte[] buffer = new byte[mBufferSizeInBytes];
                //byte[] buffer = new byte[1440];
                int bytesLen;

                try {
                    mPlayer.play();
                    //dataPlayer.play();
                    dis = new DataInputStream(new FileInputStream(filePath));

                    while ((bytesLen = dis.read(buffer)) != -1 && mStatus == Status.STATUS_START){
                        Log.d(TAG,"buffer[]" + buffer[10]);
                        mPlayer.write(buffer,0,bytesLen);
                    }

                    Log.d(TAG,"play over");
                } catch (FileNotFoundException e) {
                    Toast.makeText(mMainActivity, "file not found", Toast.LENGTH_SHORT).show();
                } catch (IOException e) {
                    Toast.makeText(mMainActivity, "IO exception", Toast.LENGTH_SHORT).show();
                }finally {
                    try {
                        if (dis != null)
                            dis.close();
                        mAudioRecorder.finishRecord();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        });
    }

    public void playChirpData(){

        if (mStatus == Status.STATUS_NO_READY || mPlayer == null){
            throw new IllegalStateException("player is not inited yet");
        }
        if (mStatus == Status.STATUS_START){
            throw new IllegalStateException("player is busy now");
        }

        mStatus = Status.STATUS_START;

        exec.execute(new Runnable() {
            @Override
            public void run() {
                DataInputStream dis = null;
                //byte[] buffer = new byte[mBufferSizeInBytes];
                byte[] buffer;
                int bytesLen;

                try {
                    dataPlayer.play();

                    while (mStatus == Status.STATUS_START){
                        buffer = doubleArrayToByteArray(ChirpSignal.upChirpWithInterval(48000,2000,0.01));
                        Log.d(TAG,"buffer play:" + buffer[20] + "buffer length:" + buffer.length);
                        dataPlayer.write(buffer,0,buffer.length);
                    }

                    Log.d(TAG,"play over");
                } finally {
                    try {
                        if (dis != null)
                            dis.close();
                        mAudioRecorder.finishRecord();
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        });
    }

    private byte[] doubleArrayToByteArray(double[] doubles){
        int len = doubles.length;//480
        byte[] byteArray = new byte[len * 8];
        byte[] temp;
        for (int i = 0;i<len;i++){
            temp = doubleToByteArray(doubles[i]);// byte[8]
            for (int j = 0;j < 8;j++){
                byteArray[i* 8 + j] = temp[j];
            }
        }
        Log.d(TAG,"byteArray in method:length : " + byteArray.length);
        return byteArray;
    }

    public static byte[] doubleToByteArray(double Value) {
        long accum = Double.doubleToRawLongBits(Value);
        byte[] byteRet = new byte[8];
        byteRet[0] = (byte) (accum & 0xFF);
        byteRet[1] = (byte) ((accum >> 8) & 0xFF);
        byteRet[2] = (byte) ((accum >> 16) & 0xFF);
        byteRet[3] = (byte) ((accum >> 24) & 0xFF);
        byteRet[4] = (byte) ((accum >> 32) & 0xFF);
        byteRet[5] = (byte) ((accum >> 40) & 0xFF);
        byteRet[6] = (byte) ((accum >> 48) & 0xFF);
        byteRet[7] = (byte) ((accum >> 56) & 0xFF);
        return byteRet;
    }


    public void stopPlay(){
        if (mStatus == Status.STATUS_READY || mStatus == Status.STATUS_NO_READY){
            Toast.makeText(mMainActivity,"play not start yet",Toast.LENGTH_SHORT).show();
            //throw new IllegalStateException("play not start yet");
            return;
        }
        mStatus = Status.STATUS_STOP;
        mPlayer.stop();
        if (mPlayer != null) {
            mPlayer.release();
            mPlayer = null;
            mStatus = Status.STATUS_NO_READY;
        }
    }

    public void initAudioTrack(){
        mBufferSizeInBytes = AudioTrack.getMinBufferSize(SAMPLE_RATE, CHANNEL, AUDIO_FORMAT);
        if (mBufferSizeInBytes <= 0) {
            throw new IllegalStateException("AudioTrack is not available " + mBufferSizeInBytes);
        }

        mPlayer = new AudioTrack(AudioManager.STREAM_MUSIC, SAMPLE_RATE, CHANNEL, AUDIO_FORMAT,
                mBufferSizeInBytes, AudioTrack.MODE_STREAM);
        dataPlayer = new AudioTrack(AudioManager.STREAM_MUSIC, SAMPLE_RATE, CHANNEL, AUDIO_FORMAT,
                mBufferSizeInBytes, AudioTrack.MODE_STREAM);

        mStatus = Status.STATUS_READY;
    }

    public boolean isPlayingFinish(){
        return mStatus != Status.STATUS_OVER;
    }

    public enum Status {
        //未开始
        STATUS_NO_READY,
        //就绪
        STATUS_READY,
        //播放
        STATUS_START,
        //停止
        STATUS_STOP,
        //结束
        STATUS_OVER

    }

}
