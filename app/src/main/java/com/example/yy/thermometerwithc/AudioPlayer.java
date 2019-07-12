package com.example.yy.thermometerwithc;

import android.content.Context;
import android.media.AudioFormat;
import android.media.AudioManager;
import android.media.AudioTrack;
import android.os.Environment;
import android.util.Log;
import android.widget.Toast;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class AudioPlayer {
    private AudioTrack mPlayer;

    private static final int SAMPLE_RATE = 48000;
    private static final int CHANNEL = AudioFormat.CHANNEL_OUT_STEREO;
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
                byte[] buffer = new byte[mBufferSizeInBytes];
                int bytesLen;

                try {
                    mPlayer.play();

                    dis = new DataInputStream(new FileInputStream(filePath));

                    while ((bytesLen = dis.read(buffer)) != -1 && mStatus == Status.STATUS_START){
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

    public void stopPlay(){
        if (mStatus == Status.STATUS_READY || mStatus == Status.STATUS_NO_READY){
            Toast.makeText(mMainActivity,"play not start yet",Toast.LENGTH_SHORT).show();
            throw new IllegalStateException("play not start yet");
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
