package com.example.yy.thermometerwithc;

import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder;
import android.os.Process;
import android.util.Log;

import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.concurrent.LinkedBlockingQueue;

public class AudioRecorder implements IAudioRecorder{

    public static final int RECORDER_SAMPLE_RATE = 48000;
    public static final int RECORDER_CHANNELS = AudioFormat.CHANNEL_OUT_STEREO;
    public static final int RECORDER_AUDIO_ENCODING = AudioFormat.ENCODING_PCM_16BIT;


    private static final int BUFFER_BYTES_ELEMENTS = 1024;
    private static final int BUFFER_BYTES_PER_ELEMENT = RECORDER_AUDIO_ENCODING;
    private static final int RECORDER_CHANNELS_IN = AudioFormat.CHANNEL_IN_STEREO;


    public static final int RECORDER_STATE_FAILURE = -1;
    public static final int RECORDER_STATE_IDLE = 0;
    public static final int RECORDER_STATE_STARTING = 1;
    public static final int RECORDER_STATE_STOPPING = 2;
    public static final int RECORDER_STATE_BUSY = 3;

    private volatile int recorderState;

    private final Object recorderStateMonitor = new Object();

    private static final String TAG = "recorder";

    private RecordingCallback recordingCallback;


    @SuppressWarnings("ResultOfMethodCallIgnored")
    private void onRecordFailure() {
        recorderState = RECORDER_STATE_FAILURE;
        finishRecord();
    }

    @Override
    public void startRecord() {
        if (recorderState != RECORDER_STATE_IDLE) {
            return;
        }

        try {
            recorderState = RECORDER_STATE_STARTING;
            startRecordThread();
        } catch (FileNotFoundException e) {
            onRecordFailure();
            e.printStackTrace();
        }
    }

    private void startRecordThread() throws FileNotFoundException {
        Log.e(TAG,"record thread run");
        new Thread(new PriorityRunnable(Process.THREAD_PRIORITY_AUDIO) {

            private void onExit() {
                synchronized (recorderStateMonitor) {
                    recorderState = RECORDER_STATE_IDLE;
                    recorderStateMonitor.notifyAll();
                }
            }


            @SuppressWarnings("ResultOfMethodCallIgnored")
            @Override
            public void runImpl() {
                Log.e(TAG,"into rumImpl");
                int bufferSize = Math.max(BUFFER_BYTES_ELEMENTS * BUFFER_BYTES_PER_ELEMENT,
                        AudioRecord.getMinBufferSize(RECORDER_SAMPLE_RATE, RECORDER_CHANNELS_IN, RECORDER_AUDIO_ENCODING));
                Log.e(TAG,"buffersize:" + bufferSize);
                int size = 1024;
                while(size < bufferSize){
                    size = size * 2;
                }
                bufferSize = size;

                AudioRecord recorder = new AudioRecord(MediaRecorder.AudioSource.MIC, RECORDER_SAMPLE_RATE, RECORDER_CHANNELS_IN, RECORDER_AUDIO_ENCODING, bufferSize);
                boolean flag = recorder == null;
                Log.e(TAG,"recoder is null:" + flag);
                if(recorder.getState() == AudioRecord.STATE_UNINITIALIZED){
                    Log.e(AudioRecorder.class.getSimpleName(),"*******************************Initialize audio recorder error");
                    return;
                }else{
                    Log.d(AudioRecorder.class.getSimpleName(),"-------------------------------Initialize AudioRecord ok");
                }
                try {
                    if (recorderState == RECORDER_STATE_STARTING) {
                        recorderState = RECORDER_STATE_BUSY;
                    }
                    recorder.startRecording();

                    short recordBuffer[] = new short[bufferSize];
                    do {
                        /**************************** 读数据 ******************************/
                        int bytesRead = recorder.read(recordBuffer, 0, bufferSize);

                        if (bytesRead > 0) {
                            /********************* 回调音频数据 *************************/
                            Log.d(TAG,"bufferSize:" + recordBuffer.length);
                            recordingCallback.onDataReady(recordBuffer,bytesRead);
                        } else {
                            Log.e(AudioRecorder.class.getSimpleName(), "error: " + bytesRead);
                            onRecordFailure();
                        }
                    } while (recorderState == RECORDER_STATE_BUSY);
                } finally {
                    recorder.release();
                }
                onExit();
            }
        }).start();
    }

    @Override
    public void finishRecord() {
        int recorderStateLocal = recorderState;
        if (recorderStateLocal != RECORDER_STATE_IDLE) {
            synchronized (recorderStateMonitor) {
                recorderStateLocal = recorderState;
                if (recorderStateLocal == RECORDER_STATE_STARTING
                        || recorderStateLocal == RECORDER_STATE_BUSY) {

                    recorderStateLocal = recorderState = RECORDER_STATE_STOPPING;
                }

                do {
                    try {
                        if (recorderStateLocal != RECORDER_STATE_IDLE) {
                            recorderStateMonitor.wait();
                        }
                    } catch (InterruptedException ignore) {
                        /* Nothing to do */
                    }
                    recorderStateLocal = recorderState;
                } while (recorderStateLocal == RECORDER_STATE_STOPPING);
            }
        }
    }

    @Override
    public boolean isRecording() {
        return recorderState != RECORDER_STATE_IDLE;
    }

    public interface RecordingCallback {
        void onDataReady(short[] data, int bytelen);
    }

    public void registerCallback(RecordingCallback callback){
        this.recordingCallback = callback;
    }

}
