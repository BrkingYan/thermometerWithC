package com.example.yy.thermometerwithc.Activity;

import android.content.Intent;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.View;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.Button;
import android.widget.EditText;
import android.widget.Spinner;

import com.example.yy.thermometerwithc.R;

import java.util.*;

public class ParamConfigActivity extends AppCompatActivity implements View.OnClickListener {


    private EditText micDistanceText;
    private String configedMicDistance;

    private Button micConfigButton;
    private Button beginMeasureBtn;

    private Spinner bandSpinner;
    private Spinner durationSpinner;

    private List<Integer> bandOptions;
    private List<Double> durationOptions;

    private ArrayAdapter<Integer> bandAdapter;
    private ArrayAdapter<Double> durationAdapter;
    private Integer selectedFreStart;
    private Double selectedDuration;
    private static final String bandKey = "band_key";
    private static final String durationKey = "duration_key";
    private static final String micKey = "mic_key";


    private static final String TAG = "ParamConfig";


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_config);

        micDistanceText = findViewById(R.id.micDistance_input);
        bandSpinner = findViewById(R.id.spinner_band);
        durationSpinner = findViewById(R.id.spinner_duration);

        // 加载Spinner数据
        loadSpinnerData();

        //创建Adapter
        bandAdapter = new ArrayAdapter<>(this,android.R.layout.simple_spinner_item,bandOptions);
        bandAdapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);
        durationAdapter = new ArrayAdapter<>(this,android.R.layout.simple_spinner_item,durationOptions);
        durationAdapter.setDropDownViewResource(android.R.layout.simple_spinner_dropdown_item);

        //绑定Adapter
        bandSpinner.setAdapter(bandAdapter);

        durationSpinner.setAdapter(durationAdapter);
        //设置Spinner的点击事件
        setOnItemSelectedListener(bandSpinner,1);
        setOnItemSelectedListener(durationSpinner,2);


        micConfigButton = findViewById(R.id.button_set_distance);
        micConfigButton.setOnClickListener(this);
        beginMeasureBtn = findViewById(R.id.solution_mixfre);
        beginMeasureBtn.setOnClickListener(this);

    }

    private void loadSpinnerData(){
        bandOptions = new ArrayList<>();
        durationOptions = new ArrayList<>();
        for (int i = 0;i<11;i++){
            bandOptions.add(2000*i);
        }
        for (int i = 0;i<=10;i++){
            double ele = 0.01*i;
            durationOptions.add(ele);
        }
    }

    private void setOnItemSelectedListener(Spinner spinner, final int spinnerID){
        spinner.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {
            @Override
            public void onItemSelected(AdapterView<?> parent, View view, int position, long id) {
                String str = parent.getItemAtPosition(position).toString();
                // 捕捉用户设置的band和duration
                if (spinnerID==1) {
                    selectedFreStart = bandOptions.get(position);
                }
                if (spinnerID==2){
                    selectedDuration = durationOptions.get(position);
                }
            }
            @Override
            public void onNothingSelected(AdapterView<?> parent) {

            }
        });
    }


    @Override
    public void onClick(View v) {
        switch (v.getId()){
            case R.id.button_set_distance:
                configedMicDistance = micDistanceText.getText().toString();
                Log.d(TAG,"distance:" + configedMicDistance);
                break;
            case R.id.solution_mixfre:
                Intent intent = new Intent(ParamConfigActivity.this,MeasureActivity.class);
                intent.putExtra(bandKey,selectedFreStart);
                intent.putExtra(durationKey,selectedDuration);
                if (configedMicDistance==null || configedMicDistance.equals("")){
                    intent.putExtra(micKey,0.138);
                }else {
                    intent.putExtra(micKey,Double.parseDouble(configedMicDistance));
                }
                startActivity(intent);
                finish();
                break;
        }
    }
}
