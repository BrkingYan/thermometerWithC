package com.example.yy.thermometerwithc.Activity;

import android.content.Intent;
import android.support.v4.app.Fragment;
import android.support.v4.app.FragmentManager;
import android.support.v4.app.FragmentTransaction;
import android.support.v7.app.AppCompatActivity;
import android.os.Bundle;
import android.util.Log;
import android.view.View;
import android.widget.Button;

import com.example.yy.thermometerwithc.Algorithm;
import com.example.yy.thermometerwithc.R;

public class MainActivity extends AppCompatActivity {


    private Button configParamButton;
    private Button mixFreButton;



    private static final String TAG = "MainActivity";


    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        configParamButton = findViewById(R.id.config_param_button);
        configParamButton.setOnClickListener(clickListener);

        //trainingButton.setOnClickListener(clickListener);
        //mixFreButton = (Button)findViewById(R.id.solution_mixfre);
        //mixFreButton.setOnClickListener(clickListener);


    }

    /**
     * 定义的clickListener
     * 当点击fstart文本框时，出现供选择的开始频率的下拉列表
     * 当点击duration文本框时，出现供选择的chirp时长的下拉列表
     * 当点击用互相关法测温的button时，启动XcrroActivity
     * 当点击用混频法测温的button时，启动MixFrequnceActivity，并将选择的fstart的值与duration的值传递。
     */
    private View.OnClickListener clickListener = new View.OnClickListener() {
        @Override
        public void onClick(View v) {
            switch (v.getId()) {
                /*case R.id.solution_mixfre:
                    Intent startIntent = new Intent(getApplicationContext(),
                            MixFrequnceActivity.class);
                    startIntent.putExtra(MixFrequnceActivity.RETURN_Fstart, tvFstart
                            .getText().toString());
                    startIntent.putExtra(MixFrequnceActivity.RETURN_Duration, tvDuration
                            .getText().toString());
                    startActivity(startIntent);break;*/
                case R.id.config_param_button:
                    Log.d(TAG,"into case");
                    Intent intent = new Intent(MainActivity.this,ParamConfigActivity.class);
                    startActivity(intent);
                    finish();
                    break;
            }
        }
    };
}
