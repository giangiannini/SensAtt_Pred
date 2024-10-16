using System.Collections;
using System.Data;
using UnityEngine;
using UnityEngine.Events;
using bmlTUX.Scripts.ExperimentParts;
using bmlTUX.Scripts.Managers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;



/// <summary>
/// Classes that inherit from Block define custom behaviour for your experiment's blocks.
/// This might be useful for instructions that differ between each block, setting up the scene for each block, etc.
///
/// This template shows how to set up a custom Block script using the toolkit's built-in functions.
///
/// You can delete any unused methods and unwanted comments. The only required part is the constructor.
///
/// You cannot edit the main execution part of Blocks since their main execution is to run their trials.
/// </summary>
public class SensAtt_PredBlock : Block {
    
    SensAtt_PredRunner myRunner;
    private PLT plt;
    private CollisionEventCall_indicator CEC_indicatorL;
    private CollisionEventCall_indicator CEC_indicatorR;
    private CollisionEventCall_indicator CEC_indicator;

    public SensAtt_PredBlock(ExperimentRunner runner, DataTable trialTable, DataRow data, int index) : base(runner, trialTable, data, index) {
        myRunner = (SensAtt_PredRunner)runner;
    }
    private Vector3 InitPositionRight = new Vector3(-0.221f, 0.7496f, 0.005f);
    private Vector3 InitPositionLeft = new Vector3(0.1786f, 0.7496f, 0.005f);
    private Vector3 TargetPositionRight = new Vector3(-0.221f, 0.7496f, 0.005f);
    private Vector3 TargetPositionLeft = new Vector3(0.1786f, 0.7496f, 0.005f);
    private string response_path; 



    protected override void PreMethod() {
        response_path = Application.dataPath + "/Response_training.txt";
        // Set stuff up for unique response in the training phase
        if ((int)Data["Block_num"] == 1)
        {
            //File.AppendAllText(response_training_path, (string)Data["ID"] + "\t" + (string)Data["initials"] + "\n");
            File.AppendAllText(response_path, "Event_Name" + "\t" + "Time" + "\t" + "\n");
            File.AppendAllText(response_path, "Training Start" + "\t" + Time.time * 1000 + "\n");
        }

        plt = myRunner.EmptyObject.GetComponent<PLT>();
        CEC_indicatorR = myRunner.RingR.GetComponent<CollisionEventCall_indicator>();
        CEC_indicatorL = myRunner.RingL.GetComponent<CollisionEventCall_indicator>();
        myRunner.Fix.SetActive(false);
        myRunner.RingR.SetActive(true);
        myRunner.RingL.SetActive(true);

        if ((int)Data["Trial_type"] == 1)
        {
            myRunner.Visual.transform.position = InitPositionRight;
            myRunner.Visual.SetActive(false);
        }
        else if ((int)Data["Trial_type"] == 2)
        {
            myRunner.Visual.transform.position = InitPositionRight;
            myRunner.Visual.SetActive(false);
        }
        else if ((int)Data["Trial_type"] == 3)
        {
            myRunner.Visual.transform.position = InitPositionLeft;
            myRunner.Visual.SetActive(false);
        }
        else if ((int)Data["Trial_type"] == 4)
        {
            myRunner.Visual.transform.position = InitPositionLeft;
            myRunner.Visual.SetActive(false);
        }
    }

    protected override IEnumerator PreCoroutine() {
        if ((int)Data["Trial_type"] == 1 || (int)Data["Trial_type"] == 2)
        {
            myRunner.Instruction_blockL.SetActive(true);
            CEC_indicator = CEC_indicatorL;
        }
        else if ((int)Data["Trial_type"] == 3 || (int)Data["Trial_type"] == 4)
        {
            myRunner.Instruction_blockR.SetActive(true);
            CEC_indicator = CEC_indicatorR;
        }
        bool waitingforfinger = true;
        //used to makes sure that the finger is on the indicator :)) 
        while (waitingforfinger)
        {
            if (CEC_indicator.Indicator_touch == true)
            {
                yield return new WaitForSeconds(2);

                // if after 1s is still on the indicator, then continue with block
                if (CEC_indicator.Indicator_touch == true)
                {
                    myRunner.Instruction_blockL.SetActive(false);
                    myRunner.Instruction_blockR.SetActive(false);
                    myRunner.Instruction_IndicatorReminder.SetActive(false); 
                    waitingforfinger = false; 
                }
                else if (CEC_indicator.Indicator_touch == false)
                {
                    myRunner.Instruction_IndicatorReminder.SetActive(true);
                }
            }
            yield return null;
        }
    }

    // Optional Post-Block code spanning multiple frames. Useful for Block debrief instructions.
    protected override IEnumerator PostCoroutine() {
        myRunner.Resp25.SetActive(true);
        myRunner.Resp50.SetActive(true);
        myRunner.Resp75.SetActive(true);
        myRunner.Response.SetActive(true);
        //myRunner.Hand.GetComponent<ResponseScript>().enabled(true);
        bool waitForResponse = true; 
        while (waitForResponse)
        {
            if (myRunner.FingerTip.GetComponent<ResponseScript>().response != 0)
            {
                myRunner.Resp25.SetActive(false);
                myRunner.Resp50.SetActive(false);
                myRunner.Resp75.SetActive(false);
                myRunner.Response.SetActive(false);

                File.AppendAllText(response_path, "Response" + "\t" + myRunner.FingerTip.GetComponent<ResponseScript>().response + "\n");

                waitForResponse = false; 
            }
            yield return null; 
        }

        if ((int)Data["Block_num"] == 6 || (int)Data["Block_num"] == 12 || (int)Data["Block_num"] == 18 || (int)Data["Block_num"] == 24)
        {
            plt.PLTsend(127); 
            myRunner.Break.SetActive(true); 
            bool waitForPause = true;
            while (waitForPause)
            {
                if (Input.GetKeyDown(KeyCode.Space))
                {
                    myRunner.Break.SetActive(false);
                    plt.PLTsend(126);
                    waitForPause = false;
                }
                yield return null; 
            }
        }
        
        



        //if ((int)Data["Block_num"] == 16 ) {
        //    myRunner.EndExp.SetActive(true);
        //    yield return new WaitForSeconds(60);
        //    myRunner.EndExp.SetActive(false);
        //}
        //yield return null; // yield return required for coroutine. Waits until next frame
    }

    // Optional Post-Block code.
    protected override void PostMethod() {
        // cleanup code (happens all in one frame at end of block)
    }

}

