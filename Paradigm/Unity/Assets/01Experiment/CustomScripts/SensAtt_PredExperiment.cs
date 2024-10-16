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
/// Classes that inherit from Experiment define custom behaviour for the start and end of your experiment.
/// This might useful for experiment setup, instructions, and debrief.
///
/// This template shows how to set up a custom experiment script using the toolkit's built-in functions.
///
/// You can delete any unused methods and unwanted comments. The only required part is the constructor.
///
/// You cannot edit the main execution part of experiments since their main execution is to run the trials and blocks.
/// </summary>
public class SensAtt_PredExperiment : Experiment {
    //
    SensAtt_PredRunner myRunner;
    private PLT plt;
    //private SG_CalibrationSequence calib;


    //
    //    // Required Constructor. Good place to set up references to objects in the unity scene
    public SensAtt_PredExperiment(ExperimentRunner runner, RunnableDesign runnableDesign) : base(runner, runnableDesign) {
         myRunner = (SensAtt_PredRunner)runner;  //cast the generic runner to your custom type.
//        // GameObject myGameObject = myRunner.MyGameObject;  // get reference to gameObject stored in your custom runner
    }

    // Optional Pre-Experiment code. Useful for pre-experiment calibration and setup.
    protected override void PreMethod() {
        //start EEG!
        plt = myRunner.EmptyObject.GetComponent<PLT>();
        plt.PLTsend(0);
        plt.PLTsend(126);
    }


    // Optional Pre-Experiment code. Useful for pre-experiment instructions.
    protected override IEnumerator PreCoroutine() {
        myRunner.RingR.SetActive(false);
        myRunner.RingL.SetActive(false);
        bool fingerinstruction = true;
        myRunner.Instruction_experiment.SetActive(true); 
        //used to makes sure that the finger is on the indicator :)) 
        while (fingerinstruction)
        {
            if (Input.GetKeyDown(KeyCode.Space))
            {
                myRunner.Instruction_experiment.SetActive(false);
                fingerinstruction = false;
            }
            yield return null;
        }
        myRunner.RingR.SetActive(true);
        myRunner.RingL.SetActive(true);
    }


    // Optional Post-Experiment code. Useful for experiment debrief instructions.
    protected override IEnumerator PostCoroutine() {
//        myRunner.EndExp.SetActive(true);
//        yield return new WaitForSeconds(10);
//        myRunner.EndExp.SetActive(false);
        yield return null; //required for coroutine
    }


    // Optional Post-Experiment code.
    protected override void PostMethod() {
        File.Move(Application.dataPath + "/Response_training.txt", "C:/Gian/GG_SensAtt_Prediction/02Data/ID" + myRunner.EmptyObject.GetComponent<ID>().ID_string + "/00Behavioural/" + "Experiment_response" + "_ID" + myRunner.EmptyObject.GetComponent<ID>().ID_string + ".txt");
        File.Move(Application.dataPath + "/Log.txt", "C:/Gian/GG_SensAtt_Prediction/02Data/ID" + myRunner.EmptyObject.GetComponent<ID>().ID_string + "/00Behavioural/" + "Experiment_Log" + "_ID" + myRunner.EmptyObject.GetComponent<ID>().ID_string + ".txt");
        File.Move(Application.dataPath + "/HandPositions.txt", "C:/Gian/GG_SensAtt_Prediction/02Data/ID" + myRunner.EmptyObject.GetComponent<ID>().ID_string + "/00Behavioural/" + "Experiment_Hand_Positions" + "_ID" + myRunner.EmptyObject.GetComponent<ID>().ID_string + ".txt");
        plt.PLTsend(127);
        // cleanup code (happens all in one frame)
    }
}

