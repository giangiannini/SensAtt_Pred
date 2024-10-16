using System.Collections;
using System.Data;
using UnityEngine;
using UnityEngine.Events;
using bmlTUX.Scripts.ExperimentParts;
using bmlTUX.Scripts.Managers;
using SG;
using System;
using System.Collections.Generic;
using System.IO;
using SG.Util;
using System.Runtime.InteropServices;

/// <summary>
/// Classes that inherit from Trial define custom behaviour for your experiment's trials.
/// Most experiments will need to edit this file to describe what happens in a trial.
///
/// This template shows how to set up a custom trial script using the toolkit's built-in functions.
///
/// You can delete any unused methods and unwanted comments. The only required parts are the constructor and the MainCoroutine.
/// </summary>
public class SensAtt_Pred_trainingTrial : Trial {

    //Reference to my experiment runner
    SensAtt_Pred_trainingRunner myRunner;

    //set up my variable path (reference from the CreateLogOnStart function).
    public string path_training;
    //set up my variable path1 (in which we will print our subject's response based on size)
    //public string path1 = Application.dataPath + "/Subj_Response.txt";

    // Required Constructor. Good place to set up references to objects in the unity scene
    public SensAtt_Pred_trainingTrial(ExperimentRunner runner, DataRow data) : base(runner, data)
    {
        myRunner = (SensAtt_Pred_trainingRunner)runner;  //cast the generic runner to your custom type.
    }

    //Instance for touch bool (NOTE: private! I declare that my script is private and cannot directly change it. At least read it :))
    private CollisionEventCall CEC;

    public PLT plt;

    public int Index_visual_PLT = 128;

    public float speed = 0.5f;

    private Vector3 InitPositionRight = new Vector3(-0.221f, 0.7496f, 0.005f);
    private Vector3 InitPositionLeft = new Vector3(0.1786f, 0.7496f, 0.005f);
    private Vector3 TargetPositionRight = new Vector3(-0.221f, 0.7496f, 0.005f);
    private Vector3 TargetPositionLeft = new Vector3(0.1786f, 0.7496f, 0.005f);
    private Vector3 TargetPosition = new Vector3(0, 0, 0);

    // Optional Pre-Trial code. Useful for setting unity scene for trials. Executes in one frame at the start of each trial
    protected override void PreMethod() {
        path_training = Application.dataPath + "/Log_training.txt";
        if ((int)Data["Trial"] == 0)
        {
            File.AppendAllText(path_training, (string)Data["ID"] + "\t" + (string)Data["initials"] + "\n");
            File.AppendAllText(path_training, "Training Start" + "\t" + Time.time * 1000 + "\n");
            File.AppendAllText(path_training, "Block \t Trial \t BlockType \t TrialType \t Stimulation \t Event_Name \t Time \n");

            myRunner.EmptyObject.GetComponent<ID>().ID_string = (string)Data["ID"]; 
        }
        // Collision Event Call reference
        CEC = myRunner.Visual.GetComponent<CollisionEventCall>();
        //Set to false at the start of each trial, it detects collision from another script and it uses it to determine when to finish trial
        CEC.touch = false;

        // Also from CEC get the event that gets called whether the finger touch the proper visual object
        UnityEvent Caller_Visual = CEC.call_visual;
        Caller_Visual.RemoveAllListeners();
        Caller_Visual.AddListener(caller_visual);

        plt = myRunner.EmptyObject.GetComponent<PLT>();


        // Also from CEC get the event that gets called whether the finger touch the proper visual object. This time it
        // is done for the haptic feedback administration.
        //UnityEvent Caller_Haptic = CEC.call_haptic;
        //Caller_Haptic.RemoveAllListeners();
        //Caller_Haptic.AddListener(caller_haptic);
        // also add some stuff for vibration signal sending

        //Set active the right set of instructions!
        if ((int)Data["Block"] == 0)
        {
            // when initial intructions, nothing should be on the screen, just the indicators. 
            // therefore, I deactivate everything if "Instructions" == 1 (initial instructions)
            myRunner.Calibration_Instr.SetActive(true);
            myRunner.Visual.SetActive(false);
            myRunner.Haptic.SetActive(false);
        }
        else
        {   
            //Make object appear in the right position if instructions are not the initial ones. 
            //I will use the normal experiment script from this point onward. 
            if ((int)Data["Trial_type"] == 1)
            {
                myRunner.Visual.transform.position = InitPositionRight;
                myRunner.Visual.SetActive(true);
                //myRunner.Haptic.SetActive(true);
                myRunner.RingL.SetActive(true);
            }
            else if ((int)Data["Trial_type"] == 2)
            {
                myRunner.Visual.transform.position = InitPositionRight;
                myRunner.Visual.SetActive(true);
                //myRunner.Haptic.SetActive(true);
                myRunner.RingL.SetActive(true);
                TargetPosition = TargetPositionLeft;
            }
            else if ((int)Data["Trial_type"] == 3)
            {
                myRunner.Visual.transform.position = InitPositionLeft;
                myRunner.Visual.SetActive(true);
                //myRunner.Haptic.SetActive(true);
                myRunner.RingR.SetActive(true);
            }
            else if ((int)Data["Trial_type"] == 4)
            {
                myRunner.Visual.transform.position = InitPositionLeft;
                myRunner.Visual.SetActive(true);
                //myRunner.Haptic.SetActive(true);
                myRunner.RingR.SetActive(true);
                TargetPosition = TargetPositionRight;
            }

            //This turns on the instructions.
            if ((int)Data["Block"] == 1 && (int)Data["TrialInBlock"] == 0)
            {
                myRunner.CueMove_instr.SetActive(true);
                myRunner.Hand.SetActive(false);
            }
            else if ((int)Data["Block"] == 2 && (int)Data["TrialInBlock"] == 0)
            {
                myRunner.CueStay_instr.SetActive(true);
                myRunner.Hand.SetActive(false);
            }
        }

        //Sets vibration on and off based on the trial
        if ((int)Data["Stimulation"] == 1)
        {
            CEC.hapticActive = true;
        }
        else if ((int)Data["Stimulation"] == 0)
        {
            CEC.hapticActive = false;
        }

        //Set fixation cross active (re-sets active every trial, and stays on for trials of the same type, otherwise is switched off by block scrips
        // in order to have a correct image presentation). 
        myRunner.Fix.SetActive(true);
    }


    // Optional Pre-Trial code. Useful for waiting for the participant to
    // do something before each trial (multiple frames). Also might be useful for fixation points etc.
    protected override IEnumerator PreCoroutine() {
        if ((int)Data["Block"] == 0)
        {
            bool keepfinger = true; 
            while(keepfinger)
            {
                if (Input.GetKeyDown(KeyCode.Space))
                {
                    myRunner.Calibration_Instr.SetActive(false); 
                    keepfinger = false; 
                }
                yield return null;
            }
            
        } else {
            if ((int)Data["TrialInBlock"] == 0)
            {
                if ((int)Data["Block"] == 1 || (int)Data["Block"] == 2)
                {
                    bool keepfinger1 = true;
                    while (keepfinger1)
                    {
                        if (Input.GetKeyDown(KeyCode.Space))
                        {
                            myRunner.CueMove_instr.SetActive(false);
                            myRunner.CueStay_instr.SetActive(false);
                            myRunner.Sequence_instr.SetActive(false);
                            myRunner.Hand.SetActive(true);
                            keepfinger1 = false;
                        }
                        yield return null;
                    }
                }
            }
        
            bool keepfinger = true;

             
            while (keepfinger)
            {
                if (myRunner.RingR.GetComponent<CollisionEventCall_indicator>().Indicator_touch == true || myRunner.RingL.GetComponent<CollisionEventCall_indicator>().Indicator_touch == true)
                {
                    myRunner.Instruction_IndicatorReminder.SetActive(false);
                    myRunner.Visual.GetComponent<Renderer>().enabled = true;
                    yield return new WaitForSeconds(1);
                    // print in log "preGO"
                    File.AppendAllText(path_training, (int)Data["Block"] + "\t" + (int)Data["TrialInBlock"] + "\t" + (string)Data["Block_type"] + "\t" + (int)Data["Trial_type"] + "\t" + (int)Data["Stimulation"] + "\t" + "preGO" + "\t" + Time.time * 1000 + "\n");

                    if ((int)Data["Trial_type"] == 4 || (int)Data["Trial_type"] == 2)
                    {
                        myRunner.PreGO_Stay.SetActive(true); 
                    }
                    else if ((int)Data["Trial_type"] == 1 || (int)Data["Trial_type"] == 3)
                    {
                        myRunner.PreGO_Move.SetActive(true);
                    }
                    yield return new WaitForSeconds(0.5f); // time for the pre cue to flash of the right color;
                    myRunner.PreGO_Move.SetActive(false);
                    myRunner.PreGO_Stay.SetActive(false); 

                    // print in log "GO"
                    File.AppendAllText(path_training, (int)Data["Block"] + "\t" + (int)Data["TrialInBlock"] + "\t" + (string)Data["Block_type"] + "\t" + (int)Data["Trial_type"] + "\t" + (int)Data["Stimulation"] + "\t" + "GO" + "\t" + Time.time * 1000 + "\n");

                    keepfinger = false; //required for coroutine
                }
                else if (myRunner.RingR.GetComponent<CollisionEventCall_indicator>().Indicator_touch == false || myRunner.RingR.GetComponent<CollisionEventCall_indicator>().Indicator_touch == false)
                {
                    myRunner.Instruction_IndicatorReminder.SetActive(true);
                }
                yield return null;
            }
        }
    }


    // Main Trial Execution Code.
    protected override IEnumerator RunMainCoroutine()
    {
        CEC.freeze = false;

        if ((int)Data["Block"] == 0)
        {
            myRunner.Hand.SetActive(false);
            bool waitingForParticipantResponse = true;
            myRunner.Training_initial_Instr1.SetActive(true);
            while (waitingForParticipantResponse)
            {
                if (Input.GetKeyDown(KeyCode.Space))
                {
                    waitingForParticipantResponse = false;
                }
            yield return null;
            }
            bool waitingForParticipantResponse1 = true;
            myRunner.Training_initial_Instr2.SetActive(true);
            while (waitingForParticipantResponse1)
            {
                if (Input.GetKeyDown(KeyCode.Space))
                {
                    myRunner.Training_initial_Instr2.SetActive(false);
                    myRunner.Training_initial_Instr1.SetActive(false);
                    waitingForParticipantResponse1 = false;
                }
                yield return null;
            }
            myRunner.Hand.SetActive(true); 
        }
        else
        {
            bool waitingForParticipantResponse = true;
            while (waitingForParticipantResponse)
            {
                // Ball moves towards the finger
                // Trial end conditions handling
                if ((int)Data["Trial_type"] == 2 || (int)Data["Trial_type"] == 4)
                {
                    if (myRunner.RingR.GetComponent<CollisionEventCall_indicator>().Indicator_touch == true || myRunner.RingL.GetComponent<CollisionEventCall_indicator>().Indicator_touch == true)
                    {
                        myRunner.Visual.transform.position = Vector3.MoveTowards(myRunner.Visual.transform.position, TargetPosition, speed * Time.deltaTime);
                    }
                    // also when in this condition, make sure that the finger doesn't move from the indicator
                    else //(myRunner.Indicator.GetComponent<CollisionEventCall_indicator>() == false)
                    {
                        myRunner.Visual.SetActive(false);
                        //myRunner.Haptic.SetActive(false); 
                        myRunner.Instruction_IndicatorReminder.SetActive(true);
                        yield return new WaitForSeconds(3);
                        myRunner.Instruction_IndicatorReminder.SetActive(false);
                        waitingForParticipantResponse = false;
                    }

                }

                //General trial end conditions handling (when participant touches somehow the ball)
                if (CEC.touch == true)
                {
                    yield return new WaitForSeconds(1);
                    myRunner.Visual.GetComponent<Renderer>().enabled = false;
                    CEC.freeze = true; 
                    waitingForParticipantResponse = false;

                }


                //condition for exiting trial by keypress
                if (Input.GetKeyDown(KeyCode.Space))
                {
                    waitingForParticipantResponse = false;
                }
                yield return null;
            }
        }
    }

    // Optional Post-Trial code. Useful for waiting for the participant to do something after each trial (multiple frames)
    protected override IEnumerator PostCoroutine() {
        yield return new WaitForSeconds((float)Data["ITIs"]);
        yield return null;
    }


    // Optional Post-Trial code. useful for writing data to dependent variables and for resetting everything.
    // Executes in a single frame at the end of each trial
    protected override void PostMethod() {
        // How to write results to dependent variables: 
        // Data["MyDependentFloatVariable"] = someFloatVariable;
    }




    ////Functions called when UnityEvent associated to Haptic Feedback signalling is sent.
    //void caller_haptic()
    //{
    //    //plt.PLTsend(Index_haptic_PLT);
    //    if ((string)Data["Block_type"] == "TouchVision")
    //    {
    //        //            this.linkedGlove = this.feedbackScript.TrackedHand.gloveHardware;
    //        ////            SGCore.Finger finger1 = (SGCore.Finger)handLocation; //can do this since the finger indices match
    //        ////            SGCore.Haptics.SG_TimedBuzzCmd buzzCmd = new SGCore.Haptics.SG_TimedBuzzCmd(finger1, impactLevel1, vibrationTime1);
    //        //            //linkedGlove.SendCmd(buzzCmd);
    //        Debug.Log("index haptic" + Time.time * 1000);
    //        File.AppendAllText(path, "Index_haptic" + "\t" + Time.time * 1000 + "\n");
    //    }
    //    if ((string)Data["Block_type"] == "TouchnoVision")
    //    {
    //        //            this.linkedGlove = this.feedbackScript.TrackedHand.gloveHardware;
    //        //            SGCore.Finger finger1 = (SGCore.Finger)handLocation; //can do this since the finger indices match
    //        //            SGCore.Haptics.SG_TimedBuzzCmd buzzCmd = new SGCore.Haptics.SG_TimedBuzzCmd(finger1, impactLevel1, vibrationTime1);
    //        //            linkedGlove.SendCmd(buzzCmd);
    //        Debug.Log("index haptic" + Time.time * 1000);
    //        File.AppendAllText(path, "Index_haptic" + "\t" + Time.time * 1000 + "\n");
    //    }
    //}

    void caller_visual()
    {
        //plt.PLTsend(Index_visual_PLT);
        Debug.Log("index visual" + Time.time * 1000);
        File.AppendAllText(path_training, (int)Data["Block"] + "\t" + (int)Data["TrialInBlock"] + "\t" + (string)Data["Block_type"] + "\t" + (int)Data["Trial_type"] + "\t" + (int)Data["Stimulation"] + "\t" + "Index_visual" + "\t" + Time.time * 1000 + "\n");
    }
}

