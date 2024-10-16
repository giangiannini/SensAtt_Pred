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
/// Classes that inherit from Trial define custom behaviour for your experiment's trials.
/// Most experiments will need to edit this file to describe what happens in a trial.
///
/// This template shows how to set up a custom trial script using the toolkit's built-in functions.
///
/// You can delete any unused methods and unwanted comments. The only required parts are the constructor and the MainCoroutine.
/// </summary>
public class SensAtt_PredTrial : Trial {

    //Reference to my experiment runner
    SensAtt_PredRunner myRunner;

    //set up my variable path (reference from the CreateLogOnStart function).
    public string path;
    public string path_responses; 
    //set up my variable path1 (in which we will print our subject's response based on size)
    //public string path1 = Application.dataPath + "/Subj_Response.txt";

    // Required Constructor. Good place to set up references to objects in the unity scene
    public SensAtt_PredTrial(ExperimentRunner runner, DataRow data) : base(runner, data) {
        myRunner = (SensAtt_PredRunner)runner;  //cast the generic runner to your custom type.
    }

    public UnityEvent Caller_Haptic;
    public UnityEvent Caller_Vision;

    //Instance for touch bool (NOTE: private! I declare that my script is private and cannot directly change it. At least read it :))
    private CollisionEventCall CEC;
    //Instance for PLT script
    public PLT plt;
    //Instance for FingerFeedback script
    //public SG_FingerFeedback fingerfeedback;

    public int StartTrial = 10;
    public int preGO = 20;
    public int GO = 30;
    public int Index_visual_PLT = 128;

    // Stuff for vibration
    //public SG_HandSection handLocation = SG_HandSection.Index;
    //public SG_HapticGlove linkedGlove;
    //public SG_HandFeedback feedbackScript;

    //public int impactLevel1 = 1000;
    //public float vibrationTime1 = 0.06f;

    //public float speed = 1f;

    private Vector3 InitPositionRight = new Vector3(-0.221f, 0.7496f, 0.005f);
    private Vector3 InitPositionLeft = new Vector3(0.1786f, 0.7496f, 0.005f);
    private Vector3 TargetPositionRight = new Vector3(-0.221f, 0.7496f, 0.005f);
    private Vector3 TargetPositionLeft = new Vector3(0.1786f, 0.7496f, 0.005f);
    private Vector3 TargetPosition = new Vector3(0, 0, 0);

    // Optional Pre-Trial code. Useful for setting unity scene for trials. Executes in one frame at the start of each trial
    protected override void PreMethod() {
        // Collision Event Call reference
        CEC = myRunner.Visual.GetComponent<CollisionEventCall>();
        //Set to false at the start of each trial, it detects collision from another script and it uses it to determine when to finish trial
        CEC.touch = false;

        // Also from CEC get the event that gets called whether the finger touch the proper visual object
        UnityEvent Caller_Visual = CEC.call_visual;
        Caller_Visual.RemoveAllListeners();
        Caller_Visual.AddListener(caller_visual);

        // Parallel Port script (PLT) reference
        plt = myRunner.EmptyObject.GetComponent<PLT>();

        // Print something in the Exp files
        // Set path for Log file
        path = Application.dataPath + "/Log.txt";
        if ((int)Data["Trial"] == 0) {
            //Append to Log file some subjects data, experiment start infos and column titles
            File.AppendAllText(path, (string)Data["ID"] + "\t" + (string)Data["initials"] + "\n");
            File.AppendAllText(path, "Experiment Start" + "\t" + Time.time*1000 + "\n");
            File.AppendAllText(path, "Block \t Trial \t BlockType \t TrialType \t Stimulation \t Event_Name \t Time \n");

            myRunner.EmptyObject.GetComponent<ID>().ID_string = (string)Data["ID"];
        }

        //Make object and force/haptic field appear
        // myRunner.Visual.SetActive(true);
        // myRunner.Visual.GetComponent<MeshRenderer>().enabled = true;
        // myRunner.Haptic.SetActive(true);

        //Make object appear in the right position
        if ((int) Data["Trial_type"] == 1) {
            myRunner.Visual.transform.position = InitPositionRight;
            myRunner.Visual.SetActive(true);
            myRunner.RingL.SetActive(true);
        } else if ((int) Data["Trial_type"] == 2) {
            myRunner.Visual.transform.position = InitPositionRight;
            myRunner.Visual.SetActive(true);
            myRunner.RingL.SetActive(true);
            TargetPosition = TargetPositionLeft;
        } else if ((int) Data["Trial_type"] == 3) {
            myRunner.Visual.transform.position = InitPositionLeft;
            myRunner.Visual.SetActive(true);
            myRunner.RingR.SetActive(true);
        } else if ((int) Data["Trial_type"] == 4) {
            myRunner.Visual.transform.position = InitPositionLeft;
            myRunner.Visual.SetActive(true);
            myRunner.RingR.SetActive(true);
            TargetPosition = TargetPositionRight;
        }

        //Sets vibration on and off based on the trial
        if ((int)Data["Stimulation"] == 1) {
            CEC.hapticActive = true;
        }else if ((int)Data["Stimulation"] == 0) {
            CEC.hapticActive = false;
        }

        //After everything's ready, print start trial inside the logfile and send input to eeg.
        File.AppendAllText(path, "Start_trial_" + (int)Data["TrialInBlock"] + "\t" + Time.time*1000 + "\n");
        plt.PLTsend(StartTrial);
        myRunner.Fix.SetActive(true);


    }


    // Optional Pre-Trial code. Useful for waiting for the participant to
    // do something before each trial (multiple frames). Also might be useful for fixation points etc.
    protected override IEnumerator PreCoroutine()
    {
        bool keepfinger = true;
        while (keepfinger) {
            if (myRunner.RingR.GetComponent<CollisionEventCall_indicator>().Indicator_touch == true || myRunner.RingL.GetComponent<CollisionEventCall_indicator>().Indicator_touch == true)
            {
                myRunner.Instruction_IndicatorReminder.SetActive(false);
                myRunner.Visual.GetComponent<Renderer>().enabled = true;
                yield return new WaitForSeconds(1);
                // print in log "preGO"
                File.AppendAllText(path, (int)Data["Block"] + "\t" + (int)Data["TrialInBlock"] + "\t" + (string)Data["Block_type"] + "\t" + (int)Data["Trial_type"] + "\t" + (int)Data["Stimulation"] + "\t" + "preGO" + "\t" + Time.time * 1000 + "\n");
                plt.PLTsend(preGO); 

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
                File.AppendAllText(path, (int)Data["Block"] + "\t" + (int)Data["TrialInBlock"] + "\t" + (string)Data["Block_type"] + "\t" + (int)Data["Trial_type"] + "\t" + (int)Data["Stimulation"] + "\t" + "GO" + "\t" + Time.time * 1000 + "\n");
                plt.PLTsend(GO); 

                keepfinger = false; //required for coroutine
            }
            else if (myRunner.RingR.GetComponent<CollisionEventCall_indicator>().Indicator_touch == false || myRunner.RingR.GetComponent<CollisionEventCall_indicator>().Indicator_touch == false)
            {
                myRunner.Instruction_IndicatorReminder.SetActive(true);
            }
            yield return null;
        }
    }


    // Main Trial Execution Code.
    protected override IEnumerator RunMainCoroutine()
    {
        CEC.freeze = false;

        bool waitingForParticipantResponse = true;
        while (waitingForParticipantResponse)
        {
            // Ball moves towards the finger
            // Trial end conditions handling
            if ((int)Data["Trial_type"] == 2 || (int)Data["Trial_type"] == 4)
            {
                if (myRunner.RingR.GetComponent<CollisionEventCall_indicator>().Indicator_touch == true || myRunner.RingL.GetComponent<CollisionEventCall_indicator>().Indicator_touch == true)
                {
                    myRunner.Visual.transform.position = Vector3.MoveTowards(myRunner.Visual.transform.position, TargetPosition, (float)Data["Speed"] * Time.deltaTime);
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


    // Optional Post-Trial code. Useful for waiting for the participant to do something after each trial (multiple frames)
    protected override IEnumerator PostCoroutine() {
        // wait for ITI
        yield return new WaitForSeconds((float)Data["ITIs"]);
        yield return null;
    }


    // Optional Post-Trial code. useful for writing data to dependent variables and for resetting everything.
    // Executes in a single frame at the end of each trial
    protected override void PostMethod() {
    }



    //Functions called when UnityEvent associated to Haptic Feedback signalling is sent.
//    void caller_haptic()
//    {
//        //plt.PLTsend(Index_haptic_PLT);
//        if ((string) Data["Block_type"] == "TouchVision"){
////            this.linkedGlove = this.feedbackScript.TrackedHand.gloveHardware;
//////            SGCore.Finger finger1 = (SGCore.Finger)handLocation; //can do this since the finger indices match
//////            SGCore.Haptics.SG_TimedBuzzCmd buzzCmd = new SGCore.Haptics.SG_TimedBuzzCmd(finger1, impactLevel1, vibrationTime1);
////            //linkedGlove.SendCmd(buzzCmd);
//            Debug.Log("index haptic" + Time.time*1000);
//            File.AppendAllText(path, "Index_haptic" + "\t" + Time.time*1000 + "\n");
//        } if ((string) Data["Block_type"] == "TouchnoVision"){
////            this.linkedGlove = this.feedbackScript.TrackedHand.gloveHardware;
////            SGCore.Finger finger1 = (SGCore.Finger)handLocation; //can do this since the finger indices match
////            SGCore.Haptics.SG_TimedBuzzCmd buzzCmd = new SGCore.Haptics.SG_TimedBuzzCmd(finger1, impactLevel1, vibrationTime1);
////            linkedGlove.SendCmd(buzzCmd);
//            Debug.Log("index haptic" + Time.time*1000);
//            File.AppendAllText(path, "Index_haptic" + "\t" + Time.time*1000 + "\n");
//        }
//    }

    void caller_visual()
    {
        Debug.Log("index visual" + Time.time * 1000);
        File.AppendAllText(path, (int)Data["Block"] + "\t" + (int)Data["TrialInBlock"] + "\t" + (string)Data["Block_type"] + "\t" + (int)Data["Trial_type"] + "\t" + (int)Data["Stimulation"] + "\t" + "Index_visual" + "\t" + Time.time * 1000 + "\n");
        //plt.PLTsend(Index_visual_PLT);
    }
}