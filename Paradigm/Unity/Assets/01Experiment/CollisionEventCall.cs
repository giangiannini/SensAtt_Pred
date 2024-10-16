using System.Collections;
using System.Data;
using UnityEngine;
using UnityEngine.Events;
using bmlTUX.Scripts.ExperimentParts;
using bmlTUX.Scripts.Managers;
using System;
using System.Collections.Generic;
using System.IO;
using SG.Util;
using System.Runtime.InteropServices;
using SGCore.Haptics;
using UnityEngine.XR;


public class CollisionEventCall : MonoBehaviour
{
    //Set up variable path for printing log and touch for signalling end trial
    public bool touch;
    private string path;
    public UnityEvent call_visual;
    //public UnityEvent call_haptic;
    public ProvaBuzz buzzing;
    public VibrateControllers buzzing1; 
    public bool hapticActive = true;
    public bool forceActive = true;
    public bool freeze = false;
    public bool freeze1 = false;
    public PLT plt; 

    //Instance for PLT script
    //private PLT plt;

    //Set up plt int signals
//    public int Thumb_visual_PLT = 41;
//    public int Index_visual_PLT = 31;
//    public int Middle_visual_PLT = 43;
//    public int Ring_visual_PLT = 44;
//    public int Pinky_visual_PLT = 45;

    void Start() {

        //        path = Application.dataPath + "/Log.txt";
                plt = GameObject.Find("Scripts").GetComponent<PLT>();
    }
    void OnTriggerEnter(Collider collision) {

        //        if (collision.name == "Thumb-Feedback") {
        //            Debug.Log("thumb visual" + Time.time*1000);
        //            File.AppendAllText(path, "Thumb_visual" + "\t" + Time.time*1000 + "\n");
        //            touch = true;
        //            plt.PLTsend(Thumb_visual_PLT);
        //            }
        //if (collision.name == "Index-Feedback") {
        //    call_visual.Invoke();
        //    call_haptic.Invoke();
        //    touch = true;
        //    if (hapticActive == true){
        //        if (freeze == false) {
        //            StartCoroutine(TimedBuzz());
        //        }
        //    }
        //    if (forceActive == true){
        //        if (freeze1 == false) {
        //            StartCoroutine(TimedForce());
        //        }
        //    }
        //}

        if (collision.name == "r_index_finger_tip_marker")
        {

            call_visual.Invoke();
            touch = true; 
            if (hapticActive == true) 
            {
                if (freeze == false)
                {
                    //plt.PLTsend(124); 
                    plt.PLTsend(252); //this sends the proper electrical stimulation (128 to the stimulator and 124 to the amplifier). 
                    //StartCoroutine(TimedBuzz_oculusController());
                    StartCoroutine(TimedBuzz_electric()); // this scrips runs to avoid double stimulations in the same trial
                }
            }
            else if (hapticActive == false)
            {
                if (freeze == false)
                {
                    plt.PLTsend(125); // this simply sends the marker in the EEG
                    StartCoroutine(TimedBuzz_electric()); // this scrips runs to avoid double stimulations in the same trial
                }
            }

        }
    }

    IEnumerator TimedBuzz() {
        freeze = true;
        buzzing.bzzzzz = true;
        yield return new WaitForSeconds(0.15f);
        buzzing.bzzzzz = false;
        yield return new WaitForSeconds(1f);
        freeze = false;
        yield return new WaitForSeconds(0.2f);
        yield return null;
    }
    IEnumerator TimedForce() {
        freeze1 = true;
        buzzing.force = true;
        yield return new WaitForSeconds(1f);
        buzzing.force = false;
        yield return new WaitForSeconds(0.5f);
        freeze1 = false;
        yield return new WaitForSeconds(0.2f);
        yield return null;
    }

    IEnumerator TimedBuzz_oculusController()
    {
        freeze = true;
        buzzing1.vibrationOn = true; 
        yield return new WaitForSeconds(0.300f);
        buzzing1.vibrationOn = false;
        yield return new WaitForSeconds(1f);
        freeze = false;
        yield return new WaitForSeconds(0.2f);
        yield return null;
    }

    IEnumerator TimedBuzz_electric()
    {
        freeze = true;
        //yield return new WaitForSeconds(1.2f);
        //freeze = false;
        yield return null; 
    }
}
