using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using OVRTouchSample;
using System;
using UnityEngine.SceneManagement;

public class KeepActive : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        //OVRInput.SetControllerVibration(100, 100, OVRInput.Controller.RTouch);
    }

    // Update is called once per frame
    void Update()
    {

        Debug.Log(OVRPlugin.systemDisplayFrequency);
        //StartCoroutine(abcd());
    }

    IEnumerator abcd() {
        OVRInput.SetControllerVibration(1, 1, OVRInput.Controller.RTouch);
        yield return new WaitForSeconds(1);
        OVRInput.SetControllerVibration(0, 0, OVRInput.Controller.RTouch);
        yield return new WaitForSeconds(100);
    }
}
