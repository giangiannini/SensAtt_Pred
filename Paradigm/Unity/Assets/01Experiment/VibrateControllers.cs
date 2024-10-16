using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.XR;


public class VibrateControllers : MonoBehaviour
{
    public bool vibrationOn = false;
    public float intensity = 5.0f;
    public float frequency = 20.0f; 
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        if (vibrationOn == true)
        {
            OVRInput.SetControllerVibration(frequency, intensity, OVRInput.Controller.LTouch);
        }
        else if (vibrationOn == false)
        {
            OVRInput.SetControllerVibration(0, 0, OVRInput.Controller.LTouch);
        }
    }
}
