using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.XR;
using UnityEngine.Experimental.XR;
using System;
using System.Linq;
using Settings = UnityEngine.XR.XRSettings;
using Node = UnityEngine.XR.XRNode;


public class ControllerStatus : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        Debug.Log(OVRInput.GetConnectedControllers());
    }

    public void CheckForControllers()
    {

//        bool xrControllerFound = false;
//        bool leftFound = false;
//        bool rightFound = false;
//
//        List<InputDevice> leftHandDevices = new List<InputDevice>();
//        InputDevices.GetDevicesWithCharacteristics(InputDeviceCharacteristics.Left, leftHandDevices);
//
//        if (leftHandDevices.Count > 0)
//        {
//            xrControllerFound = true;
//            Debug.Log(leftHandDevices);
//            leftFound = true;
//        }
//
//        List<InputDevice> rightHandDevices = new List<InputDevice>();
//        InputDevices.GetDevicesWithCharacteristics(InputDeviceCharacteristics.Right, rightHandDevices);
//        if (rightHandDevices.Count > 0)
//        {
//            xrControllerFound = true;
//                        Debug.Log(rightHandDevices);
//
//            rightFound = true;
//        }
//
//        if (!xrControllerFound)
//        {
//           // deal with this
//        }
//
////        leftController.SetActive(leftFound);
////        rightController.SetActive(rightFound);
    }
}
