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

public class PLT : MonoBehaviour
{
    [DllImport("inpoutx64.dll")]
    private static extern UInt32 IsInpOutDriverOpen();

    //Inport the function Out32 that is the main function that we'll use to send signals to Parallel Port.
    [DllImport("inpoutx64.dll")]
    private static extern void Out32(int PortAddress, int Data);

    [DllImport("inpoutx64.dll", EntryPoint = "IsInpOutDriverOpen")]
    private static extern UInt32 IsInpOutDriverOpen_x64();

    [DllImport("inpoutx64.dll", EntryPoint = "Out32")]
    private static extern void Out32_x64(short Address, short Data);

    public int PortAddress = 16376;

    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    //Function to send the Parallel Port infos through

    //Declared public otherwise it would have had some problems being called from the main trial scripts
    public void PLTsend(int Data){
        StartCoroutine(LPTWrite(Data));
    }

    IEnumerator LPTWrite(int Data){
        Out32(PortAddress, Data);
        yield return new WaitForSeconds((float)0.001);
        Out32(PortAddress, 0);
    }

}
