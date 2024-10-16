using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Data;
using UnityEngine.Events;
using bmlTUX.Scripts.ExperimentParts;
using SG;
using SG.Util;
using SGCore.Haptics;


public class ProvaBuzz : MonoBehaviour
{
    public SG.SG_HapticGlove glove;
    SGCore.HapticGlove hapticGlove = null;
    SGCore.Haptics.SG_BuzzCmd buzzCmd = SG_BuzzCmd.Off;
    SGCore.Haptics.SG_FFBCmd forcefeedCmd = SG_FFBCmd.Off;
    bool setup = false;
    public bool bzzzzz = false;
    public bool force = false;


    private void SetupAfterConnect()
    {
        hapticGlove = (SGCore.HapticGlove)glove.InternalGlove;
        Debug.Log(hapticGlove);
        if (hapticGlove == null)
        {
            Debug.Log("riprovo a connettermi al guanto");
        }
        else
        {
            setup = true;
        }
    }

    // Start is called before the first frame update
    void Start()
    {
        //hapticGlove = (SGCore.HapticGlove)glove.InternalGlove;
    }

    // Update is called once per frame
    void Update()
    {
        if (!setup)
        {
            this.SetupAfterConnect();
        }
        //SGCore.Haptics.SG_BuzzCmd buzz = new SGCore.Haptics.SG_BuzzCmd(SGCore.Finger.Index, 100);
        //SGCore.Haptics.SG_BuzzCmd.Off buzz = new SGCore.Haptics.SG_BuzzCmd(0,100,0,0,0);
        if (bzzzzz == true)
        {
            int[] buzz = new int[5];
            buzz[0] = 0;
            buzz[1] = 100;
            buzz[2] = 0;
            buzz[3] = 0;
            buzz[4] = 0;
            buzzCmd = new SG_BuzzCmd(buzz);
            Debug.Log("Bzzzz");
            hapticGlove.SendHaptics(buzzCmd);
        }
        else if (bzzzzz == false)
        {
            int[] buzz = new int[5];
            buzz[0] = 0;
            buzz[1] = 0;
            buzz[2] = 0;
            buzz[3] = 0;
            buzz[4] = 0;
            buzzCmd = new SG_BuzzCmd(buzz);
            hapticGlove.SendHaptics(buzzCmd);
        }
        if (force == true)
        {
            int[] forcefeed = new int[5];
            forcefeed[0] = 0;
            forcefeed[1] = 100;
            forcefeed[2] = 0;
            forcefeed[3] = 0;
            forcefeed[4] = 0;
            forcefeedCmd = new SG_FFBCmd(forcefeed);
            Debug.Log("Forceeee");
            hapticGlove.SendHaptics(forcefeedCmd);
        }
        else if (force == false)
        {
            int[] forcefeed = new int[5];
            forcefeed[0] = 0;
            forcefeed[1] = 0;
            forcefeed[2] = 0;
            forcefeed[3] = 0;
            forcefeed[4] = 0;
            forcefeedCmd = new SG_FFBCmd(forcefeed);
            hapticGlove.SendHaptics(forcefeedCmd);
        }
    }

    void OnApplicationQuit()
    {
        if (this.hapticGlove != null)
        {
            this.hapticGlove.StopHaptics(); //end Haptics(!)
        }
    }
}
