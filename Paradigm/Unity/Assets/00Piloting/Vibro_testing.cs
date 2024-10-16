using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.XR;


public class Vibro_testing : MonoBehaviour
{
    public VibrateControllers buzzing1;
    public PLT plt;
    public bool start = false;
    public bool freeze = false;
    public int signal = 1; 

    // Start is called before the first frame update
    void Start()
    {
        //int num = rnd.Next();


        //plt.PLTsend(126); 
    }

    // Update is called once per frame
    void Update()
    {
        if (start == true)
        {
            if (freeze == false)
            {
                freeze = true; 
                StartCoroutine(Vibration_interval());
            }
        }
    }

    IEnumerator Vibration_interval()
    {
        int iterationCount = 20; 
        for (int i = 0; i < iterationCount; i++)
        {
            // Code to be repeated.
            //plt.PLTsend(20);
            //buzzing1.vibrationOn = true;
            //yield return new WaitForSeconds(0.300f);
            //buzzing1.vibrationOn = false;
            plt.PLTsend(signal); 
            yield return new WaitForSeconds(1.5f);
            yield return new WaitForSeconds(Random.Range(0.0f, 1.0f));
            yield return null;
        }
    }
}
