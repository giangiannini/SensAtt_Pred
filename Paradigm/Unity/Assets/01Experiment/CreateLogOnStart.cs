using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using UnityEngine.Events;
using SG;
    
public class CreateLogOnStart : MonoBehaviour
{
    /// <summary>
    /// This small script will serve to set up the log file that is generated in the Assets folder :)
    /// The file here is only set up and will be updated trial by trial in the TUX_SGRunner script.
    /// </summary>
    public string path;
    public string path_responses;
    public SG_CalibrationSequence calibrationSequence;

    // Start is called before the first frame update
    public void Start() {
        path = Application.dataPath + "/Log.txt";
        path_responses = Application.dataPath + "Responses.txt";

        //UnityEvent prova = GameObject.Find("SG_Left_Hand").GetComponent<SG_HapticGlove>().CalibrationStateChanged;
        File.WriteAllText(path, "Log Start \n");
        File.AppendAllText(path, System.DateTime.Now + "\t" + Time.time*1000 + "\n");

        File.WriteAllText(path_responses, "Response Start \n");
        File.AppendAllText(path_responses, System.DateTime.Now + "\t" + Time.time * 1000 + "\n"); 
        //string titles = "Counter," + "Event," + "Time" + "\n";
        //File.AppendAllText(path, titles);

        //UnityEvent startcalib = calibrationSequence.CalibHasStarted;
        //UnityEvent endcalib = calibrationSequence.CalibHasFinished;
        //startcalib.AddListener(StartCalib);
        //endcalib.AddListener(EndCalib);

    }


    //void StartCalib()
    //{
    //    Debug.Log("Starting calibration");
    //    File.AppendAllText(path, "Calibration start" + "\t" + Time.time*1000 + "\n");
    //}
    
    //void EndCalib()
    //{
    //    Debug.Log("Finishing calibration");
    //    File.AppendAllText(path, "Calibration end" + "\t" + Time.time*1000 + "\n");
    //}

}
