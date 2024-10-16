using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;
using SG;

public class Reset_Calib : MonoBehaviour
{
    public SG_CalibrationSequence calibrationSequence;
    public GameObject IndexCollider;

   void Update() {
        if(Input.GetKeyDown(KeyCode.N)) {
            calibrationSequence = GameObject.Find("CalibrationLayer").GetComponent<SG_CalibrationSequence>();
            calibrationSequence.StartCalibration(true);
        }

        if(Input.GetKeyDown(KeyCode.O)) {
            Debug.Log(IndexCollider.transform.position.x);
            Debug.Log(IndexCollider.transform.position.y);
            Debug.Log(IndexCollider.transform.position.z);
        }
//        if(Input.GetKeyDown(KeyCode.I)) {
//                if(BlankScreen.activeInHierarchy){
//                    canvas.SetActive(false);
//                    Debug.Log("Instruction stop");
//                    File.AppendAllText(GetComponent<CreateLogOnStart>().path, "Instruction end" + "," + Time.time*1000 + "\n");
//                }
//                else {
//                    canvas.SetActive(true);
//                    Debug.Log("Instruction start");
//                    File.AppendAllText(GetComponent<CreateLogOnStart>().path, "Instruction start" + "," + Time.time*1000 + "\n");
//                }
//        }
   }

}
