using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.IO;

public class PrintPositions : MonoBehaviour
{
    private string path2;
    public GameObject Index;

    // Start is called before the first frame update
    void Start()
    {
        path2 = Application.dataPath + "/HandPositions.txt";
        File.WriteAllText(path2, "Log Start \n");
        File.AppendAllText(path2, System.DateTime.Now + "\t" + Time.time*1000 + "\n");
        File.AppendAllText(path2, "Time (ms)" + "\t" + "Index" + "\t \t \t \t \t \t" + "\n");
    }

    // Update is called once per frame
    void Update()
    {
        File.AppendAllText(path2, Time.time*1000 + "\t");
        //File.AppendAllText(path2, Wrist.transform.position.x + "\t" + Wrist.transform.position.y + "\t" + Wrist.transform.position.z + "\t" + Wrist.transform.rotation.x + "\t" + Wrist.transform.rotation.y + "\t" + Wrist.transform.rotation.z + "\t");
        //File.AppendAllText(path2, Thumb.transform.position.x + "\t" + Thumb.transform.position.y + "\t" + Thumb.transform.position.z + "\t" + Thumb.transform.rotation.x + "\t" + Thumb.transform.rotation.y + "\t" + Thumb.transform.rotation.z + "\t");
        File.AppendAllText(path2, Index.transform.position.x + "\t" + Index.transform.position.y + "\t" + Index.transform.position.z + "\t" + Index.transform.rotation.x + "\t" + Index.transform.rotation.y + "\t" + Index.transform.rotation.z + "\n");
        //File.AppendAllText(path2, Middle.transform.position.x + "\t" + Middle.transform.position.y + "\t" + Middle.transform.position.z + "\t" + Middle.transform.rotation.x + "\t" + Middle.transform.rotation.y + "\t" + Middle.transform.rotation.z + "\t");
        //File.AppendAllText(path2, Ring.transform.position.x + "\t" + Ring.transform.position.y + "\t" + Ring.transform.position.z + "\t" + Ring.transform.rotation.x + "\t" + Ring.transform.rotation.y + "\t" + Ring.transform.rotation.z + "\t");
        //File.AppendAllText(path2, Pinky.transform.position.x + "\t" + Pinky.transform.position.y + "\t" + Pinky.transform.position.z + "\t" + Pinky.transform.rotation.x + "\t" + Pinky.transform.rotation.y + "\t" + Pinky.transform.rotation.z + "\n");
    }
}
