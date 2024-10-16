using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Active : MonoBehaviour
{
    public GameObject Visual;
    public GameObject Haptic;
    public GameObject EndExp;
    public GameObject BlankScreen;
    public GameObject Instruction;
    public GameObject StartBlock;
    public GameObject EndBlock;

    // The game starts and some objects must be set inactive otherwise they will pop up in the screen during calibration
    void Start()
    {
       Visual.SetActive(false);
       Haptic.SetActive(false);
       EndExp.SetActive(false);
       BlankScreen.SetActive(false);
       Instruction.SetActive(false);
       StartBlock.SetActive(false);
       EndBlock.SetActive(false);
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
