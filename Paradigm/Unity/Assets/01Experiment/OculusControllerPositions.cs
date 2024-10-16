using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.XR;

public class OculusControllerPositions : MonoBehaviour
{
    public GameObject VirtualHand;
    private Vector3 multiplicationVector = new Vector3(-1f, 0f, -1f);
    public Vector3 startingposition = new Vector3(0f, 0.77f, 0f);

    public Vector3 multiplicationVector1 = new Vector3(1, 1, 1);
    public Vector3 startingrotation = new Vector3(0, 90, 0);


    //public GameObject Hand; 
    // Start is called before the first frame update
    void Start()
    {
        //startingposition = VirtualHand.transform.position;
        //Debug.Log(startingposition); 
    }

    // Update is called once per frame
    void Update()
    {
        //Debug.Log(OVRInput.GetLocalControllerPosition(OVRInput.Controller.RTouch));
        VirtualHand.transform.position = startingposition + Vector3.Scale(OVRInput.GetLocalControllerPosition(OVRInput.Controller.RTouch), multiplicationVector);
        Quaternion prova = OVRInput.GetLocalControllerRotation(OVRInput.Controller.RTouch);
        VirtualHand.transform.rotation = new Quaternion(0, prova.y, 0, prova.w);
        //VirtualHand.transform.rotation = Quaternion.Euler(startingrotation) * OVRInput.GetLocalControllerRotation(OVRInput.Controller.RTouch) * Quaternion.Euler(multiplicationVector1);
        //Debug.Log(Quaternion.Euler(startingrotation));
        //VirtualHand.transform.eulerAngles = startingrotation + Vector3.Scale(multiplicationVector1, multiplicationVector1);

        //Hand.transform.position = startingposition + Vector3.Scale(OVRInput.GetLocalControllerPosition(OVRInput.Controller.RTouch), multiplicationVector);
        //Debug.Log(startingposition);
    }
}
