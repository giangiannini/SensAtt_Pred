using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class KeepPosition : MonoBehaviour
{
    public GameObject RightHandAnchor;

    private Vector3 InitPosition = new Vector3(0f, 0f, -0.055f);
    private Vector3 InitRotation = new Vector3(0, 0, 0);

    public Vector3 LastPosition;
    public Quaternion LastRotation;

    // Start is called before the first frame update
    void Start()
    {

    }

    // Update is called once per frame
    void Update()
    {
        // First save the last position of our hand before it fucks up
        if (RightHandAnchor.transform.localPosition != InitPosition) {
            LastPosition = RightHandAnchor.transform.localPosition;
            LastRotation = RightHandAnchor.transform.localRotation;
            //Debug.Log(LastPosition);
        }
        if (RightHandAnchor.transform.localPosition == InitPosition) {
            RightHandAnchor.transform.localPosition = LastPosition;
            RightHandAnchor.transform.localRotation = LastRotation;
            //Debug.Log(LastPosition);
        }

    }
}
