using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ResponseScript : MonoBehaviour
{
    public bool responding = false;
    public int response = 0;
    // Start is called before the first frame update
    void Start()
    {

    }

    IEnumerator OnTriggerEnter(Collider collision)
    {
        responding = true;
        yield return new WaitForSeconds(1);
        if (responding == true)
        {
            if (collision.name == "Resp25")
            {
                response = 25;
                Debug.Log(response);
            }
            else if (collision.name == "Resp50")
            {
                response = 50;
                Debug.Log(response);
            }
            else if (collision.name == "Resp75")
            {
                response = 75;
                Debug.Log(response);
            }
        }
        //StartCoroutine(CheckResponse());
    }

    void OnTriggerExit(Collider collision)
    {
        responding = false;
        response = 0;
        //Debug.Log(response); 
    }

}
