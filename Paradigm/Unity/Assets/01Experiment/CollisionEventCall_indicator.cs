using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CollisionEventCall_indicator : MonoBehaviour
{
    public bool Indicator_touch = false; 
    // Start is called before the first frame update
    void Start()
    {
        
    }

    void OnTriggerEnter(Collider collision)
    {
        //        if (collision.name == "Thumb-Feedback") {
        //            Debug.Log("thumb visual" + Time.time*1000);
        //            File.AppendAllText(path, "Thumb_visual" + "\t" + Time.time*1000 + "\n");
        //            touch = true;
        //            plt.PLTsend(Thumb_visual_PLT);
        //            }
        if (collision.name == "Index-Feedback" || collision.name == "Hand" || collision.name == "r_index_finger_tip_marker") 
        {
            Indicator_touch = true; 
        }
    }
    void OnTriggerExit(Collider collision)
    {
        if(collision.name == "Index-Feedback" || collision.name == "Hand" || collision.name == "r_index_finger_tip_marker")
        {
            Indicator_touch = false; 
        }
    
    }
}
