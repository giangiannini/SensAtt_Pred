using bmlTUX.Scripts.ExperimentParts;
// ReSharper disable once RedundantUsingDirective
using UnityEngine;

/// <summary>
/// This class is the main communication between the toolkit and the Unity scene. Drag this script onto an empty gameObject in your Unity scene.
/// In the gameObject's inspector you need to drag in your design file and any custom scripts.
/// </summary>
public class SensAtt_Pred_trainingRunner : ExperimentRunner {
    // Here is where you make a list of objects in your unity scene that need to be referenced by your scripts.
    //public GameObject ReferenceToGameObject;
    public GameObject Visual;
    public GameObject Haptic;
    public GameObject RingL;
    public GameObject RingR;
    public GameObject Hand;
    public GameObject FingerTip; 

    // Empty object containing most of the scripts of the experiment
    public GameObject EmptyObject;

    // Stuff for experiment canvases
    public GameObject Instruction_experiment;
    public GameObject Instruction_experiment1;
    public GameObject Instruction_blockL;
    public GameObject Instruction_blockR;
    public GameObject Instruction_IndicatorReminder;
    public GameObject Instruction_calibration;
    public GameObject Warning;
    public GameObject Fix;
    public GameObject PreGO_Stay;
    public GameObject PreGO_Move; 
    public GameObject StartBlock;
    public GameObject Response;
    public GameObject Break;

    // Colors of the stimuli
    public Material SkinTone;
    public Material Transparent;
    public Material Black;
    public Material Green;
    public Material Red;

    // Stuff for response
    public GameObject Resp25;
    public GameObject Resp50;
    public GameObject Resp75;

    // Stuff for training phase
    public GameObject Calibration_Instr;
    public GameObject Training_initial_Instr1;
    public GameObject Training_initial_Instr2;
    public GameObject CueMove_instr;
    public GameObject CueStay_instr;
    public GameObject Sequence_instr;

}





