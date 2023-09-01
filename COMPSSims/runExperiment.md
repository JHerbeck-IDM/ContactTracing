# How to Run

This is for cloning Adam's Eswatini sims and rerunning them with additional output.
It is based on Jonathan's work here: https://github.com/InstituteforDiseaseModeling/emodpy-tbhiv/blob/fancier-clone-and-run/examples/clone_and_run_by_files/example.py
This probably will not work right out of the box. Specifically it will have trouble
finding the Eradication executable. It will either need to be renamed in the
runExperiment.py script or I created a symbolic link to give it the right name.

```
python3 example.py <experiment_id>
```

Experiment ids known to work well:

* exp_id = '78a2a1ad-7643-ed11-a9fc-b88303911bc1' # recent, bloedow
* exp_id = 'bf11ad69-9e42-ed11-a9fc-b88303911bc1' # TBHIV_SIM (India-like)
* exp_id = '4739c9b7-cd26-ec11-9ecd-9440c9bee941' # HIV_SIM
* exp_id = '6b09923f-0844-ed11-a9fc-b88303911bc1' # Malaria_Sim

