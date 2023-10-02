import sys
import pathlib # for a join
import shutil
import os
import threading
import copy

from COMPS import Client
from COMPS import Data
from COMPS.Data import Experiment, Simulation, QueryCriteria, AssetCollection, AssetCollectionFile

# emodpy
from emodpy.emod_task import EMODTask
import emodpy.emod_task as emod_task
from idmtools.core.platform_factory import Platform

import manifest

comps_host = 'https://comps.idmod.org'
Client.login(comps_host)
sweep_size = 1
threads = []

def get_sim_files_from_comps( sim_id ):
    import sys, os, shutil

    from COMPS import Client
    from COMPS.Data import Simulation, AssetCollection, QueryCriteria

    sim = Simulation.get(sim_id)

    print('Found sim {0}'.format(str(sim.id)))
    print(sim)

    # example of code to just re-run failed sims from an experiment
    #if "Failed" not in str(sim.state):
    #    return None

    sim.refresh(QueryCriteria().select_children(['files', 'hpc_jobs']))

    # clear the output dir if it already exists
    p = manifest.assets_dir
    sp = manifest.special_assets_dir
    if os.path.exists(p):
        shutil.rmtree(p, True)
    os.makedirs(p)
    #if os.path.exists(sp):
        #shutil.rmtree(sp, True)
    #os.makedirs(sp)

    print('Getting input files', end='')
    notes = { "SIF": None, "Python": None, "Config": None, "Name": sim.name }

    for f in sim.files:
        if f.file_type.lower() == 'input':
            with open(os.path.join(p, f.file_name), 'wb') as outfile:
                outfile.write(f.retrieve())
            print('.', end='')
            if "config" in f.file_name: 
                notes["Config"] = f.file_name
    print()

    job = sim.hpc_jobs[-1]

    if job and job.configuration:
        if job.configuration.asset_collection_id:
            print('Getting asset-collection files', end='')
            os.makedirs(os.path.join(p, 'Assets'))

            ac = AssetCollection.get(job.configuration.asset_collection_id, QueryCriteria().select_children(['assets']))

            for asset in ac.assets:
                # TBD: there are some special assets:
                # Eradication, SIF, and python go in SpecialAssets
                if "Eradication" in asset.file_name or asset.file_name.endswith( ".sif" ) or asset.relative_path=="python":
                    put_path = sp
                else:
                    put_path=p

                if asset.file_name.endswith( ".sif" ):
                    # Special case. Don't download file.
                    ac_id = get_sif_by_id( asset )
                    #notes["SIF"] = asset.file_name
                    with open(os.path.join(put_path, "sif.id"), 'w') as outfile:
                        outfile.write(str(ac_id))
                    notes["SIF"] = "sif.id"
                    continue
                if asset.relative_path == "python":
                    notes["Python"] = True
                #if "Eradication" in asset.file_name and asset.file_name.endswith( ".exe"):
                    #raise ValueError( f"Executable appears to be windows {asset.file_name}. Not yet supported." )

                ap = os.path.join(put_path, asset.relative_path or '')
                if not os.path.exists(ap):
                    os.makedirs(ap)
                try:
                    file_path = os.path.join(ap, asset.file_name)
                    if not os.path.exists( file_path ):
                        with open(file_path, 'wb') as outfile:
                            outfile.write(asset.retrieve())
                except Exception as ex:
                    print( str( ex ) )
                    import pdb
                    pdb.set_trace()
                print('.', end='')

            print()

        cmdline = job.configuration.executable_path + ' ' + job.configuration.simulation_input_args
    print( notes )
    return notes
#node_map = {}

def get_sif_by_id( sim_id ):

    s = Simulation.get(sim_id, QueryCriteria().select_children(['configuration']).add_extra_params({'coalesceconfig':True}))
    ac = AssetCollection.get(s.configuration.asset_collection_id, QueryCriteria().select_children(['assets']))
    print(ac.id)
    ## WARNING - hackery ahead... we don't have a direct connection with the original AC that just contained a SIF, so
    # instead we create a *new* AC that just contains the SIF file that we found from the merged AC, and we should get
    # back the same AC we had originally (since if two people "create" identical ACs, the second time, we just get back
    # the ID of the first AC rather than creating a new, identical AC)
    sif_asset = next(af for af in ac.assets if af.file_name.endswith('.sif'))
    new_ac = AssetCollection()
    new_ac.add_asset(AssetCollectionFile(file_name=sif_asset.file_name, md5_checksum=sif_asset.md5_checksum))
    new_ac.save()
    print(new_ac.id)
    return new_ac.id

def update_sim_random_seed(simulation, sim_id):
    simulation.task.config["Run_Number"] = sim_id
    return {"Original_Sim_Id": sim_id}

from idmtools.assets import Asset, AssetCollection  #
from idmtools.builders import SimulationBuilder
from idmtools.entities.experiment import Experiment

platform = Platform("Calculon", priority="Normal", num_retries=3)

def run_test_in_thread(experiment_name="EMOD Clone & Run", python_path=None, config_path=None, sif_path=None):
    # Create a platform
    #platform = Platform("SLURM", node_group="idm_a", priority="Highest")

    # WIP: To support use case of legacy used old windows binary and we want to use new linux binary
    #manifest.ep4 = "EP4"
    #manifest.sif="dtk_rocky_36_emodapi.id"
    # create EMODTask
    task = EMODTask.from_files(
            eradication_path=os.path.join( manifest.special_assets_dir, "Eradication" ),
            ep4_path=python_path,
            #custom_reports_path=manifest.plugins_folder,
            config_path=os.path.join( manifest.assets_dir, config_path )
        )

    task.common_assets.add_directory(assets_directory=manifest.assets_dir)

    print( sif_path )
    if sif_path:
        #task.set_sif( os.path.join( "SpecialAssets", real_manifest_sif ) )
        task.set_sif( os.path.join( sif_path ) )

    builder = SimulationBuilder()
    builder.add_sweep_definition( update_sim_random_seed, experiment_name )

    # create experiment from builder
    experiment  = Experiment.from_builder(builder, task, name=experiment_name)

    def wait_until_done( experiment, task ):
        # The last step is to call run() on the ExperimentManager to run the simulations.
        experiment.run(wait_until_done=True, platform=platform)

        # Check result
        task.handle_experiment_completion( experiment )
        EMODTask.get_file_from_comps( experiment.uid, [ "InsetChart.json" ] )
        task.cache_experiment_metadata_in_sql( experiment.uid )
        #import emod_api.channelreports.plot_icj_means as plotter
        #chan_data = plotter.collect( str( experiment.uid ), "Infected" )
        #plotter.display( chan_data, False, "Infected", str( experiment.uid ) )
    t = threading.Thread( target=wait_until_done, args=(experiment,task,) )
    threads.append(t)
    t.start()


def run_sim( sim_id ):
    notes = get_sim_files_from_comps( sim_id )
    real_manifest_sif = None # manifest.sif
    if notes["SIF"]:
        real_manifest_sif = os.path.join( manifest.special_assets_dir, notes["SIF"] )
    real_config = notes["Config"]
    real_python = manifest.ep4
    if not notes["Python"]:
        real_python = None
    print( real_manifest_sif )
    print( real_config )
    print( real_python )


    run_test_in_thread( str(sim_id), python_path=real_python, config_path=real_config, sif_path=real_manifest_sif )

def run_experiment_from_cloned_tasks( sims ):
    # Create an experiment
    experiment = Experiment(name='Eswatini with relationships')
    counter = 0
    for sim in sims:
        counter += 1
        # if counter == 10:
        #     break
        sim_id = sim.id
        notes = get_sim_files_from_comps( sim_id )
        if not notes:
            continue
        real_manifest_sif = None # manifest.sif
        if notes["SIF"]:
            real_manifest_sif = os.path.join( manifest.special_assets_dir, notes["SIF"] )
        real_config = notes["Config"]
        real_python = manifest.ep4
        if not notes["Python"]:
            real_python = None
        print( real_manifest_sif )
        print( real_config )
        print( real_python )
        task = EMODTask.from_files(
            eradication_path=os.path.join( manifest.special_assets_dir, "Eradication" ),
            ep4_path=real_python,
            #custom_reports_path=manifest.plugins_folder,
            config_path=os.path.join( manifest.assets_dir, real_config )
        )
        task.common_assets.add_directory(assets_directory=manifest.assets_dir)
        
        camp_json_name = "Assets/" + task.config["Campaign_Filename"]
        task.config["Campaign_Filename"] = camp_json_name
        if task.config["Custom_Reports_Filename"] != "":
            task.config["Custom_Reports_Filename"] = "Assets/" + task.config["Custom_Reports_Filename"]
            
        
        task.config["Report_Relationship_End"] = 1
        task.config["Report_Relationship_Start"] = 1
        task.config["Report_Relationship_Start_Individual_Properties"] = ["Risk"]
        task.config["Report_Event_Recorder"] = 1
        task.config["Report_Event_Recorder_Events"] = [
                "Births",
                "DiseaseDeaths",
                "NonDiseaseDeaths",
                "OpportunisticInfectionDeath",
                "NewInfectionEvent",
                "ARTStaging0",
                "ARTStaging1",
                "ARTStaging2",
                "ARTStaging3",
                "ARTStaging4",
                "ARTStaging5",
                "ARTStaging6",
                "ARTStaging8",
                "ARTStaging9",
                "HIV_Positive_at_ANC",
                "OnART0",
                "OnPreART0",
                "OnPreART4",
                "HIVInfectionStageEnteredLatent",
                "HIVInfectionStageEnteredAIDS",
                "HIVInfectionStageEnteredOnART",
        ]
        
        sim = task.to_simulation()
        sim.name = str(sim)
        sim.tags = {'sim_id':str(sim_id)}
        experiment.simulations.append(sim) 
   

        # run experiment
    experiment.run(platform=platform, wait_until_done=True)

exp_id = sys.argv[1]
print( f"Experiment ID = {exp_id}." )
exp = Data.Experiment.get(exp_id)
sims = [ sim for sim in exp.get_simulations() ]
print( f"Found {len(sims)} in experiment." )
run_experiment_from_cloned_tasks( sims )
