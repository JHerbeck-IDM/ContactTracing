import time, os, shutil, sys

from create_asset_collection import create_asset_collection

from COMPS import Client
from COMPS.Data import WorkItem, WorkItemFile, QueryCriteria, Priority
from COMPS.Data.WorkItem import WorkItemState, WorkerOrPluginKey, RelationType

compshost = 'https://comps.idmod.org'
compsenv = 'Calculon'

if len(sys.argv) != 2:
    print('\r\nUsage:\r\n\t{0} C:\\path\\to\\create\\workitem\\from'.format(sys.argv[0]))
    exit()

path_to_wi = os.path.normpath(sys.argv[1])

if not os.path.exists(path_to_wi) or not os.path.isdir(path_to_wi):
    print('Path \'{0}\' doesn\'t exist or is not a directory'.format(path_to_wi))
    exit()

Client.login(compshost)

ac_id = None
path_to_ac = os.path.join(path_to_wi, 'Assets')

if os.path.exists(path_to_ac) and os.path.isdir(path_to_ac):
    print('Found \'Assets\' dir; creating AC')

    ac_id = create_asset_collection(path_to_ac)

# Create a work-item (locally)
wi = WorkItem('Create Singularity Image (HIV Contact Tracing)', WorkerOrPluginKey('ImageBuilderWorker', '1.0.0.0_RELEASE'), compsenv, asset_collection_id=ac_id, priority=Priority.AboveNormal)

tags = {
        'WorkItem type': 'ImageBuilderWorker'
}

wi.set_tags(tags)

workorder_string = \
"""{
  "WorkItem_Type": "ImageBuilderWorker",
  "Build": {
    "Type": "singularity",
    "Input": "Singularity.def",
    "Output": "HIV_contact_tracing.sif",
    "Tags": { }
  }
}"""

wi.add_work_order(data=bytes(workorder_string, 'utf-8'))

additional_files = [ ('input', f) for f in os.listdir(path_to_wi) if os.path.isfile(os.path.join(path_to_wi,f)) ]

# add the linked files
for tup in additional_files:
    print('Adding linked file: {0}'.format(tup[1]))
    wi.add_file(WorkItemFile(tup[1], tup[0], ''), file_path=os.path.join(path_to_wi, tup[1]))

# Save the work-item to the server
wi.save()

wi.refresh()

print('Created work-item {0}'.format(wi.id))
print('Commissioning...')

wi.commission()

print('Commissioned')
print('Refreshing WorkItem state until it completes')

print('State -> {}'.format(wi.state.name))

while wi.state not in [WorkItemState.Succeeded, WorkItemState.Failed, WorkItemState.Canceled]:
    time.sleep(5)
    wi.refresh()
    print('State -> {}'.format(wi.state.name))

print('Worker complete: ' + wi.state.name)

if wi.state is not WorkItemState.Succeeded:
    exit(-1)

# retrieve AC id
img_ac = wi.get_related_asset_collections()[0]
print(f'Created AC containing singularity image: {img_ac.id}')
