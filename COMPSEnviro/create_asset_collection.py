import os
import sys
from hashlib import md5

from COMPS import Client
from COMPS.Data import AssetCollection, AssetCollectionFile

compshost = 'https://comps.idmod.org'

def create_asset_collection(path_to_ac):
    path_to_ac = os.path.normpath(path_to_ac)

    if not os.path.exists(path_to_ac) or not os.path.isdir(path_to_ac):
        print('Path \'{0}\' doesn\'t exist or is not a directory'.format(path_to_ac))
        exit()

    tags = {
        'Example AssetCollection': None,
        'Name': os.path.basename(path_to_ac)
    }

    ac = AssetCollection()
    ac.set_tags(tags)

    # First try creating it without uploading any files (just by md5sum)
    for (dirpath, dirnames, filenames) in os.walk(path_to_ac):
        for fn in filenames:
            rp = os.path.relpath(dirpath, path_to_ac) if dirpath != path_to_ac else ''
            print('Adding {0}'.format(os.path.join(rp, fn)))

            with open(os.path.join(dirpath, fn), 'rb') as f:
                md5calc = md5()
                while True:
                    datachunk = f.read(8192)
                    if not datachunk:
                        break
                    md5calc.update(datachunk)
                md5_checksum_str = md5calc.hexdigest()

            acf = AssetCollectionFile(fn, rp, md5_checksum=md5_checksum_str, tags={'Executable':None} if os.path.splitext(fn)[1] == '.exe' else None)
            ac.add_asset(acf)

    missing_files = ac.save(return_missing_files=True)

    # If COMPS responds that we're missing some files, then try creating it again,
    # uploading only the files that COMPS doesn't already have.
    if missing_files:
        print('Uploading missing files: [' + ','.join([ str(u) for u in missing_files]) + ']')

        ac2 = AssetCollection()
        ac2.set_tags(tags)

        for acf in ac.assets:
            if acf.md5_checksum in missing_files:
                rp = acf.relative_path
                fn = acf.file_name
                acf2 = AssetCollectionFile(fn, rp, tags=acf.tags)
                ac2.add_asset(acf2, os.path.join(path_to_ac, rp, fn))
            else:
                ac2.add_asset(acf)

        ac2.save()
        ac = ac2

    print('done - created AC ' + str(ac.id))

    return ac.id


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print('\r\nUsage:\r\n\t{0} C:\\path\\to\\create\\asset\\collection\\from'.format(sys.argv[0]))
        exit()

    Client.login(compshost)

    create_asset_collection(sys.argv[1])