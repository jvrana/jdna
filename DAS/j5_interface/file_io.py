'''
Project: jdna
File: file_io
Author: Justin
Date: 2/9/17

Description: 

'''
import zipfile
import os

def compress_directory(directory):

    with zipfile.ZipFile(directory + '.zip', 'w') as zf:
        for root, _, filenames in os.walk(os.path.basename(directory)):
            for name in filenames:
                name = os.path.join(root, name)
                name = os.path.normpath(name)
                zf.write(name, name)