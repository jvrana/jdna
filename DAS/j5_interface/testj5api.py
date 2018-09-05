'''
Project: jdna
File: testj5api
Author: Justin
Date: 2/9/17

Description: 

'''
from j5api import *
import json

credentials = None
with open('j5_credentials.json', 'rU') as handle:
    credentials = json.load(handle)

api = J5API(**credentials)
print((api.run_assembly("NewDesign")))