#!/usr/bin/python                                                                                                          
import os
import sys

JOBDIR="/ebio/ag-neher/share/users/rneher/RecAda/CoalescentRecombination/"
command = '/ebio/ag-neher/share/programs/bin/python '+ JOBDIR+"src/"+sys.argv[1]
for arg in sys.argv[2:]:
    command+=" "+arg

print command
os.system(command)


