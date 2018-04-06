import os
import subprocess
import pandas

cwd = os.getcwd()
print(cwd)

scriptsdir = os.path.dirname(__file__)

process = subprocess.Popen(
    ['python2', '{0}/fit_all_ec50_data.py'.format(scriptsdir)],
    stdout=subprocess.PIPE
)

out, err = process.communicate()

print(out)
print(err)