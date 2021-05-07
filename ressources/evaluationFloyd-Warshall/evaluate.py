from subprocess import STDOUT, check_output, TimeoutExpired, CalledProcessError
from os import listdir
from os.path import isfile, join
import subprocess
import sys
import time

dataFolder="data/"
binFolder="bin/"
srcFolder="src/"

# Change number of processors
#np = 3
np =4
datasets=[2,3]


def check_output(resultFile, output) :
    resultat=""
    with open(resultFile) as f:
        for line in f:
            # print("output :", output)
            #print("line :", line)
            resultat+=line.replace('\n', '').replace(' ', '')
        #print("output :", output)
        #print("result :", resultat)
        #return str(resultat==output)
        if resultat.strip().replace('\n', '')==output:
            return "2"
        else:
            return "1"

def runDataset(name, i) :
    output=""
    try:
        timeout=20
        if (i>5) :
            timeout=60
        binFile=binFolder+name
        print("Running ", binFile)
        if not isfile(binFile) :
            return "0"
        output = subprocess.check_output(["mpirun", "-np", str(np), "-mca", "btl","^openib", binFolder+name, dataFolder+"mat_"+str(i)], stderr=STDOUT, timeout=timeout,universal_newlines=True).strip().replace('\n', '').replace(' ', '')
        #print(output)

        return check_output(dataFolder+"/result_"+str(i), output)

    except  CalledProcessError as e :
        #we still have to check the output because
        #some student call returns non 0 value in case of success !

        output=e.output.strip().replace('\n', '')
        #print(output, " return code =", e.returncode)
        if (e.returncode != 0) :
            return "0"
        return check_output(dataFolder+"/result_"+str(i), output)
    except (TimeoutExpired) as e :
        return "0"



def runProject(name) :
    result=""
    for d in datasets :
        tmp = runDataset(name, d)
        result = result + ";" + tmp
        #print("Testing " + name + " " + tmp, file=sys.stderr )
        try:
            subprocess.run(["killall", name])
        except  CalledProcessError as e :
            pass
        time.sleep(2)
        try :
            subprocess.run(["killall", "-9", name])
        except  CalledProcessError as e :
            pass
    return result

def compileAndRunProject(f):
    print("Testing ", f[:-2], " ",file=sys.stderr)
    resultCompilation=compile(f)
    print(f + ";" + str(resultCompilation) + runProject(f))

def runAllProjects() :
    onlyfiles = [f for f in listdir(binFolder) if isfile(join(binFolder, f))]
    for f in onlyfiles:
        print(runProject(f))

def compile(f) :
    name=binFolder + f
    #print(f)
    try:
        output = subprocess.check_output(["mpicc","-std=c99", "-o" , name,  srcFolder+ f+".c" ,"-lm", "-fopenmp"], stderr=STDOUT, universal_newlines=True).strip().replace('\n', '')
        #resultCompilation=(output=="")
        resultCompilation=0
        if (output==""):
            resultCompilation=1
        # print(output)
    except  CalledProcessError as e :
        # print(e)
        resultCompilation=0
    return resultCompilation

def compileAndRunProjects() :
    onlyfiles = [f for f in sorted(listdir(srcFolder), key=lambda s: s.lower()) if isfile(join(srcFolder, f))]
    for f in onlyfiles:
        print("Testing ", f[:-2], " ",file=sys.stderr)
        resultCompilation=compile(f[:-2])
        print(f[:-2] + ";" + str(resultCompilation)  + runProject(f[:-2]))

if len(sys.argv) > 1 :
    compileAndRunProject(sys.argv[1])
else :
    compileAndRunProjects()
