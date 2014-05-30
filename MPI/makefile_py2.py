import fileinput
import shutil
import subprocess
import os
import uuid

#build_ref=str(uuid.uuid4)
#os.mkdir(build_ref)
#os.chdir(build_ref)

c_names=[]

def build():
    global all
    ### duplicate function file
    shutil.copy2( f_name, f_build )

    #apply changes
    for line in fileinput.input( f_build, inplace=True ):
        print line.replace(C_placeholder, C),

    for line in fileinput.input( f_build, inplace=True ):
        print line.replace(Q_placeholder, Q),

    for line in fileinput.input( f_build, inplace=True ):
        print line.replace(R_placeholder, R),

    for line in fileinput.input( f_build, inplace=True ):
        print line.replace(r_placeholder, r),

    for line in fileinput.input( f_build, inplace=True ):
        print line.replace(img_d_placeholder, img_d),

    for line in fileinput.input( f_build, inplace=True ):
        print line.replace(D_placeholder, D),
    #ADD MORE


    ### duplicate mpi file
    shutil.copy2( mpi_name, mpi_build )

    #apply changes
    for line in fileinput.input( mpi_build, inplace=True ):
        print line.replace(mpi_f_placeholder, f_build),

    #compile new mpic++
    p = subprocess.Popen('mpic++ '+mpi_build+" -o "+mpi_c+" -std=c++11", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print line

    #see exceptions
    retval = p.wait()
    print retval

##### GENERAL VALUES

f_name="doTheMath.cpp" #name of the virgin function script
mpi_name="MPIbasejump.cpp" #name of the virgin mpi launcher file


C_placeholder="@C@"
Q_placeholder="@Q@"
R_placeholder="@R@"
r_placeholder="@r@"
img_d_placeholder="@img_d@"
D_placeholder="@D@"


mpi_f_placeholder="@function@"

##### Build 1 :: ...description...
f_build="f1.cpp"
mpi_build="mpi1.cpp"

mpi_c="mpi1.o"
c_names.append(mpi_c)

R="0.01"
r="30"
C="0.01"
Q="30"
img_d="1.0"
D="40"

build()


##### Build 2 ...


##### WRITE OUT BASH LAUNCHER
NODES=4


launcher = open("launcher.sh" , "w+")
launcher.close()
launcher = open("launcher.sh" , "a")
launcher.write("""
#!/bin/sh \n
#PBS -V \n
#PBS -N mag_pendulum_job \n
#PBS -l nodes="""+str(NODES)+""" \n
#PBS -l alberto.anzellotti=10:00:00 \n
#PBS -m bea \n
#PBS -M alberto.anzellotti@science.unitn.it \n
cd $PBS_O_WORKDIR \n
""")
for step in c_names:
	launcher.write("/usr/local/bin/mpirun -np "+str(NODES)+" /home/alberto.anzellotti/"+step+"\n")
	launcher.write("python3 collector.py \n")

launcher.close()
