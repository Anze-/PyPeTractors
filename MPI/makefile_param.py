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
        print(line.replace(C_placeholder, C), end='')

    for line in fileinput.input( f_build, inplace=True ):
        print(line.replace(Q_placeholder, Q), end='')

    for line in fileinput.input( f_build, inplace=True ):
        print(line.replace(R_placeholder, R), end='')

    for line in fileinput.input( f_build, inplace=True ):
        print(line.replace(r_placeholder, r), end='')

    for line in fileinput.input( f_build, inplace=True ):
        print(line.replace(img_d_placeholder, img_d), end='')

    for line in fileinput.input( f_build, inplace=True ):
        print(line.replace(D_placeholder, D), end='')
    #ADD MORE


    ### duplicate mpi file
    shutil.copy2( mpi_name, mpi_build )

    #apply changes
    for line in fileinput.input( mpi_build, inplace=True ):
        print(line.replace(mpi_f_placeholder, f_build), end='')

    #compile new mpic++
    p = subprocess.Popen('mpic++ '+mpi_build+" -o "+mpi_c+" -std=c++11", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print(line)

    #see exceptions
    retval = p.wait()
    print (retval)

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


R="0.03"
r="20"
C="0.01"
Q="3"
img_d="1.0"
D="60"


while float(R) < 0.055:
	Q="3"
	while float(Q) < 30:
		f_build="f"+Q+R+".cpp"
		mpi_build="mpi"+Q+R+".cpp"

		mpi_c="mpi"+Q+R+".o"
		c_names.append(mpi_c)
		build()
		a = float(Q) + 3
		Q = str(a)

	Q="30"

	while float(Q) < 120:
		f_build="f"+Q+R".cpp"
		mpi_build="mpi"+Q+R+".cpp"

		mpi_c="mpi"+Q+R+".o"
		c_names.append(mpi_c)
		build()
		a = float(Q) + 30
		Q = str(a)
	a = float(R) + 0.01
	R = str(a)

##### WRITE OUT BASH LAUNCHER
NODES=4


launcher = open("launcher.sh" , "w+")
launcher.close()
launcher = open("launcher.sh" , "a")
launcher.write("""
#!/bin/sh \n
""")
for step in c_names:
	launcher.write("mpirun -np 8 ./"+step+"\n")
	launcher.write("python3 collector.py \n")

launcher.close()
