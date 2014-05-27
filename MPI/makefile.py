import fileinput
import shutil
import subprocess



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
f_build="f1.cpp"
mpi_build="mpi1.cpp"

mpi_c="mpi1.o"

R="0.01"
r="30"
C="0.01"
Q="30"
img_d="1.0"
D="40"

build()

##### Build 2 ...
