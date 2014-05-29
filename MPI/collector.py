import glob
import os
import uuid
import time
#os.chdir("/mydir")
nums=[]
files=[]
for file in glob.glob("*.pp"):
    nums.append(int(file[0:-3]))
nums.sort()
for num in nums:
	files.append(str(num)+".pp")
print (files)
pname=str(time.time())+"__"+str(uuid.uuid4())+".ppm"
picture = open(pname , 'w+')
picture.close()
picture = open(pname, 'a')
picture.write("P3 \n")
size=str(len(files)-1)
picture.write(size+" "+size+" 1 \n")
for file in files:
	part = open(file, "r")
	picture.write(part.read())
	part.close()
	os.remove(file) #remove when debugging

picture.close()
