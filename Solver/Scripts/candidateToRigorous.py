import re
import sys
import time

counter = 1

outtxt = sys.argv[2]

aa = sys.argv[3]
timeStep = re.findall( ".*?{", aa )
timeStep = re.sub( ".*?:\s", "", timeStep[0] )
timeStep = re.sub( "\s{", "", timeStep )

bb = re.sub( ".*?{", "", aa )
a = re.sub("}", "", bb)
e = re.split(",", a)

outtxt += "*{"
for i in e:
    outtxt += "[" + i + "," + i +"],"
outtxt = outtxt[:-1]
outtxt += "}"
print outtxt
