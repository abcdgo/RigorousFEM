#!/usr/bin/python

import re
import sys
import time
import math
from decimal import Decimal

import os

if( len(sys.argv) > 1 ):
    file_in = open(sys.argv[1], "r")
else:
    file_in = open("a.txt", "r")


def getOverlap(a, b):
    return max( Decimal(0), min(a[1], b[1]) - max(a[0], b[0]) )

def allIn(a, b):
    if( a[0] <= b[0] and a[1] >= b[1] ): return True
    else:
        return False

aa = file_in.readline().rstrip(',\n')
timeStep = re.findall( ".*?{", aa )
timeStep = re.sub( ".*?:\s", "", timeStep[0] )
timeStep = re.sub( "\s{", "", timeStep )
timeStep = re.sub('\[', "", timeStep)
timeStep = re.split(",", timeStep)
timeStep = timeStep[0]
timeStep = Decimal(timeStep)


bb = re.sub( ".*?{", "", aa )
a = re.sub("}", "", bb)
b = re.sub('\[', "", a)
c = re.sub('\]', "", b)
start = re.split(",", c)


aa = file_in.readline().rstrip('\n')

file_in.close()

timeStep1 = re.findall( ".*?{", aa )
timeStep1 = re.sub( ".*?:\s", "", timeStep1[0] )
timeStep1 = re.sub( "\s{", "", timeStep1 )
timeStep1 = re.sub('\[', "", timeStep1)
timeStep1 = re.split(",", timeStep1)
timeStep1 = timeStep1[0]
timeStep1 = Decimal(timeStep1)

bb = re.sub( ".*?{", "", aa )
a = re.sub("}", "", bb)
b = re.sub('\[', "", a)
c = re.sub('\]', "", b)
period = re.split(",", c)

print "Start time: " + str(timeStep) + " after period time: " + str(timeStep1)
if( timeStep == timeStep1 ):
    print "ERROR Start time == End time"

print "Start interval \t\t\t\t\t after period interval"
for i in range( len(start) / 2 ):
    over = getOverlap( [ Decimal(start[2*i]), Decimal(start[2*i+1]) ], [ Decimal(period[2*i]), Decimal(period[2*i+1]) ] )
    if( over != Decimal(0) ):
        if( allIn( [ Decimal(start[2*i]), Decimal(start[2*i+1]) ], [ Decimal(period[2*i]), Decimal(period[2*i+1]) ] ) == True ):
            zawiera = "Tak"
        else:
            zawiera = "Nie"
    else:
        zawiera = "Nie"
    print "[" + start[2*i] + " ; " + start[2*i+1] + "]\t[" + period[2*i] + " ; " + period[2*i+1] + "]" + " \t Przeciecie: " + str( over ) + "\tCzy sie w pelni zawiera: " + zawiera


