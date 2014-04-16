#allows system calls and command line arguments
import sys
from subprocess import call
#interval for what trials to run
interval = sys.argv
interval.remove("graphing/graphSim.py")
#function to change the time steps in demoConfiguration.txt
def changeTimeSteps(timeSteps, molecules, size):
	f = open('resources/demoConfiguration.txt','r')
	lines = f.readlines()
	lines[1] = str(size) + '\n'
	lines[2] = str(size) + '\n'
	lines[3] = str(size) + '\n'
	lines[9] = str(timeSteps) + '\n'
	lines[11] = str(molecules) + '\n'
	f.close()
	f = open('resources/demoConfiguration.txt','w')
	for item in lines:
		f.write(item)
	f.close()
def getTime(fileName):
	inFile = open(fileName,'r')
	line = inFile.readlines()
	return line[0]
#create charts
f = open('graphing/simResults.txt', 'w')
f.write('# Steps\tSerial\tParallel\n')
for i in interval:
	s = str(i)
	changeTimeSteps(s, 400, 40.0)
	f.write(s+'\t')
	call(["./bin/metrosim", "resources/demoConfiguration.txt", "-s"])
	f.write(getTime('graphing/lastRun.txt')+'\t')
	call(["./bin/metrosim", "resources/demoConfiguration.txt", "-p"])
	f.write(getTime('graphing/lastRun.txt')+'\n')
f.close()
call(["rm", "graphing/lastRun.txt"])
call(["gnuplot", "graphing/gnuplotScript"])
