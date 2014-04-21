#This program will run and graph the simulation, comparing serial and parallel
import sys
import ConfigParser
from subprocess import call
#If not the independent variable, these values will be used
MOLECULES = 2112
SIZE = 40
STEPS = 10000
DENSITY = .033 #the density of water
#function to change the time steps in demoConfiguration.txt
def changeConfigFile(timeSteps, molecules, size):
	f = open('resources/demo.config','r')
	lines = f.readlines()
	lines[1] = str(size) + '\n'
	lines[2] = str(size) + '\n'
	lines[3] = str(size) + '\n'
	lines[9] = str(timeSteps) + '\n'
	lines[11] = str(molecules) + '\n'
	f.close()
	f = open('resources/demo.config','w')
	for item in lines:
		f.write(item)
	f.close()
#function to get the run-time from the last run
def getRunTime():
	config = ConfigParser.RawConfigParser()
	config.read('run.results')
	line = config.get('Results','Run-Time')
	line = line.split()
	return line[0]
#function to print out the help information
def printHelp():
	print "***HELP***"
	print "Usage: python graphing/graphSim.py xaxis interval"
	print "xaxis is what parameter should be on the xaxis: 'steps', 'size', or 'density"
	print "interval is the xvalues separated by spaces"
	print "density is in number of molecules/volume"
	print "Example: python graphing/graphSim density .01 .02 .03 .04\n"
#Begin program
input = sys.argv
#check for invalud inputs
if len(input) < 3:
	printHelp()
	exit(0)
input.remove("graphing/graphSim.py")
xaxis = input.pop(0)
if xaxis != "steps" and xaxis != "size" and xaxis != "density":
	print "Error: invalid format"
	printHelp()
	exit(0)
#capitalize for ease of printing
xaxis = xaxis.capitalize()
#run simulation and prepare output
f = open('graphing/simResults.txt', 'w')
f.write(xaxis + '\tSerial\tParallel\n')
for i in input:
	if xaxis == "Steps":
		changeConfigFile(i, MOLECULES, SIZE)
	if xaxis == "Size":
		changeConfigFile(STEPS, int(i)*int(i)*int(i)*DENSITY, i)
	if xaxis == "Density":
		changeConfigFile(STEPS, float(i)*SIZE*SIZE*SIZE, SIZE) 
	f.write(i+'\t')
	call(["./bin/metrosim", "resources/demo.config", "-s"])
	f.write(getRunTime()+'\t')
	call(["./bin/metrosim", "resources/demo.config", "-p"])
	f.write(getRunTime()+'\n')
f.close()
#set up gnuplot
f = open('graphing/gnuplotScript','w')
f.write("set terminal png\nset output 'graph.png'\nset key left top\n")
f.write("set title 'Serial/Parallel Time Comparison'\n")
f.write("set ylabel 'Time (s)'\nset xlabel '" + xaxis)
f.write("'\nplot 'graphing/simResults.txt' using 1:2 w lp title 'Serial', 'graphing/simResults.txt' using 1:3 w lp title 'Parallel'")
f.close()
#call gnuplot
call(["gnuplot", "graphing/gnuplotScript"])