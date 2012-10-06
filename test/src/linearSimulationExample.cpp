#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <sstream>
#include "../../Utilities/src/Opls_Scan.h"
#include "../../Utilities/src/Config_Scan.h"
#include "../../Utilities/src/metroUtil.h"
#include "../../Utilities/src/Zmatrix_Scan.h"
#include "../../Utilities/src/State_Scan.h"
#include "linearSim.h"
#include "SimBox.h"
   
using namespace std;



stringstream ss;

double randomFloat(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}



/**
  ./bin/linearSim flag path/to/config/file
    flags:
        -z run from z matrix file spcecified in the configuration file
        -s run from state input file specified in the configuration file
*/
int main(int argc, char ** argv){
    writeToLog("",START);
    clock_t startTime, endTime;
    startTime = clock();
	
    if(argc != 2){
	
		// If not, produce an error message.
		
        printf("Error.  Expected usage: bin/linearSim {configPath}\n");
        exit(1);
    }
	
	// Is our config_path argument valid?
   
	if (argv[1] == NULL)
	{
	
		// If not, produce an error message.
	
		ss << "configuration file argument null"<<endl;
        cout <<ss.str()<< endl; writeToLog(ss);
		exit(1);
	}
	
    // Copy the config_path argument.
    
	string configPath(argv[1]);

    // Scan in the configuration file properties.
	
    Config_Scan configScan(configPath);
    configScan.readInConfig();
	
	// Declare a flag string.
	
	string flag;
	
	// Do we have a statefile path configured?
	
	if (!configScan.getStatePath().empty())
	{
		flag = string("-s");
	}
	
	// No statefile path, so we assume z-matrix path.
	
	else 
	{
		flag = string("-z");
	}

    //Environment for the simulation
    SimBox box(configPath.c_str());
    
    Environment *enviro=box.getEnviro();
    unsigned long simulationSteps = configScan.getSteps();
    //Molecule *molecules;
    LinearSim sim(&box,100);
    

    ss << "\nBeginning simulation with: " << endl;
    ss << "\tmolecules "<< enviro->numOfMolecules << endl;
    ss << "\tatoms: "<< enviro->numOfAtoms << endl;
    ss << "\tsteps: "<< simulationSteps << endl;
    cout << ss.str() <<endl;
    writeToLog(ss);
    
    ss << "Assigning Molecule Positions..." << endl;
    writeToLog(ss);
	  box.generatefccBox(box.getMolecules(), box.getEnviro());
   
    ss << "Finished Assigning Molecule Positions" << endl;
    writeToLog(ss);
		
    box.WriteStateFile("initialState");

    ss << "Starting first step in Simulation" <<endl;
    writeToLog(ss);

    char fileName[50];	 
    for(int move = 0; move < simulationSteps; move=move+100){
    	sim.runLinear(100);

    	//Print the state every 100 moves.

    	sprintf(fileName, "%dState.state", move+100);
    	string fileNameStr(fileName);
    	box.WriteStateFile(fileNameStr);

			ss << "Step Number: "<< move+100 <<  endl;
			ss << "Current Energy: " << sim.getcurrentEnergy() << endl;
			cout << ss.str();

    	ss << "Accepted: "<< sim.getaccept() << endl <<"Rejected: "<< sim.getreject() << endl;
    	writeToLog(ss);
    }

    sprintf(fileName, "FinalState.state");
    string fileNameStr(fileName);
    box.WriteStateFile(fileNameStr);
    box.writePDB("pdbFile.pdb");
    
    ss << "Steps Complete"<<endl;        
    ss << "Final Energy: " << sim.getcurrentEnergy() << endl;
    ss << "Accepted Moves: " << sim.getaccept() << endl;
    ss << "Rejected Moves: " <<  sim.getreject() << endl;
    ss << "Acceptance Rate: " << (int) ((float) sim.getaccept()/ (float) simulationSteps) << "%" << endl;

    cout << ss.str();
	  writeToLog(ss);
	  
    endTime = clock();

    double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;
    ss << "\nSimulation Complete \nRun Time: " << diffTime << endl;
    cout << ss.str() <<endl;
    writeToLog(ss);      
    writeToLog("",END);
}
