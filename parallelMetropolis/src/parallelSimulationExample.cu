#include <cuda.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "../../Utilities/src/Opls_Scan.h"
#include "../../Utilities/src/Config_Scan.h"
#include "../../Utilities/src/metroUtil.h"
#include "../../Utilities/src/Zmatrix_Scan.h"
#include "../../Utilities/src/State_Scan.h"
#include "../../Utilities/src/geometricUtil.h"
#include "GPUSimBox.cuh"
#include "parallelSim.cuh"

/**
Will run simulations of any system consisting of a single molecule type.  Can run from either z matrix
file or from a state file.

Date: 10/31/2012
Author(s): Xiao Zhang(David),Matthew Hardwick, Riley Spahn, Seth Wooten, Alexander Luchs, and Orlando Acevedo
*/

// boltzman constant
//const double kBoltz = 1.987206504191549E-003;

//const double maxRotation = 10.0; // degrees

using namespace std;

stringstream ss;

/**
  ./bin/parallelSim flag path/to/config/file
    flags:
        -z run from z matrix file spcecified in the configuration file
        -s run from state input file specified in the configuration file
*/

int main(int argc, char ** argv)
{
    char statename[255];
    writeToLog("",START);
    clock_t startTime, endTime;
    startTime = clock();    
       
    if(argc != 2)
    {
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
    // Copy the config_path argument.
   	string configPath(argv[1]);
   	Config_Scan configScan(configPath);
    configScan.readInConfig();

    unsigned int seed = configScan.getrandomseed();

    if (seed==0)
    {
        seed = (unsigned int)time(NULL);
    }
    srand(seed);//init rand sequnce for box positions

    //Environment for the simulation
    GPUSimBox box(configScan);
    SimBox *hostbox=box.getSimBox();
    ss << "Using seed number in Simulation:"<<seed<<endl;
    cout << ss.str();writeToLog(ss);
    
    srand(seed);//init rand sequnce for simulation
     
    Environment *enviro=hostbox->getEnviro();
    unsigned long simulationSteps = configScan.getSteps();
    
    hostbox->WriteStateFile("InitState.state");
    
    ss << "Starting first step in Simulation" <<endl;
    writeToLog(ss);

    //Molecule *molecules;
    ParallelSim sim(&box,100);   
        
    ss << "\nBeginning simulation with: " << endl;
    ss << "\tmolecules "<< hostbox->getEnviro()->numOfMolecules << endl;
    ss << "\tatoms: "<< hostbox->getEnviro()->numOfAtoms << endl;
    ss << "\tsteps: "<< simulationSteps << endl;
    cout << ss.str() <<endl;
    writeToLog(ss);
    //runLinear(molecules, &enviro, simulationSteps, configScan.getStateOutputPath(),
    configScan.getPdbOutputPath();
    
    double initEnergy=sim.calcEnergyWrapper(sim.getGPUSimBox());
   	ss << "Step Number: "<< 0 <<  endl;
	ss << "Current Energy: " << initEnergy << endl;
	cout << ss.str();
   	writeToLog(ss);
    
    for(int i=0;i<simulationSteps;i+=100)
    {
        sim.runParallel(100);        
		ss << "Step Number: "<< i+100 <<  endl;
		ss << "Current Energy: " << sim.getcurrentEnergy() << endl;
		cout << ss.str();
		ss << "Accepted: "<< sim.getaccepted() << endl <<"Rejected: "<< sim.getrejected() << endl;
    	writeToLog(ss);

        sprintf(statename,"%dState.state",i+100);
        hostbox->WriteStateFile(statename);
    }

    ss << "Steps Complete"<<endl;        
    ss << "Final Energy: " << sim.getcurrentEnergy() << endl;
    ss << "Accepted Moves: " << sim.getaccepted() << endl;
    ss << "Rejected Moves: " << sim.getrejected() << endl;
    ss << "Acceptance Rate: " << (int) ((float) sim.getaccepted()/ (float) simulationSteps*100) << "%" << endl;
	cout << ss.str();writeToLog(ss);
         
    endTime = clock();
    double diffTime = difftime(endTime, startTime) / CLOCKS_PER_SEC;
    ss << "\nSimulation Complete \nRun Time: " << diffTime << endl;
    cout << ss.str() <<endl;
    writeToLog(ss);
	
}
