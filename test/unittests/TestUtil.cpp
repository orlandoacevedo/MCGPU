#include "TestUtil.h"

// Generates a config file with the given attributes
void createConfigFile(std::string MCGPU, std::string fileName,
		std::string primaryAtomIndexString,
		ConfigFileData settings) {
	std::ofstream configFile;
	std::string configFilePath (std::string(MCGPU + settings.working_path + "/" + fileName));

	configFile.open(configFilePath.c_str());
	std::stringstream cfg;
	cfg << ""
		<< "# Size of periodic box (x, y, z in angstroms)" << std::endl
		<< "x=" << settings.sizeX << std::endl
		<< "y=" << settings.sizeY << std::endl
		<< "z=" << settings.sizeZ << std::endl
		<< "# Temperature in Kelvin\n"
		<< "temp=" << settings.tempKelvin << std::endl
		<< "# Max translation\n"
		<< "max-translation=" << settings.maxTranslation << std::endl
		<< "# Number of steps\n"
		<< "steps=" << settings.numSteps << std::endl
		<< "# Number of molecules\n"
		<< "molecules=" << settings.numMolecules << std::endl
		<< "# Path to opla.par file\n"
		<< "opla.par=" << MCGPU << settings.oplsaa_par_path << std::endl
		<< "# Path to z matrix file\n"
		<< "z-matrix=" << MCGPU << settings.z_matrix_file << std::endl
		<< "# Path to state input\n"
		<< "state-input=" << MCGPU << settings.working_path << std::endl
		<< "# Path to state output directory\n"
		<< "state-output=" << MCGPU << settings.working_path << std::endl
		<< "# Path to pdb output directory\n"
		<< "pdb-output=" << MCGPU << settings.working_path << std::endl
		<< "# Cutoff distance in angstroms\n"
		<< "cutoff=" << settings.cutoffDistance << std::endl
		<< "# Max rotation\n"
		<< "max-rotation=" << settings.maxRotation << std::endl
		<< "# Random Seed Input\n"
		<< "random-seed=" << settings.randomSeed << std::endl
		<< "# Primary Atom Index\n"
		<< "primary-atom=" << primaryAtomIndexString;
		configFile << cfg.str();
		configFile.close();
}

std::string getMCGPU_path () {
	char *path = (char*) malloc(4096);
	size_t size = 4096;
	path = getcwd(path, size);
	//std::string directory = get_current_dir_name();
	std::string directory(path);
	std::string mc ("MCGPU");
	std::size_t found = directory.rfind(mc);

	if(found != std::string::npos) {
		directory = directory.substr(0, found + 6);
	}
	std::string MCGPU = directory + "/";
	return MCGPU;
}

bool inDebugMode() {
	char *path = (char*) malloc(4096);
	size_t size = 4096;
	path = getcwd(path, size);
	//std::string directory = get_current_dir_name();
	std::string directory(path);
	//std::string directory = get_current_dir_name();
	std::string mcdebug ("MCGPU/bin/debug");
	std::size_t found = directory.find(mcdebug);

	if(found != std::string::npos) {
		return true;
	}
	return false;
}

std::string buildCommand(std::string MCGPU, std::string configFile, std::string outputName, bool series, bool neighborlist, bool errorExpected, std::string working_path) {
	//Setting up standard build command.
	std::stringstream ss;
	std::string metrosimCommand;
	if(inDebugMode())
		metrosimCommand = "bin/debug/metrosim ";
	else
		metrosimCommand = "bin/metrosim ";

	ss << MCGPU << metrosimCommand << MCGPU << working_path << "/" << configFile << " ";

	if(series) {
		ss << "-s ";				//If testing in series, give the flag.
	} else {
		ss << "-p ";				//If testing in parallel, give the corresponding flag.
	}

	if(neighborlist) {
		ss << "-l 10 ";				//Add the neighborlist flag if applicable.
	}

	if(errorExpected) {
		ss << "-i 10000 > " << MCGPU << "bin/" << outputName << " 2>&1 ";	//If we expect an error, pipe cerr to a textfile where it can be read.
	} else {
		ss << "--name " << outputName << " -i 10000 ";							//If we do not expect an error, simply give the name for the results file.
	}
	std::string output = ss.str();
  // std::cout << "RUNNING: " << output << std::endl;
	return output;
}

double getEnergyResult(std::string MCGPU, std::string resultsFile) {
	std::string dir;
	if(inDebugMode())
		dir = "bin/debug/";
	else
		dir = "bin/";

	std::ifstream infile(std::string(MCGPU + dir + resultsFile).c_str());
	std::size_t found;

	for(std::string line; getline(infile, line);) {
		std::string str2 ("Final-Energy");
		std::string result;
		found = line.find(str2);
		if(found != std::string::npos) {
			result = line.substr(15);
			return strtod(result.c_str(), NULL);
		}
	}

	return -1;
}

std::string getErrorResult(std::string MCGPU, std::string errorFile) {
	std::ifstream infile(std::string(MCGPU + "bin/" + errorFile).c_str());
	std::size_t found;

    for(std::string line; getline(infile, line);) {
        std::string str2 ("Error");
        found = line.find(str2);
        if (found != std::string::npos) {
            return line.substr(7,13);
        }
    }
	return "ERROR: COULD NOT PARSE ERROR FILE!";
}
