#include "TestUtil.h"

void createConfigFile(std::string MCGPU, std::string fileName, std::string primaryAtomIndexString, ConfigFileData settings) {
	std::ofstream configFile;
	//std::cout << "CREATING CONFIG FILE: " << MCGPU << settings.working_path << "/" << fileName;
	std::string configFilePath (std::string (MCGPU + settings.working_path + "/" + fileName));

	configFile.open(configFilePath.c_str());
	std::stringstream cfg;
	cfg << ""
		<< "#size of periodic box (x, y, z in angstroms)" << std::endl
		<< settings.sizeX << std::endl
		<< settings.sizeY << std::endl
		<< settings.sizeZ << std::endl
		<< "#temperature in Kelvin\n"
		<< settings.tempKelvin << std::endl
		<< "#max translation\n"
		<< settings.maxTranslation << std::endl
		<< "#number of steps\n"
		<< settings.numSteps << std::endl
		<< "#number of molecules\n"
		<< settings.numMolecules << std::endl
		<< "#path to opla.par file\n"
		<< MCGPU << settings.oplsaa_par_path << std::endl
		<< "#path to z matrix file\n"
		<< MCGPU << settings.z_matrix_file << std::endl
		<< "#path to state input\n"
		<< MCGPU << settings.working_path << std::endl
        << "#path to state output\n"
        << MCGPU << settings.working_path << std::endl
        << "#pdb output path\n"
        << MCGPU << settings.working_path << std::endl
        << "#cutoff distance in angstroms\n"
        << settings.cutoffDistance << std::endl
        << "#max rotation\n"
        << settings.maxRotation << std::endl
        << "#Random Seed Input\n"
        << settings.randomSeed << std::endl
        << "#Primary Atom Index\n"
        << primaryAtomIndexString;
        configFile << cfg.str();
        configFile.close();
}

std::string getMCGPU_path () {
	std::string directory = get_current_dir_name();
	std::string mc ("MCGPU");
	std::size_t found = directory.find(mc);

	if(found != std::string::npos) {
		directory = directory.substr(0, found + 6);
	}
	std::string MCGPU = directory + "/";
	return MCGPU;
}

bool inDebugMode() {
	std::string directory = get_current_dir_name();
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
		ss << "-s --threads 12 ";	//If testing in series, give the flag and specify a the number of threads.
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
	std::ifstream infile(std::string(MCGPU + "bin/" + resultsFile).c_str());
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
