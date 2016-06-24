#include "FileUtilities.h"

#include <exception>
#include <stdexcept>
//#include <sstream>
#include "Parsing.h"
#include "StructLibrary.h"
#include "Metropolis/Box.h"
#include "Metropolis/SimulationArgs.h"

void SBScanner::processBond(string line) {
	string atom1, atom2;
	if (line.substr(1, 1) == " ") {
		atom1 = line.substr(0, 1);
	} else {
		atom1 = line.substr(0, 2);
	}

	if (line.substr(4, 1) == " ") {
		atom2 = line.substr(3, 1);
	} else {
		atom2 = line.substr(3, 2);
	}

	double forceK, bondDist;
	string rest_of_line = line.substr(6, line.size()-6);
	stringstream ss(rest_of_line);
	ss >> forceK >> bondDist;
	bondDataMap[atom1][atom2].kBond = forceK;
	bondDataMap[atom2][atom1].kBond = forceK;
	bondDataMap[atom1][atom2].eqBondDist = bondDist;
	bondDataMap[atom2][atom1].eqBondDist = bondDist;
}

void SBScanner::processAngle(string line) {
  string end1 = line.substr(0, 2);
	string midAtom = line.substr(3, 2);
	string end2 = line.substr(6, 2);

  if (line.substr(1, 1) == " ") {
    end1 = end1.substr(0, 1);
  }

  if (line.substr(4, 1) == " ") {
    midAtom = midAtom.substr(0, 1);
  }

  if (line.substr(7, 1) == " ") {
    end2 = end2.substr(0, 1);
  }

	string rest_of_line = line.substr(8, line.size()-8);
	double angleK, angle;
	stringstream ss(rest_of_line);
	ss >> angleK >> angle;

  angleDataMap[end1][midAtom][end2].kAngle = angleK;
	angleDataMap[end2][midAtom][end1].kAngle = angleK;
	angleDataMap[end1][midAtom][end2].eqAngle = angle;
	angleDataMap[end2][midAtom][end1].eqAngle = angle;
}

void SBScanner::processLine(string line) {
	if(line.size() > 5 && line.at(5) == '-') {
		processAngle(line);
  } else {
		processBond(line);
  }
}

SBScanner::SBScanner() {

}

SBScanner::~SBScanner() {
	bondDataMap.clear();
	angleDataMap.clear();
}

bool SBScanner::readInSB(string filename) {
	fileName = filename;

	if (filename.empty()) {
		std::cerr << "Error: readInOpls(): empty filename given" << std::endl;
		return false;
	}

	int numOfLines = 0;

  ifstream sbScanner(filename.c_str());
  if (!sbScanner.is_open()) {
    std::cerr << "Error: readInOpls(): could not open file (" << filename
              << ")" << std::endl;
    return false;
  }
  string line;

  while (sbScanner.good()) {
    numOfLines++;
    getline(sbScanner, line);

    //check if it is a commented line,
    //or if it is a title lines
    try {
      if(line.at(0) != '*' && line.at(0) != '\"' && line.at(0) != '#') {
        processLine(line);
      }
    } catch (std::out_of_range& e) {
      // Eat the exception and continue
    }
  }
  sbScanner.close();
	return true;
}

Real SBScanner::getKBond(string atom1, string atom2) {
	if(bondDataMap.count(atom1) == 0) return -1;
	if(bondDataMap[atom1].count(atom2) == 0) return -1;
	return bondDataMap[atom1][atom2].kBond;
}

Real SBScanner::getEqBondDist(string atom1, string atom2) {
	if(bondDataMap.count(atom1) == 0) return -1;
	if(bondDataMap[atom1].count(atom2) == 0) return -1;
	return bondDataMap[atom1][atom2].eqBondDist;
}

Real SBScanner::getKAngle(string endpoint1, string middleAtom, string endpoint2) {
	if(angleDataMap.count(endpoint1) == 0) return -1;
	if(angleDataMap[endpoint1].count(middleAtom) == 0) return -1;
	if(angleDataMap[endpoint1][middleAtom].count(endpoint2) == 0) return -1;
	return angleDataMap[endpoint1][middleAtom][endpoint2].kAngle;
}

Real SBScanner::getEqAngle(string endpoint1, string middleAtom, string endpoint2) {
	if(angleDataMap.count(endpoint1) == 0) return -1;
	if(angleDataMap[endpoint1].count(middleAtom) == 0) return -1;
	if(angleDataMap[endpoint1][middleAtom].count(endpoint2) == 0) return 0;
	return angleDataMap[endpoint1][middleAtom][endpoint2].eqAngle;
}
