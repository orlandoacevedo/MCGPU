#ifndef SIMBOX_BUILDER_H
#define SIMBOX_BUILDER_H

#include "SimBox.h"
#include "Utilities/FileUtilities.h"

class SimBoxBuilder {

private:

  /**
   * sb holds the simulation box object that is being created by the builder
   *     instance.
   */
  SimBox* sb;

  /**
   * sbData points to the SBScanner object that will be used to retrieve data
   *     read in from the OPLSAA.sb / OPLSUA.sb file.
   */
  SBScanner* sbData;

  /**
   * Initializes basic environment variables, such as the box's temperature,
   *     cutoff distance, and dimensions.
   *
   * @param environment The box's environment.
   */
  void initEnvironment(Environment* environment);

  /**
   * Adds molecules from the box to the simulation box. Currently only handles
   *     atoms, not bonds or angles.
   *
   * @param molecules A dynamic array containing all of the molecules in the box.
   */
  void addMolecules(Molecule* molecules);

  /**
   * Adds primary indexes, read in from the config file, to every molecule in
   *     the simulation box.
   *
   * @param primaryAtomIndexArray Points to a vector, which holds pointers to
   *     other vectors containing the primary indexes for each type of molecule.
   */
  void addPrimaryIndexes(std::vector< std::vector<int>* >* primaryAtomIndexArray);

  /**
   * Initializes the NLC of the Simulation Box.
   */
  void fillNLC();

public:

  /**
   * Constructor for SimBoxBuilder.
   *
   * @param useNLC True if this Simulation run will use the NLC, false otherwise.
   * @param sbData Points to SBScanner, which retrieves information about bonds
   *     and angles from oplsaa.sb.
   */
  SimBoxBuilder(bool useNLC, SBScanner* sbData);

  /**
   * Driver function for SimBoxBuilder. Constructs and returns a simulation box
   *     based on a box passed in.
   *
   * @param box Points to a box containing every aspect of the Simulation, as
   *     read in, to be converted into the simulation box.
   */
  SimBox* build(Box* box);

};

#endif
