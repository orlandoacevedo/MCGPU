#include "configurationTest.h"

void testConfigScan()
{
    cout << "Configuration completed" << endl;
    string configPath = "bin/configurationTest.txt"; 
    string oplsPath = "path/to/opls/file";
    string zMatrixPath = "path/to/zMatrix/file";
    string stateInputPath = "path/to/state/input";
    string stateOutputPath = "path/to/State/output";
    string pdbOutputPath = "path/to/pdb/output";
    Config_Scan cs (configPath);
    cs.readInConfig();
    //test configuration path
    assert(configPath.compare(cs.getConfigPath()) == 0);
    //test opls path
    assert(oplsPath.compare(cs.getOplsusaparPath()) == 0);
    //test zmatrix path
    assert(zMatrixPath.compare(cs.getZmatrixPath()) == 0);
    //test state input path
    assert(stateInputPath.compare(cs.getStatePath()) == 0);
    //test state ouput path
    assert(stateOutputPath.compare(cs.getStateOutputPath()) == 0);
    //test pdb output
    assert(pdbOutputPath.compare(cs.getPdbOutputPath()) == 0);

    //test box dimensions
    Environment* enviro = cs.getEnviro();
    double x = 10.342;
    double y = 1234.45;
    double z = 100.3;
    double temperature = 293;
    double maxTranslation = .5;
    int numberOfSteps = 10000;
    int numberOfMolecules = 500;
    
    assert(enviro->x == x);
    assert(enviro->y == y);
    assert(enviro->z == z);

    assert(enviro->maxTranslation == maxTranslation);
    assert(enviro->numOfMolecules == numberOfMolecules);
    assert(enviro->temperature == temperature);
    assert(cs.getSteps() == numberOfSteps); 

    cout << "Configuration completed Complete\n" << endl;
}

