#include "parallelTest.cuh"

void testKeepMoleculeInBox(){
    double X = 10.0;
    double Y = 10.0;
    double Z = 10.0;
    double temp = 100.0;
    double maxTrans = .5;
    int numOfAtoms = 3;

    Environment enviro = createEnvironment(X, Y, Z, maxTrans, temp, numOfAtoms);

    //test molecule compeletely outside of box.
    Atom a1 = createAtom(1, 11, 10.3, 5);
    Atom a2 = createAtom(2, 12.4, 1.2, 5);
    Atom a3 = createAtom(3, 8.1, 2, 1.5);
    Atom *atoms = (Atom *)malloc(sizeof(Atom) * numOfAtoms);
    atoms[0] = a1;
    atoms[1] = a2;
    atoms[2] = a3;

    Molecule molec;
    molec.atoms = atoms;
    molec.numOfAtoms = numOfAtoms;
    
    double expectedA1[] = {2.9, 9.1, 5};
    double expectedA2[] = {4.3, 0, 5};
    double expectedA3[] = {0, .8, 1.5};

    double **answers = (double **)malloc(sizeof(double) * 9);
    answers[0] = expectedA1;
    answers[1] = expectedA2;
    answers[2] = expectedA3;

    keepMoleculeInBox(&molec, &enviro);
    
    double precision = .0001;
    for(int i = 0; i < numOfAtoms; i++){
        double expectedX = answers[i][0];
        double actualX = molec.atoms[i].x;
        
        double expectedY = answers[i][1];
        double actualY = molec.atoms[i].y;
        
        double expectedZ = answers[i][2];
        double actualZ = molec.atoms[i].z;

        assert(fabs(expectedX - actualX) < precision);
        assert(fabs(expectedY - actualY) < precision);
        assert(fabs(expectedZ - actualZ) < precision);

    }

    cout << "keepMoleculeInBox() passed tests" << endl;
    printAtoms(molec.atoms, 3);
}

void setupGetIndexTest(){
    int numberOfBlocks = 3;
    int threadsPerBlock = 2;
    int totalTests = numberOfBlocks * threadsPerBlock;

    int *xValues;
    int *yValues;
    int *yValues_device;
    int *xValues_device;
    
    size_t xSize = totalTests * sizeof(int);
    
    yValues = (int *) malloc(xSize);
    xValues = (int *)malloc(xSize);
    cudaMalloc((void **) &yValues_device, xSize);
    cudaMalloc((void **) &xValues_device, xSize);
    
    testGetXKernel <<<numberOfBlocks, threadsPerBlock>>>(xValues_device, totalTests);

    cudaMemcpy(xValues, xValues_device, xSize, cudaMemcpyDeviceToHost);

    assert(xValues[0] == 1);
    assert(xValues[1] == 2);
    assert(xValues[2] == 2);
    assert(xValues[3] == 3);
    assert(xValues[4] == 3);
    assert(xValues[5] == 3);

    printf("getXFromIndex Correct\n");

    //test getYFromIndex
    testGetYKernel <<<numberOfBlocks, threadsPerBlock>>> (xValues_device,
            yValues_device, totalTests);

    cudaMemcpy(yValues, yValues_device, xSize, cudaMemcpyDeviceToHost);

    assert(yValues[0] == 0);
    assert(yValues[1] == 0);
    assert(yValues[2] == 1);
    assert(yValues[3] == 0);
    assert(yValues[4] == 1);
    assert(yValues[5] == 2);

    printf("getYFromIndex Correct.\n");

    cudaFree(xValues_device);
    cudaFree(yValues_device);
    free(yValues);
    free(xValues);
}

bool compareDouble(double a, double b, double limit){
    if((a - b) / b < limit)
        return true;
    else
        return false;
}

void setupMakePeriodic(){
    srand(time(NULL));
    int numberOfTests = 128;
    double *box;;
    
    double *inputs_host;
    double *inputs_device;
    double *outputs_host;
    double *dev_box;
    size_t inputSize = sizeof(double) * numberOfTests;

    box = (double *) malloc(sizeof(double));
    *box = 10.0;
    inputs_host = (double *) malloc(inputSize);
    outputs_host = (double *) malloc(inputSize);
    cudaMalloc((void **) &inputs_device, inputSize);
    cudaMalloc((void **) &dev_box, sizeof(double));
    
    //generate random numbers
    for(int i = 0; i < numberOfTests; i++){
        inputs_host[i] = ((double) (rand() % 100));
    }

    //copy data to device
    cudaMemcpy(inputs_device, inputs_host, inputSize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_box, box, sizeof(double), cudaMemcpyHostToDevice);
    
    int threadsPerBlock = numberOfTests / 2;
    int blocks = numberOfTests / threadsPerBlock +
        (numberOfTests % threadsPerBlock == 0 ? 0 : 1);

    testMakePeriodicKernel <<< blocks, threadsPerBlock >>> (inputs_device,
            dev_box, numberOfTests);

    cudaMemcpy(outputs_host, inputs_device, inputSize, cudaMemcpyDeviceToHost);

    //check that values are the same as known correct function
    for(int i = 0; i < numberOfTests; i++){
        double test_output = make_periodic(inputs_host[i], *box);
        assert(outputs_host[i] == test_output);
    }

    printf("makePeriodic passed Tests\n");

    free(inputs_host);
    free(outputs_host);
    cudaFree(inputs_device);
}


void testWrapBox(){
 srand(time(NULL));
    int numberOfTests = 128;
    double box;
    
    double *testDoubles;
    size_t inputSize = sizeof(double) * numberOfTests;

    box = 10.f;
    testDoubles = (double *) malloc(inputSize);
    
    //generate random numbers
    for(int i = 0; i < numberOfTests; i++){
        testDoubles[i] = (double) rand() / (double) RAND_MAX;
    }

     //check that values are the same as known correct function
    for(int i = 0; i < numberOfTests; i++){
        double test_output = wrap_into_box(testDoubles[i], box);
        assert(wrapBox(testDoubles[i], box) == test_output);
    }

    printf("wrapBox passed Tests\n");

    free(testDoubles);
}

void setupCalc_lj(){
    double kryptonSigma = 3.624;
    double kryptonEpsilon = 0.317;
    int numberOfAtoms = 2;

    Atom *atoms = new Atom[numberOfAtoms];
    double *energy = (double *) malloc(sizeof(double));
    *energy = 1000.f;
    Atom *atoms_device;
    Environment *enviro_device;
    double *energy_device;

    cudaMalloc((void **) &atoms_device, sizeof(Atom) * numberOfAtoms);
    cudaMalloc((void **) &enviro_device, sizeof(Environment));
    cudaMalloc((void **) &energy_device, sizeof(double));

    Environment stableEnviro = createEnvironment(10, 10, 10, .5,
            298.15, numberOfAtoms);

    Environment *enviro = &stableEnviro;
    generatePoints(atoms, enviro);
    atoms[0].sigma = kryptonSigma;
    atoms[0].epsilon = kryptonEpsilon; 
    atoms[1].sigma = kryptonSigma;
    atoms[1].epsilon = kryptonEpsilon;

    cudaMemcpy(atoms_device, atoms, sizeof(Atom) * numberOfAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy(enviro_device, enviro, sizeof(Environment), cudaMemcpyHostToDevice);

    testCalcLJ<<<1,1>>>(atoms_device, enviro_device, energy_device);

    cudaMemcpy(energy, energy_device, sizeof(double), cudaMemcpyDeviceToHost);

    double baseEnergy = calculate_energy(atoms, enviro);
    assert((int)(*energy * pow(10.f, 6.f)) == (int)( baseEnergy * pow(10.f,6.f))); 
    printf("\nparallelEnergy = %2.10f\nlinearEnergy = %2.10f\n", *energy, baseEnergy); 
    printf("Calc_lj is correct\n");
    free(atoms);
    free(energy);
    cudaFree(atoms_device);
    cudaFree(enviro_device);
    cudaFree(energy_device);
}


void testGeneratePoints(){
    //init atoms, environment
    int numberOfAtoms = 1000;
    Atom *atoms = (Atom *) malloc(numberOfAtoms * sizeof(Atom));

    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, 0, 0, 0);
    }
    Environment enviro = createEnvironment(10.0, 20.0, 35.0, 1.0, 298.15, numberOfAtoms);

    generatePoints(atoms, &enviro);

    //assert that all atoms positions are in range of the box
    for (int i = 0; i < numberOfAtoms; i++){
        double dim_x = atoms[i].x;
        double dim_y = atoms[i].y;
        double dim_z = atoms[i].z;

        assert(dim_x >= 0 && dim_x <= (enviro.x) &&
               dim_y >= 0 && dim_y <= (enviro.y) &&
               dim_z >= 0 && dim_z <= (enviro.z));
    }

    printf("testGeneratePoints (atoms) successful.\n");

    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, i % 2, i % 3, i % 4);
    }

    int numberOfMolecules = 250;

    Molecule *molecules = (Molecule *)malloc(sizeof(Molecule)*numberOfMolecules);
    for (int i = 0; i < numberOfMolecules; i++){
        Bond *blankBonds = NULL;
        Angle *blankAngles = NULL;
        Dihedral *blankDihedrals = NULL;
        int atomCount = numberOfAtoms / numberOfMolecules;
        Atom *molAtoms = (Atom *) malloc(sizeof(Atom)*atomCount);
        for (int j = 0; j < atomCount; j++){
            molAtoms[j] = atoms[i*atomCount + j];
        }
        molecules[i] = createMolecule(-1, molAtoms, blankAngles, blankBonds, blankDihedrals, atomCount, 0, 0, 0);
    }
    
    generatePoints(molecules, &enviro);

    for (int i = 0; i < numberOfMolecules; i++){
        for (int j = 0; j < molecules[i].numOfAtoms; j++){
            double dim_x = molecules[i].atoms[j].x;
            double dim_y = molecules[i].atoms[j].y;
            double dim_z = molecules[i].atoms[j].z;
            
             assert(dim_x >= 0 && dim_x <= (enviro.x) &&
               dim_y >= 0 && dim_y <= (enviro.y) &&
               dim_z >= 0 && dim_z <= (enviro.z));
        }
    }

    printf("testGeneratePoints (molecules) successful.\n");

    free(atoms);
    free(molecules);
}

void testCalcEnergy(){
    // the sigma value of nitrogen
    double nSigma = 3.250;
    // the epsilon value of nitrogen
    double nEpsilon = 0.17;

    // charge of nitrogen
    double nCharge = -0.850;

    struct timeval le_tvBegin, le_tvEnd, pl_tvBegin, pl_tvEnd;

    //Generate enviorment and atoms
    int numberOfAtoms = 1000;
    Environment stableEnviro = createEnvironment(5.0, 10.0, 15.0, 1.0, 298.15, numberOfAtoms);

    Environment *enviro = &stableEnviro;

    Atom *atoms = new Atom[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; i++){
        atoms[i] = createAtom(i, 0.0, 0.0, 0.0, nSigma, nEpsilon);
        atoms[i].charge = nCharge;
    }

    generatePoints(atoms, enviro);
    
    //calculate energy linearly
    gettimeofday(&le_tvBegin,NULL); //start clock for execution time

    double te_linear = calculate_energy(atoms, enviro);

    gettimeofday(&le_tvEnd,NULL); //stop clock for execution time
    long le_runTime = timevaldiff(&le_tvBegin,&le_tvEnd); //get difference in time in milli seconds

    //calculate energy in parallel
    gettimeofday(&pl_tvBegin,NULL); //start clock for execution time

    double te_parallel =  calcEnergyWrapper(atoms, enviro);	 

    gettimeofday(&pl_tvEnd,NULL); //start clock for execution time
    long pl_runTime = timevaldiff(&pl_tvBegin,&pl_tvEnd); //get difference in time in milli seconds


    //Print out Results
    printf("Number of elements: %d\n", numberOfAtoms);
    printf("Linear Total Energy:   %f \n", te_linear);
    printf("In %lu ms\n", le_runTime);
    printf("Parallel Total Energy: %f \n", te_parallel);
    printf("In %lu ms\n", pl_runTime);
    assert(compareDouble(te_linear, te_parallel, .05));
    printf("testCalcEnergy successful.\n");

    
}

void testCalcEnergyWithMolecules(){
    // the sigma value of nitrogen
    double nSigma = 3.250;
    // the epsilon value of nitrogen
    double nEpsilon = 0.17;
    // charge of nitrogen
    double nCharge = -0.850;

    // the sigma value of krypton used in the LJ simulation
    double kryptonSigma = 3.624;
    // the epsilon value of krypton used in the LJ simulation
    double kryptonEpsilon = 0.317;
    // charge of krypton
    double kryptonCharge = 0.0;

    struct timeval le_tvBegin, le_tvEnd, pl_tvBegin, pl_tvEnd;

    //Generate enviorment and atoms
    int numberOfAtoms = 500;
    Environment stableEnviro = createEnvironment(5.0, 10.0, 15.0, 1.0, 298.15, numberOfAtoms);

    Environment *enviro = &stableEnviro;

    Atom *atoms = new Atom[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; i++){
        if ((i % 5) < 3){
            atoms[i] = createAtom(i, 0.0, 0.0, 0.0, kryptonSigma, kryptonEpsilon);
            atoms[i].charge = kryptonCharge;
        }
        else{
            atoms[i] = createAtom(i, 0.0, 0.0, 0.0, nSigma, nEpsilon);
            atoms[i].charge = nCharge;
        }
    }
    enviro->numOfMolecules = 100;

    Molecule *molecules;
    molecules = (Molecule *)malloc(sizeof(Molecule) * numberOfAtoms);

    int atomCount = 0;
    for(int i = 0; i < enviro->numOfMolecules; i++){
        molecules[i].numOfAtoms = 5;
        molecules[i].atoms = (Atom *)malloc(sizeof(Atom) * 5);
        molecules[i].id = atomCount;
        Hop *hops = (Hop *)malloc(sizeof(Hop) * 2);
        hops[0] = createHop(atomCount, atomCount+3, 3);
        hops[1] = createHop(atomCount, atomCount+4, 4);
        molecules[i].hops = hops;
        molecules[i].numOfHops = 2;

        for (int j = 0; j < molecules[i].numOfAtoms; j++){
            molecules[i].atoms[j] = atoms[atomCount];
            atomCount++;
        }
    }

    generatePoints(molecules, enviro);

    for (int i = 0; i < enviro->numOfMolecules; i++){
        for (int j = 0; j < molecules[i].numOfAtoms; j++){
            atoms[i * 5 + j] = molecules[i].atoms[j];
        }
    }

    //calculate energy linearly
    gettimeofday(&le_tvBegin,NULL); //start clock for execution time

    double te_linear = calculate_energy(atoms, enviro, molecules);

    gettimeofday(&le_tvEnd,NULL); //stop clock for execution time
    long le_runTime = timevaldiff(&le_tvBegin,&le_tvEnd); //get difference in time in milli seconds

    //calculate energy in parallel
    gettimeofday(&pl_tvBegin,NULL); //start clock for execution time

    double te_parallel =  calcEnergyWrapper(molecules, enviro);	 

    gettimeofday(&pl_tvEnd,NULL); //start clock for execution time
    long pl_runTime = timevaldiff(&pl_tvBegin,&pl_tvEnd); //get difference in time in milli seconds

    //Print out Results
    printf("Number of elements: %d\n", numberOfAtoms);
    printf("Linear Total Energy:   %f \n", te_linear);
    printf("In %lu ms\n", le_runTime);
    printf("Parallel Total Energy: %f \n", te_parallel);
    printf("In %lu ms\n", pl_runTime);
    assert(compareDouble(te_linear, te_parallel, .05));
    printf("testCalcEnergyWithMolecules successful.\n");
}

void testGetMoleculeFromIDWrapper(){
    int numberOfAtoms = 11;
    int numberOfMolecules = 3;
    
    Atom *atoms;
    DeviceMolecule *molecules;
    Environment enviro;
    int *answers;

    Atom *atoms_device;
    DeviceMolecule *molecules_device;
    int *answers_device;

    enviro.numOfAtoms = numberOfAtoms;
    enviro.numOfMolecules = numberOfMolecules;

    atoms = (Atom *)malloc(sizeof(Atom) * numberOfAtoms);
    molecules = (DeviceMolecule *)malloc(sizeof(DeviceMolecule) *numberOfMolecules);
    answers = (int *)malloc(sizeof(int) * numberOfAtoms);

    cudaMalloc((void **) &atoms_device, sizeof(Atom) * numberOfAtoms);
    cudaMalloc((void **) &molecules_device, sizeof(DeviceMolecule) * numberOfMolecules);
    cudaMalloc((void **) &answers_device, sizeof(int) * numberOfAtoms);

    enviro.numOfAtoms = numberOfAtoms;
    enviro.numOfMolecules = numberOfMolecules;
    
    for(int i = 0; i < numberOfAtoms; i++){
        atoms[i].id = i;
    }
    molecules[0].id = 0;
    molecules[1].id = 2;
    molecules[2].id = 6;

    cudaMemcpy(atoms_device, atoms, sizeof(Atom) * numberOfAtoms, cudaMemcpyHostToDevice);
    cudaMemcpy(molecules_device, molecules, sizeof(DeviceMolecule) * numberOfMolecules, cudaMemcpyHostToDevice);

    int numberOfBlocks = 1;
    int threadsPerBlock = 128;
    testGetMoleculeFromID<<<numberOfBlocks,threadsPerBlock>>>(atoms_device,
            molecules_device, enviro, numberOfAtoms, answers_device);
   
    cudaMemcpy(answers, answers_device, sizeof(int) * numberOfAtoms, cudaMemcpyDeviceToHost);

    assert(answers[0] == 0);
    assert(answers[1] == 0);
    assert(answers[2] == 2);
    assert(answers[3] == 2);
    assert(answers[4] == 2);
    assert(answers[5] == 2);
    assert(answers[6] == 6);
    assert(answers[7] == 6);
    assert(answers[8] == 6);
    assert(answers[9] == 6);
    assert(answers[10] == 6);
   
    printf("getMoleculeFromID passed tests\n");

    free(atoms);
    free(molecules);
    free(answers);

    cudaFree(atoms_device);
    cudaFree(molecules_device);
    cudaFree(answers_device);
}


void testCalcBlendingWrapper(){
    double *d1, *d2, *d1_device, *d2_device, *answers, *answers_device;
    int numberOfTests = 5;
    size_t doubleSize = sizeof(double) * numberOfTests;

    d1 = (double *)malloc(doubleSize);
    d2 = (double *)malloc(doubleSize);
    answers = (double *)malloc(doubleSize);
    
    cudaMalloc((void **) &d1_device, doubleSize);
    cudaMalloc((void **) &d2_device, doubleSize);
    cudaMalloc((void **) &answers_device, doubleSize);

    d1[0] = 0.f;
    d2[0] = 0.f;

    d1[1] = 4.5;
    d2[1] = 2.32;
    
    d1[2] = 52.34;
    d2[2] = 5.f;

    d1[3] = 1.f;
    d2[3] = 7.f;

    d1[4] = 34.56;
    d2[4] = 12.7;
    
    cudaMemcpy(d1_device, d1, doubleSize, cudaMemcpyHostToDevice);
    cudaMemcpy(d2_device, d2, doubleSize, cudaMemcpyHostToDevice);

    int blocks = 1;
    int threadsPerBlock = 64;

    testCalcBlending <<<blocks, threadsPerBlock>>>(d1_device, d2_device, answers_device, numberOfTests);

    cudaMemcpy(answers, answers_device, doubleSize, cudaMemcpyDeviceToHost);

    for(int i = 0 ; i < numberOfTests; i++){
        double expected = sqrt(d1[i] * d2[i]);
        assert(answers[i] / sqrt(d1[i] * d2[i]) < 0.01 || answers[i] == expected);
    }

    printf("calcBlending passed tests.\n");
    
    free(d1);
    free(d2);
    free(answers);
    cudaFree(d1_device);
    cudaFree(d2_device);
    cudaFree(answers_device);
}

void testGetFValueWrapper(){
    Environment *enviro, *dev_enviro;
    DeviceMolecule *molecules, *dev_molecules;
    Atom *mol1_atoms, *mol2_atoms, *atom1List, *atom2List, *dev_atom1List, *dev_atom2List;
    Hop *dev_hops;
    double *fvalues, *dev_fvalues;

    int numberOfTests = 5;

    Environment stable_enviro = createEnvironment(5.0,5.0,5.0,1.0,270.0,5);
    enviro = &stable_enviro;
    mol1_atoms = (Atom *)malloc(sizeof(Atom)*5);
    mol2_atoms = (Atom *)malloc(sizeof(Atom));
    atom1List = (Atom *)malloc(sizeof(Atom)*numberOfTests);
    atom2List = (Atom *)malloc(sizeof(Atom)*numberOfTests);

    fvalues = (double *)malloc(sizeof(double)*4);
    
    for (int i = 0; i < 5; i++){
        mol1_atoms[i] = createAtom(i,1.0,1.0,1.0);
    }
    for (int i = 0; i < 5; i++){
        atom1List[i] = mol1_atoms[0];
        if (i < 4)
            atom2List[i] = mol1_atoms[i+1];
    }

    mol2_atoms[0] = createAtom(5,1.0,1.0,1.0);
    atom2List[4] = mol2_atoms[0];

    Hop* hops = (Hop *)malloc(sizeof(Hop)*2);
    hops[0] = createHop(0,3,3);
    hops[1] = createHop(0,4,4);

    molecules = (DeviceMolecule *)malloc(sizeof(DeviceMolecule)*2);
    molecules[0] = createDeviceMolecule(0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 2);
    molecules[1] = createDeviceMolecule(5, 5, 1, 0, 0, 0, 0, 0, 0, 2, 0);

    cudaMalloc((void **) &dev_enviro, sizeof(Environment));
    cudaMalloc((void **) &dev_molecules, sizeof(DeviceMolecule)*2);
    cudaMalloc((void **) &dev_atom1List, sizeof(Atom)*numberOfTests);
    cudaMalloc((void **) &dev_atom2List, sizeof(Atom)*numberOfTests);
    cudaMalloc((void **) &dev_fvalues, sizeof(double)*numberOfTests);
    cudaMalloc((void **) &dev_hops, sizeof(Hop)*2);

    cudaMemcpy(dev_enviro, enviro, sizeof(Environment), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_molecules, molecules, sizeof(DeviceMolecule)*2, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_atom1List, atom1List, sizeof(Atom)*numberOfTests, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_atom2List, atom2List, sizeof(Atom)*numberOfTests, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_hops, hops, sizeof(Hop)*2, cudaMemcpyHostToDevice);

    int blocks = 1;
    int threadsPerBlock = 64;

    testGetFValue <<<blocks, threadsPerBlock>>>(dev_atom1List, dev_atom2List, dev_molecules, dev_enviro, dev_fvalues, numberOfTests, dev_hops);

    cudaMemcpy(fvalues, dev_fvalues, sizeof(double)*numberOfTests, cudaMemcpyDeviceToHost);

    double *expected = (double *)malloc(sizeof(double)*numberOfTests);
    expected[0] = 0.0;
    expected[1] = 0.0;
    expected[2] = 0.5;
    expected[3] = 1.0;
    expected[4] = 1.0;
    for(int i = 0 ; i < numberOfTests; i++){
        assert(expected[i] == fvalues[i]);
    }

    printf("testGetFValue passed tests.\n");
   
    free(mol1_atoms);
    free(mol2_atoms);
    free(atom1List);
    free(atom2List);
    free(fvalues);
    free(molecules);
    free(hops);
    cudaFree(dev_enviro);
    cudaFree(dev_molecules);
    cudaFree(dev_atom1List);
    cudaFree(dev_atom2List);
    cudaFree(dev_fvalues);
    cudaFree(dev_hops);
}

Atom findMaxRotation(Atom pivot, Atom toRotate, double rotation){
    toRotate.x -= pivot.x;
    toRotate.y -= pivot.y;
    toRotate.z -= pivot.z;

    rotateAboutX(toRotate, rotation);
    rotateAboutY(toRotate, rotation);
    rotateAboutZ(toRotate, rotation);

    toRotate.x += pivot.x;
    toRotate.y += pivot.y;
    toRotate.z += pivot.z;

    return toRotate;
}

void testRotateMolecule(){
    srand(time(NULL));
    
    //Testing on a molecule that is not totaly unlike water
    double bondDistance = 0.9584; // angstroms
    double maxRotation = 10.0; // degrees
    int numOfAtoms = 3;
    int numOfAngles = 1;
    int numOfBonds = 2;
    int numOfDihedrals = 0;

    Atom oxygen = createAtom(1, 0, 0, 0);
    Atom hydrogen1 = createAtom(2, 0, bondDistance, 0);
    Atom hydrogen2 = createAtom(3, bondDistance, 0, 0);

    Atom *atoms = (Atom *)malloc(sizeof(Atom) * 3);
    atoms[0] = oxygen;
    atoms[1] = hydrogen1;
    atoms[2] = hydrogen2;
    
    vector<Atom> atomVector;
    atomVector.push_back(oxygen);
    atomVector.push_back(hydrogen1);
    atomVector.push_back(hydrogen2);
    Bond b1 = createBond(1,2, bondDistance, false);
    Bond b2 = createBond(1,3, bondDistance, false);
    
    Bond *bonds = (Bond *)malloc(sizeof(Bond) * 2);
    bonds[0] = b1;
    bonds[1] = b2;

    Angle a1 = createAngle(2,3,90,false);
    Angle *angles = (Angle *)malloc(sizeof(Angle));
    angles[0] = a1;

    Dihedral *dihedrals = (Dihedral *)malloc(sizeof(Dihedral) * 0);

    Molecule molec;
    molec = createMolecule(1, atoms, angles, bonds, dihedrals,
            numOfAtoms, numOfAngles, numOfBonds, numOfDihedrals);
    
    int testNumber = 10;
    
    printf("Testing rotateMolecule\n");

    for(int i = 0 ; i < testNumber; i++){
        //pick atom to rotate about.  Cycle through all of them
        Atom toRotate = atoms[1];
        
        
        rotateMolecule(molec, toRotate, maxRotation);
        
        //test that rotation is within limit
        Atom newAtom1 = atoms[2];
        Atom origAtom1 = getAtom(atomVector, newAtom1.id);
        
        double angleChange1 = getAngle(newAtom1, toRotate, origAtom1);
        printf("Atom1 angle change = %f\n", angleChange1);

        Atom newAtom2 = atoms[0];
        Atom origAtom2 = getAtom(atomVector, newAtom2.id);
        double angleChange2 = getAngle(newAtom2, toRotate, origAtom2);
        
        printf("Atom2 angle change = %f\n", angleChange2);
        
        Atom maxAtom1 = findMaxRotation(toRotate, newAtom1, maxRotation);
        Atom maxAtom2 = findMaxRotation(toRotate, newAtom2, maxRotation);
        double maxAngle1 = getAngle(maxAtom1, toRotate, origAtom1);
        double maxAngle2 = getAngle(maxAtom2, toRotate, origAtom2);
       
        assert(angleChange1 <= maxAngle1);
        assert(angleChange2 <= maxAngle2);


        //reset atoms
        molec.atoms[0] = oxygen; 
        molec.atoms[1] = hydrogen1;
        molec.atoms[2] = hydrogen2;
    }
    printf("rotateMolecule passed tests.\n");
}

void testCalcChargeWrapper(){
    printf("Testing calcCharge()\n");
    
    int numberOfTests = 10;
    
    // data on the host
    Atom *atoms1_h;
    Atom *atoms2_h;
    Environment *enviro_h;
    double *answers_h;

    // data on the device
    Atom *atoms1_d;
    Atom *atoms2_d;
    Environment *enviro_d;
    double *answers_d;

    // get sizes of data
    size_t atomSize = sizeof(Atom) * numberOfTests;
    size_t enviroSize = sizeof(Environment);
    size_t answerSize = sizeof(double) * numberOfTests;

    // mallocate on host
    atoms1_h = (Atom *)malloc(atomSize);
    atoms2_h = (Atom *)malloc(atomSize);
    enviro_h = (Environment *)malloc(enviroSize);
    answers_h = (double *) malloc(answerSize);

    // mallocate on device
    cudaMalloc((void **) &atoms1_d, atomSize);
    cudaMalloc((void **) &atoms2_d, atomSize);
    cudaMalloc((void **) &enviro_d, enviroSize);
    cudaMalloc((void **) &answers_d, answerSize);

    double xSize = 10;
    double ySize = xSize;
    double zSize = ySize;

    //generate atoms for test
    srand(time(NULL));
    for(int i = 0; i < numberOfTests; i++){
        atoms1_h[i].x = (double) rand() / (double) RAND_MAX * xSize;
        atoms2_h[i].x = (double) rand() / (double) RAND_MAX * xSize;
        
        atoms1_h[i].y = (double) rand() / (double) RAND_MAX * ySize;
        atoms2_h[i].y = (double) rand() / (double) RAND_MAX * ySize;
        
        atoms1_h[i].z = (double) rand() / (double) RAND_MAX * zSize;
        atoms2_h[i].z = (double) rand() / (double) RAND_MAX * zSize;
   
        atoms1_h[i].charge = (double) rand() / (double) RAND_MAX * 2 - 1;
        atoms2_h[i].charge = (double) rand() / (double) RAND_MAX * 2 - 1; 
    }

    enviro_h->x = xSize;
    enviro_h->y = ySize;
    enviro_h->z = zSize;
    enviro_h->numOfAtoms = numberOfTests;

    //transfer data to the device
    cudaMemcpy(atoms1_d, atoms1_h, atomSize, cudaMemcpyHostToDevice);
    cudaMemcpy(atoms2_d, atoms2_h, atomSize, cudaMemcpyHostToDevice);
    cudaMemcpy(enviro_d, enviro_h, enviroSize, cudaMemcpyHostToDevice);

    //call test function
    int numOfBlocks = 1;
    int threadsPerBlock = 64;
    
    testCalcCharge<<<numOfBlocks, threadsPerBlock>>>(atoms1_d, atoms2_d, answers_d, enviro_d);

    //transfer answers from device to host
    cudaMemcpy(answers_h, answers_d, answerSize, cudaMemcpyDeviceToHost);

    //TEST ANSWERS
    for(int i = 0; i < numberOfTests; i++){
        double expected = calc_charge(atoms1_h[i], atoms2_h[i], *enviro_h);
        assert((expected - answers_h[i]) / expected < .01);
    }

    printf("calcCharge passed tests.\n");

    free(atoms1_h);
    free(atoms2_h);
    free(enviro_h);
    free(answers_h);

    cudaFree(atoms1_d);
    cudaFree(atoms2_d);
    cudaFree(enviro_d);
    cudaFree(answers_d);
}

int main(){
    testCopyMolecules();
    testKeepMoleculeInBox();
    testRotateMolecule();
    testCalcChargeWrapper();
    testCalcBlendingWrapper();
    testGetMoleculeFromIDWrapper();
    testWrapBox();
    setupCalc_lj();
    setupGetIndexTest();
    setupMakePeriodic();
    testGeneratePoints();
    testCalcEnergy();
    testCalcEnergyWithMolecules();
    
    return 0;
}

