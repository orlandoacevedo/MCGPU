#include "metroLinearUtil.h"

int getXFromIndex(int idx){
    int c = -2 * idx;
    int discriminant = 1 - 4 * c;
    int qv = (-1 + sqrtf(discriminant)) / 2;
    return qv + 1;
}

int getYFromIndex(int x, int idx){
    return idx - (x * x - x) / 2;
}

double makePeriodic(double x, double box){
    
    while(x < -0.5 * box){
        x += box;
    }

    while(x > 0.5 * box){
        x -= box;
    }

    return x;

}

double wrapBox(double x, double box){

    while(x > box){
        x -= box;
    }
    while(x < 0){
        x += box;
    }

    return x;
}

void keepMoleculeInBox(Molecule *molecule, Environment *enviro){		
		for (int j = 0; j < molecule->numOfAtoms; j++){
		//X axis
			wrapBox(molecule->atoms[j].x, enviro->x);
		//Y axis
			wrapBox(molecule->atoms[j].y, enviro->y);
		//Z axis
			wrapBox(molecule->atoms[j].z, enviro->z);
		}
}

double calc_lj(Atom atom1, Atom atom2, Environment enviro){
    //store LJ constants locally
    double sigma = calcBlending(atom1.sigma, atom2.sigma);
    double epsilon = calcBlending(atom1.epsilon, atom2.epsilon);
    
    //calculate difference in coordinates
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;

    //calculate distance between atoms
    deltaX = makePeriodic(deltaX, enviro.x);
    deltaY = makePeriodic(deltaY, enviro.y);
    deltaZ = makePeriodic(deltaZ, enviro.z);
    const double r2 = (deltaX * deltaX) +
                      (deltaY * deltaY) + 
                      (deltaZ * deltaZ);

	// Original code
    if (r2 == 0.0){
        return 0.0;
    }
    else{
    	//calculate terms
    	const double sig2OverR2 = pow(sigma, 2) / r2;
   		const double sig6OverR6 = pow(sig2OverR2, 3);
    	const double sig12OverR12 = pow(sig6OverR6, 2);
    	const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
        return energy;
    }

/*
    //Hard code Rcutoff of 9 Ang (81 = Rcutoff^2) as a test
    if (r2 == 0.0){
        return 0.0;
    }
    else if (r2 < 81.0){
    	//calculate terms
    	const double sig2OverR2 = pow(sigma, 2) / r2;
   		const double sig6OverR6 = pow(sig2OverR2, 3);
    	const double sig12OverR12 = pow(sig6OverR6, 2);
    	const double energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
        return energy;
    }
    else{
        return 0.0;
    }
*/
}

void assignAtomPositions(double *dev_doublesX, double *dev_doublesY, double *dev_doublesZ, Molecule *molec, Environment *enviro){
    //Translates each Molecule a random X,Y,and Z direction
	 //By translating every atom in that molecule by that translation

    //for each Molecule...
	 for(int i=0; i<enviro->numOfMolecules; i++){
	     for(int a=0; a<molec[i].numOfAtoms;a++){
		      Atom myAtom  =  molec[i].atoms[a];
		      myAtom.x =  dev_doublesX[i] * enviro->x + myAtom.x;
				myAtom.y =  dev_doublesY[i] * enviro->y + myAtom.y;
				myAtom.z =  dev_doublesZ[i] * enviro->z + myAtom.z;
		  }
		   keepMoleculeInBox(&molec[i],enviro);
    }
}

void generatePoints(Molecule *molecules, Environment *enviro){

    srand((unsigned int) time(NULL));
	 //for each Molecule assign a new XYZ
    for (int i = 0; i < enviro->numOfMolecules; i++){
        double baseX = ( (double) rand() / RAND_MAX) * enviro->x;
        double baseY = ( (double) rand() / RAND_MAX) * enviro->y;
        double baseZ = ( (double) rand() / RAND_MAX) * enviro->z;
        for (int j = 0; j < molecules[i].numOfAtoms; j++){
            molecules[i].atoms[j].x += baseX;
            molecules[i].atoms[j].y += baseY;
            molecules[i].atoms[j].z += baseZ;
        }

        keepMoleculeInBox(&(molecules[i]), enviro);
    }
}

void generatefccBox(Molecule *molecules, Environment *enviro){
	
	double cells, dcells, cellL, halfcellL;
	
	//Determine the number of unit cells in each coordinate direction
	dcells = pow(0.25 * (double) enviro->numOfMolecules, 1.0/3.0);
	cells = (int)(dcells + 0.5);
		
	//Check if numOfMolecules is a non-fcc number of molecules
	//and increase the number of cells if necessary
	while((4 * cells * cells * cells) < enviro->numOfMolecules)
			cells++;
			
	//Determine length of unit cell
	cellL = enviro->x/ (double) cells;
	halfcellL = 0.5 * cellL;
	
	//Construct the unit cell
	for (int j = 0; j < molecules[0].numOfAtoms; j++){
	molecules[0].atoms[j].x += 0.0;
    molecules[0].atoms[j].y += 0.0;
    molecules[0].atoms[j].z += 0.0;
	}
	
	for (int j = 0; j < molecules[1].numOfAtoms; j++){
	molecules[1].atoms[j].x += halfcellL;
    molecules[1].atoms[j].y += halfcellL;
    molecules[1].atoms[j].z += 0.0;
    }
    
    for (int j = 0; j < molecules[2].numOfAtoms; j++){	
	molecules[2].atoms[j].x += 0.0;
    molecules[2].atoms[j].y += halfcellL;
    molecules[2].atoms[j].z += halfcellL;
    }
    
    for (int j = 0; j < molecules[3].numOfAtoms; j++){
    molecules[3].atoms[j].x += halfcellL;
    molecules[3].atoms[j].y += 0.0;
    molecules[3].atoms[j].z += halfcellL;
    }
    
	//Init all other molecules to initial coordinates
	//Build the lattice from the unit cell by repeatedly translating
	//the four vectors of the unit cell through a distance cellL in
	//the x, y, and z directions
	for(int i = 4; i < enviro->numOfMolecules; i++){
		for (int j = 0; j < molecules[i].numOfAtoms; j++){
			molecules[i].atoms[j].x += 0.0;
    		molecules[i].atoms[j].y += 0.0;
   	 		molecules[i].atoms[j].z += 0.0;
   	 	}		
	}
	
	int offset = 0;
	for(int z = 1; z <= cells; z++)
		for(int y = 1; y <= cells; y++)
			for(int x = 1; x <= cells; x++){
				for(int a = 0; a < 4; a++){
					int i = a + offset;
					if(i < enviro->numOfMolecules){								
						for (int j = 0; j < molecules[i].numOfAtoms; j++){
							molecules[i].atoms[j].x = molecules[a].atoms[j].x + cellL * (x-1);
							molecules[i].atoms[j].y = molecules[a].atoms[j].y + cellL * (y-1);
							molecules[i].atoms[j].z = molecules[a].atoms[j].z + cellL * (z-1);
						}
					}
				}
				offset += 4;
			}
	
	//Shift center of box to the origin
	for(int i = 0; i < enviro->numOfMolecules; i++){
		for (int j = 0; j < molecules[i].numOfAtoms; j++){
			molecules[i].atoms[j].x -= halfcellL;
			molecules[i].atoms[j].y -= halfcellL;
			molecules[i].atoms[j].z -= halfcellL;
		}
	}
}

double calcEnergyWrapper(Molecule *molecules, Environment *enviro){
    //Copy atoms out of molecules
    Atom *atoms = (Atom *) malloc(sizeof(Atom) * enviro->numOfAtoms);
    int atomIndex = 0;
    for(int i = 0; i < enviro->numOfMolecules; i++){
        Molecule currentMolecule = molecules[i];
        for(int j = 0; j < currentMolecule.numOfAtoms; j++){
            atoms[atomIndex] = currentMolecule.atoms[j];
            atomIndex++;
        }
    }

    //pass to original wrapper
    double totalEnergy = calcEnergyWrapper(atoms, enviro, molecules);
    free(atoms);

    return totalEnergy;
}

double calcEnergyWrapper(Atom *atoms, Environment *enviro, Molecule *molecules){
    //setup storage
    double totalEnergy = 0.0;
    double *energySum_device;

    //determine number of energy calculations
    int N =(int) ( pow( (float) enviro->numOfAtoms,2)-enviro->numOfAtoms)/2;	 
    size_t energySumSize = N * sizeof(double);
	double* energySum = (double*) malloc(energySumSize);

    //calulate all energies
    calcEnergy(atoms, enviro, energySum);
    
    for(int i = 0; i < N; i++){
        //apply fudge factor
        if (molecules != NULL){ 
            int atomXid = getXFromIndex(i);
            int atomYid = getYFromIndex(atomXid, i);
            energySum[i] = energySum[i] * getFValue(&(atoms[atomXid]), &(atoms[atomYid]), molecules, enviro); 
        }

        totalEnergy += energySum[i];
    }

    free(energySum);
    return totalEnergy;
}

void calcEnergy(Atom *atoms, Environment *enviro, double *energySum){
    double lj_energy,charge_energy, fValue, nonbonded_energy;

    //determine number of calculations
    int N =(int) ( pow( (float) enviro->numOfAtoms,2)-enviro->numOfAtoms)/2;

    //for each calculation
    for(int idx=0; idx<N; idx++){
        //calculate the x and y positions in the Atom array
        int xAtom_pos, yAtom_pos;
        xAtom_pos = getXFromIndex(idx);
        yAtom_pos = getYFromIndex(xAtom_pos, idx);

        Atom xAtom, yAtom;
        xAtom = atoms[xAtom_pos];
        yAtom = atoms[yAtom_pos];

        //determine the lennard-jones and charge sum between the two atoms
        if (xAtom.sigma < 0 || xAtom.epsilon < 0 || yAtom.sigma < 0 || yAtom.epsilon < 0){
            energySum[idx] = 0.0;
        }
        else{
            //lj_energy = calc_lj(xAtom,yAtom,*enviro);
            //charge_energy = calcCharge(xAtom, yAtom, enviro);
            nonbonded_energy = calcNonBondEnergy(xAtom, yAtom, enviro);

            //store the sum in array
            //energySum[idx] = (lj_energy + charge_energy);
            energySum[idx] = (nonbonded_energy);
        }
	 }
}

double calcCharge(Atom atom1, Atom atom2, Environment *enviro){
    // conversion factor below for units in kcal/mol
    const double e = 332.06;
 
    //calculate difference in coordinates
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;

    //calculate distance between atoms
    deltaX = makePeriodic(deltaX, enviro->x);
    deltaY = makePeriodic(deltaY, enviro->y);
    deltaZ = makePeriodic(deltaZ, enviro->z);

    const double r2 = (deltaX * deltaX) +
                      (deltaY * deltaY) + 
                      (deltaZ * deltaZ);

 
	// Original Code   
    if (r2 == 0.0){
        return 0.0;
    }
    else{
        const double r = sqrt(r2);
        return (atom1.charge * atom2.charge * e) / r;
    }
    
/*
    //Hard code Rcutoff of 9 Ang (81 = Rcutoff^2) as a test
    if (r2 == 0.0){
        return 0.0;
    }
    else if (r2 < 81.0){
        const double r = sqrt(r2);
        return (atom1.charge * atom2.charge * e) / r;
    }
    else{
        return 0.0;
    }
*/
}

double calcNonBondEnergy(Atom atom1, Atom atom2, Environment *enviro){
    //store LJ constants locally and define terms in kcal/mol
    const double e = 332.06;
    double sigma = calcBlending(atom1.sigma, atom2.sigma);
    double epsilon = calcBlending(atom1.epsilon, atom2.epsilon);
    
    //square cutoff value for easy comparison
    double cutoffSQ = enviro->cutoff * enviro->cutoff;
        
    //calculate difference in coordinates
    double deltaX = atom1.x - atom2.x;
    double deltaY = atom1.y - atom2.y;
    double deltaZ = atom1.z - atom2.z;
  
  	//calculate distance between atoms
    deltaX = makePeriodic(deltaX, enviro->x);
    deltaY = makePeriodic(deltaY, enviro->y);
    deltaZ = makePeriodic(deltaZ, enviro->z);
    
    const double r2 = (deltaX * deltaX) +
                	  (deltaY * deltaY) + 
                	  (deltaZ * deltaZ);

    //check if atoms overlap
    if (r2 == 0.0){
        return 0.0;
    }
    else if (r2 < cutoffSQ){
    	//calculate LJ energies
    	const double sig2OverR2 = (sigma * sigma) / r2;
        const double sig6OverR6 = sig2OverR2 * sig2OverR2 * sig2OverR2;
    	const double sig12OverR12 = sig6OverR6 * sig6OverR6;
    	const double lj_energy = 4.0 * epsilon * (sig12OverR12 - sig6OverR6);
    	
    	//calculate Coulombic energies
    	const double r = sqrt(r2);
    	double charge_energy = (atom1.charge * atom2.charge * e) / r;
    	
    	//return total nonbonded energy
    	return (lj_energy + charge_energy);
    	
	}
	else{
        return 0.0;
    }
}

double calcEnergyWrapper_NLC(Molecule *molecules, Environment *enviro){
    //setup storage
    double totalEnergy = 0.0;
    //Copy atoms out of molecules
    Atom *atoms = (Atom *) malloc(sizeof(Atom) * enviro->numOfAtoms);
    int atomIndex = 0;
    for(int i = 0; i < enviro->numOfMolecules; i++){
        Molecule currentMolecule = molecules[i];
        for(int j = 0; j < currentMolecule.numOfAtoms; j++){
            atoms[atomIndex] = currentMolecule.atoms[j];
            atomIndex++;
        }
    }

    totalEnergy = calcEnergy_NLC(atoms, enviro, molecules);
    free(atoms);

    return totalEnergy;
}

double calcEnergy_NLC(Atom *atoms, Environment *enviro, Molecule *molecules){
	// Variables for linked-cell neighbor list	
	int lc[3];            /* Number of cells in the x|y|z direction */
	double rc[3];         /* Length of a cell in the x|y|z direction */
	int head[NCLMAX];     /* Headers for the linked cell lists */
	int mc[3];			  /* Vector cell */
	int lscl[NMAX];       /* Linked cell lists */
	int mc1[3];			  /* Neighbor cells */
	double rshift[3];	  /* Shift coordinates for periodicity */
	const double Region[3] = {enviro->x, enviro->y, enviro->z};  /* MD box lengths */
	int c1;				  /* Used for scalar cell index */
	double dr[3];		  /* Pair vector dr = atom[i]-atom[j] */
	double rrCut = enviro->cutoff * enviro->cutoff;	/* Cutoff squared */
	double rr;			  /* Distance between atoms */
	double nonbonded_energy;	/* Holds current nonbonded energy */
	double totalEnergy = 0.0;	/* Total nonbonded energy x fudge factor */
			
	// Compute the # of cells for linked cell lists
	for (int k=0; k<3; k++) {
		lc[k] = Region[k] / enviro->cutoff; 
		rc[k] = Region[k] / lc[k];
	}
		
  /* Make a linked-cell list, lscl--------------------------------------------*/
	int lcyz = lc[1]*lc[2];
	int lcxyz = lc[0]*lcyz;
		
	// Reset the headers, head
	for (int c = 0; c < lcxyz; c++) 
		head[c] = EMPTY;

	// Scan atoms to construct headers, head, & linked lists, lscl
	for (int i = 0; i < enviro->numOfAtoms; i++) {
		mc[0] = atoms[i].x / rc[0]; 
		mc[1] = atoms[i].y / rc[1];
		mc[2] = atoms[i].z / rc[2];
		
		// Translate the vector cell index, mc, to a scalar cell index
		int c = mc[0]*lcyz + mc[1]*lc[2] + mc[2];

		// Link to the previous occupant (or EMPTY if you're the 1st)
		lscl[i] = head[c];

		// The last one goes to the header
		head[c] = i;
	} /* Endfor molecule i */

  /* Calculate pair interaction-----------------------------------------------*/
		
	// Scan inner cells
	for (mc[0] = 0; mc[0] < lc[0]; (mc[0])++)
	for (mc[1] = 0; mc[1] < lc[1]; (mc[1])++)
	for (mc[2] = 0; mc[2] < lc[2]; (mc[2])++) {

		// Calculate a scalar cell index
		int c = mc[0]*lcyz + mc[1]*lc[2] + mc[2];
		// Skip this cell if empty
		if (head[c] == EMPTY) continue;

		// Scan the neighbor cells (including itself) of cell c
		for (mc1[0] = mc[0]-1; mc1[0] <= mc[0]+1; (mc1[0])++)
		for (mc1[1] = mc[1]-1; mc1[1] <= mc[1]+1; (mc1[1])++)
		for (mc1[2] = mc[2]-1; mc1[2] <=mc [2]+1; (mc1[2])++) {
			// Periodic boundary condition by shifting coordinates
			for (int a = 0; a < 3; a++) {
				if (mc1[a] < 0)
					rshift[a] = -Region[a];
				else if (mc1[a] >= lc[a])
					rshift[a] = Region[a];
				else
					rshift[a] = 0.0;
			}
			// Calculate the scalar cell index of the neighbor cell
			c1 = ((mc1[0] + lc[0]) % lc[0]) * lcyz
			    +((mc1[1] + lc[1]) % lc[1]) * lc[2]
			    +((mc1[2] + lc[2]) % lc[2]);
			// Skip this neighbor cell if empty
			if (head[c1] == EMPTY) continue;

			// Scan atom i in cell c
			int i = head[c];
			while (i != EMPTY) {

				// Scan atom in cell c1
				int j = head[c1];
				while (j != EMPTY) {

					// Avoid double counting of pairs
					if (i < j) {
						// Pair vector dr = atom[i]-atom[j]
						rr = 0.0;
						dr[0] = atoms[i].x - (atoms[j].x + rshift[0]);
						dr[1] = atoms[i].y - (atoms[j].y + rshift[1]);
						dr[2] = atoms[i].z - (atoms[j].z + rshift[2]);
						rr = (dr[0] * dr[0]) + (dr[1] * dr[1]) + (dr[2] * dr[2]);
						}
						
						// Calculate energy if rij < Cutoff
						if (rr < rrCut) {
							Atom xAtom, yAtom;
        					xAtom = atoms[i];
        					yAtom = atoms[j];
        						if (xAtom.sigma < 0 || xAtom.epsilon < 0 || yAtom.sigma < 0 || yAtom.epsilon < 0){
            						nonbonded_energy = 0.0;
        						}
        						else{
    								nonbonded_energy = calcNonBondEnergy(xAtom, yAtom, enviro);
    								nonbonded_energy = nonbonded_energy * getFValue(&xAtom, &yAtom, molecules, enviro);
    							}
    						totalEnergy += nonbonded_energy;
						} /* Endif i<j */

						j = lscl[j];
					} /* Endwhile j not empty */

					i = lscl[i];
				} /* Endwhile i not empty */
			} /* Endfor neighbor cells, c1 */
		} /* Endfor central cell, c */
		return totalEnergy;
}

double calcBlending(double d1, double d2){
    return sqrt(d1 * d2);
}

Molecule* getMoleculeFromAtomID(Atom *a1, Molecule *molecules, Environment *enviro){
    int atomId = a1->id;
    int currentIndex = enviro->numOfMolecules - 1;
    int molecId = molecules[currentIndex].id;
    while(atomId < molecId && currentIndex > 0){
        currentIndex -= 1;
        molecId = molecules[currentIndex].id;
    }

    return &(molecules[currentIndex]);
}

double getFValue(Atom *atom1, Atom *atom2, Molecule *molecules, Environment *enviro){
    Molecule *m1 = getMoleculeFromAtomID(atom1, molecules, enviro);
    Molecule *m2 = getMoleculeFromAtomID(atom2, molecules, enviro);
    
    if(m1->id != m2->id)
        return 1.0;
    else{
        int hops = hopGE3(atom1->id, atom2->id, m1);
        if (hops == 3)
            return 0.5;
        else if (hops > 3)
            return 1.0;
        else
            return 0.0;
    }
}

int hopGE3(int atom1, int atom2, Molecule *molecule){
    for(int x=0; x< molecule->numOfHops; x++){
        Hop *myHop = &(molecule->hops[x]);
        //compare atoms to each hop struct in molecule
		if((myHop->atom1==atom1 && myHop->atom2==atom2) || (myHop->atom1==atom2 && myHop->atom2==atom1)){
            return myHop->hop;
        }
	}
	 return 0;
}

double Energy_LRC(Molecule *molec, Environment *enviro){
	double Ecut = 0.0;		/* Holds LJ long-range cutoff energy correction */
	double NMOL = enviro->numOfMolecules;	/* Number of molecules in simulation */
	// Need to update to accept more than 1 molecule from z-matrix
	int a = 0;
	int NATM = molec[a].numOfAtoms;	/* Number of atoms in molecule */
	// Volume below will have to accept new values in NPT simulations
	double Vnew = enviro->x * enviro->y * enviro->z;	/* Volume of box in Ang^3 */
	
	// Setup arrays to hold all sigma and epsilon atom values
	double Sigma[NATM], Epsilon[NATM];
	for(int i = 0; i < NATM; i++){
		if (molec[a].atoms[i].sigma < 0 || molec[a].atoms[i].epsilon < 0){
			Sigma[i] = 0.0;
			Epsilon[i] = 0.0;
		}
		else{
			Sigma[i] = molec[a].atoms[i].sigma;
			Epsilon[i] = molec[a].atoms[i].epsilon;
		}
	}
		
	// (4 * epsilon * sigma^6 or 12)^0.5
	double A6[NATM], A12[NATM];
	double sig2, sig6, sig12;
	for(int i = 0; i < NATM; i++){
		sig2 = pow(Sigma[i], 2);
        sig6 = pow(sig2, 3);
    	sig12 = pow(sig6, 2);
		A6[i] = sqrt(4 * Epsilon[i] * sig6);
		A12[i] = sqrt(4 * Epsilon[i] * sig12);
	}	
	
	double RC3 = 1.00 / pow(enviro->cutoff, 3);		/* 1 / cutoff^3 */
	double RC9 = pow(RC3, 3);		/* 1 / cutoff^9 */
	// Loop over all atoms in a pair
	for(int i = 0; i < NATM; i++)
		for(int j = 0; j < NATM; j++)
			Ecut += (2*PI*NMOL*NMOL/(3.0*Vnew)) * (A12[i]*A12[j]*RC9/3.0 - A6[i]*A6[j]*RC3);
	
	return Ecut;
}
