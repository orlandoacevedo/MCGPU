/*!\file
  \brief Structures and functions used to read and write configuration files.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 */
#include "State_Scan.h"

void printState(Environment *enviro, Molecule *molecules, int numOfMolecules, string filename){
    ofstream outFile;
    outFile.open(filename.c_str());
    //print the environment
    outFile << enviro->x << " " << enviro->y << " " << enviro->z << " " << enviro->numOfAtoms
        << " " << enviro->temperature << endl;
    outFile << endl; // blank line
    for(int i = 0; i < numOfMolecules; i++){
        Molecule currentMol = molecules[i];
        outFile << currentMol.id << endl;
        outFile << "= Atoms" << endl; // delimiter
    
        //write atoms
        for(int j = 0; j < currentMol.numOfAtoms; j++){
            Atom currentAtom = currentMol.atoms[j];
            outFile << currentAtom.id << " "
                << currentAtom.x << " " << currentAtom.y << " " << currentAtom.z
                << " " << currentAtom.sigma << " " << currentAtom.epsilon  << " "
                << currentAtom.charge << endl;
        }
        outFile << "= Bonds" << endl; // delimiter
        
        //write bonds
        for(int j = 0; j < currentMol.numOfBonds; j++){
            Bond currentBond = currentMol.bonds[j];
            outFile << currentBond.atom1 << " " << currentBond.atom2 << " "
                << currentBond.distance << " ";
            if(currentBond.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;

        }
        outFile << "= Dihedrals" << endl; // delimiter
        for(int j = 0; j < currentMol.numOfDihedrals; j++){
            Dihedral currentDi = currentMol.dihedrals[j];
            outFile << currentDi.atom1 << " " << currentDi.atom2 << " "
                << currentDi.value << " ";
            if(currentDi.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;
        }

        outFile << "=Hops" << endl;

        for(int j = 0; j < currentMol.numOfHops; j++){
            Hop currentHop = currentMol.hops[j];

            outFile << currentHop.atom1 << " " << currentHop.atom2 << " "
                << currentHop.hop << endl;
        }
        
        
        outFile << "= Angles" << endl; // delimiter

        //print angless
        for(int j = 0; j < currentMol.numOfAngles; j++){
            Angle currentAngle = currentMol.angles[j];

            outFile << currentAngle.atom1 << " " << currentAngle.atom2 << " "
                << currentAngle.value << " ";
            if(currentAngle.variable)
                outFile << "1" << endl;
            else
                outFile << "0" << endl;
        }


        //write a == line
        outFile << "==" << endl;
    }
    outFile.close();
}

Environment getEnvironmentFromLine(string line){
    Environment enviro;
    //tokenize input line
    char *tokens; 
    char *charLine = (char *) malloc(sizeof(char) * line.size());
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    //extract data from line
    int tokenNumber = 0;
    double x, y, z;
    int numOfAtoms;
    while(tokens != NULL){
        switch(tokenNumber){
            case 0:
                enviro.x = atof(tokens);
                break;
            case 1:
                enviro.y = atof(tokens);
                break;
            case 2:
                enviro.z = atof(tokens);
                break;
            case 3:
                enviro.numOfAtoms = atoi(tokens);
                break;
            case 4:
                enviro.temperature = atof(tokens);
                break;
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }

    return enviro;
}

Atom getAtomFromLine(string line){
    Atom returnAtom = createAtom(-1, -1, -1, -1, -1, -1);
    char *tokens;
    char *charLine = (char *)malloc(sizeof(char) * line.size());
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    
    //read in atoms
    int tokenNumber = 0;
    while(tokens != NULL){
        switch(tokenNumber){
            case 0: // id
                returnAtom.id = atoi(tokens);
                break;
            case 1: // x
                returnAtom.x = atof(tokens);
                break;
            case 2: // y
                returnAtom.y = atof(tokens);
                break;
            case 3: // z
                returnAtom.z = atof(tokens);
                break;
            case 4: // sigma
                returnAtom.sigma = atof(tokens);
                break;
            case 5: // epsilon
                returnAtom.epsilon = atof(tokens);
                break;
            case 6: //charge
                returnAtom.charge = atof(tokens);
                break;
    
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }
    return returnAtom;
}

Bond getBondFromLine(string line){
    Bond bond = createBond( -1, -1, -1, false);

    char *tokens;
    char *charLine = (char *)malloc(sizeof(char) * line.size());
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    strcpy(charLine, line.c_str());

    int tokenNumber = 0;
    while(tokens != NULL){
        switch(tokenNumber){
            case 0: // atom 1
                bond.atom1 = atoi(tokens);
                break;
            case 1: // atom 2
                bond.atom2 = atoi(tokens);
                break;
            case 2: // value
                bond.distance = atof(tokens);
                break;
            case 3: // variable
                if(atoi(tokens) == 1)
                    bond.variable = true;
                else
                    bond.variable = false;
                break;
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }
    return bond;
}

Angle getAngleFromLine(string line){
    Angle angle = createAngle(-1, -1, -1, false);

    char *tokens;
    char *charLine = (char *)malloc(sizeof(char) * line.size());
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    strcpy(charLine, line.c_str());
    
    int tokenNumber = 0;
    while(tokens != NULL){
       switch(tokenNumber){
           case 0:
               angle.atom1 = atoi(tokens);
               break;
            case 1:
               angle.atom2 = atoi(tokens);
               break;
            case 2:
               angle.value = atof(tokens);
               break;
            case 3:
               if(atoi(tokens) == 1)
                   angle.variable = true;
               else
                   angle.variable = false;
       }
       tokens = strtok(NULL, " ");
       tokenNumber++;
    }
    return angle;
}

Dihedral getDihedralFromLine(string line){
    Dihedral dihedral = createDihedral(-1, -1, -1, false);

    char *tokens;
    char *charLine = (char *)malloc(sizeof(char) * line.size());
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    int tokenNumber = 0;
    while(tokens != NULL){
        switch(tokenNumber){
            case 0:
                dihedral.atom1 = atoi(tokens);
                break;
            case 1:
                dihedral.atom2 = atoi(tokens);
                break;
            case 2:
                dihedral.value = atof(tokens);
                break;
            case3:
                if(atoi(tokens) == 1)
                    dihedral.variable = true;
                else
                    dihedral.variable = false;
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }
    return dihedral;

}

//returns a hop based on the information on the line
Hop getHopFromLine(string line){
    Hop hop = createHop(-1, -1, -1);

    char *tokens;
    char *charLine = (char *)malloc(sizeof(char) * line.size());
    
    strcpy(charLine, line.c_str());
    tokens = strtok(charLine, " ");
    int tokenNumber = 0;
    while(tokens != NULL){
        switch(tokenNumber){
            case 0:
                hop.atom1 = atoi(tokens);
                break;
            case 1:
                hop.atom2 = atoi(tokens);
                break;
            case 2:
                hop.hop = atoi(tokens);
                break;
        }
        tokens = strtok(NULL, " ");
        tokenNumber++;
    }
    return hop;
}

Environment readInEnvironment(string filename){
    ifstream inFile (filename.c_str());
    string line;
    Environment enviro;

    if(inFile.is_open()){
        getline(inFile, line);
        enviro = getEnvironmentFromLine(line);    
    }
    else
        return enviro;

    return enviro;
        
}


vector<Molecule> readInMolecules(string filename){
    vector<Molecule> molecules;
    ifstream inFile(filename.c_str());
    string line;
    
    if(inFile.is_open()){
        vector<Bond> bonds;
        vector<Angle> angles;
        vector<Atom> atoms;
        vector<Dihedral> dihedrals;
        vector<Hop> hops;

        //the third line starts the first molecule
        getline(inFile, line); // envrionment
        getline(inFile, line); //blank
        
        Molecule currentMolec;
        int section = 0; // 0 = id, 1 = atom, 2 = bond, 3 = dihedral, 4 = hop, 5 = angle
        while(inFile.good()){
            //printf("bonds: %d\nangles: %d\natoms: %d\ndihedrals: %d\n\n",
              //      bonds.size(), angles.size(), atoms.size(), dihedrals.size());
            getline(inFile, line);
            switch(section){
                case 0: // id
                    if(line.compare("=") == 0)
                        section++;
                    else{
                        currentMolec.id = atoi(line.c_str());
                    }
                    break;
                case 1: // atom
                    if(line.compare("=") == 0)
                        section++;
                    else{
                       atoms.push_back(getAtomFromLine(line)); 
                    }
                    break;
                case 2: // bond
                    if(line.compare("=") == 0)
                        section++;
                    else{
                       bonds.push_back(getBondFromLine(line)); 
                    }
                    break;
                case 3: // dihedral
                    if(line.compare("=") == 0)
                        section++;
                    else{
                       dihedrals.push_back(getDihedralFromLine(line)); 
                    }
                    break;
                case 4: // hop
                    if(line.compare("=") == 0)
                        section++;
                    else{
                        hops.push_back(getHopFromLine(line));
                    }
                case 5: // angle
                    if(line.compare("==") == 0){
                        section = 0;
                        
                        //convert all vectors to arrays
                        Bond *bondArray = (Bond *) malloc(sizeof(Bond) * bonds.size());
                        Angle *angleArray = (Angle *) malloc(sizeof(Angle) * angles.size());
                        Atom *atomArray = (Atom *) malloc(sizeof(Atom) * atoms.size());
                        Dihedral *dihedralArray = (Dihedral *) malloc(sizeof(Dihedral) * dihedrals.size());
                        Hop *hopArray = (Hop *) malloc(sizeof(Hop) * hops.size());

                        for(int i = 0; i < bonds.size(); i++)
                            bondArray[i] = bonds[i];
                        for(int i = 0; i < angles.size(); i++)
                            angleArray[i] = angles[i];
                        for(int i = 0; i < atoms.size(); i++)
                            atomArray[i] = atoms[i];
                        for(int i = 0; i < dihedrals.size(); i++)
                            dihedralArray[i] = dihedrals[i];
                        for(int i = 0; i < hops.size(); i++)
                            hopArray[i] = hops[i];
                       
                        //assign arrays to molecule
                        currentMolec.atoms = atomArray;
                        currentMolec.numOfAtoms = atoms.size();
                        
                        currentMolec.bonds = bondArray;
                        currentMolec.numOfBonds = bonds.size();
                        
                        currentMolec.angles = angleArray;
                        currentMolec.numOfAngles = angles.size();

                        currentMolec.dihedrals = dihedralArray;
                        currentMolec.numOfDihedrals = dihedrals.size();

                        currentMolec.hops = hopArray;
                        currentMolec.numOfHops = hops.size();
                        
                        //add molecule to array of returned molecules
                        molecules.push_back(currentMolec); 
                        
                        Molecule newMolec;
                        currentMolec = newMolec;

                        dihedrals.clear();
                        atoms.clear();
                        bonds.clear();
                        angles.clear();
                        hops.clear();

                    }
                    else{
                       angles.push_back(getAngleFromLine(line)); 
                    }
                    break;
            }
        }
    }
    else{
        return molecules;
    }
    
   return molecules;

}


