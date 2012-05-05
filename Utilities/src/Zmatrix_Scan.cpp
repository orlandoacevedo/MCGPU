/*!\file
  \brief Structures and functions used to read and write Z-matrix files.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 */
#include "Zmatrix_Scan.h"

Zmatrix_Scan::Zmatrix_Scan(string filename, Opls_Scan* oplsScannerRef){
    fileName = filename;
    oplsScanner = oplsScannerRef;
    startNewMolecule = false;
}

Zmatrix_Scan::~Zmatrix_Scan(){}

int Zmatrix_Scan::scanInZmatrix(){
    stringstream output;
    int numOfLines=0;
    ifstream zmatrixScanner(fileName.c_str());

    if( !zmatrixScanner.is_open() )
        return -1;
    else {
        string line; 
    while( zmatrixScanner.good() )
    {
        numOfLines++;
        getline(zmatrixScanner,line);

        Molecule workingMolecule;

        //check if it is a commented line,
        //or if it is a title line
        try{
            if(line.at(0) != '#' && numOfLines > 1)
            parseLine(line,numOfLines);
        }
        catch(std::out_of_range& e){}

        if (startNewMolecule){
            Atom* atomArray;
            Bond* bondArray;
            Angle* angleArray;
            Dihedral* dihedralArray;
            
            atomArray = (Atom*) malloc(sizeof(Atom) * atomVector.size());
            bondArray = (Bond*) malloc(sizeof(Bond) * bondVector.size());
            angleArray = (Angle*) malloc(sizeof(Angle) * angleVector.size());
            dihedralArray = (Dihedral*) malloc(sizeof(Dihedral) * dihedralVector.size());

            for (int i = 0; i < atomVector.size(); i++){
                atomArray[i] = atomVector[i];
            }
            for (int i = 0; i < bondVector.size(); i++){
                bondArray[i] = bondVector[i];
            }
            for (int i = 0; i < angleVector.size(); i++){
                angleArray[i] = angleVector[i];
            }
            for (int i = 0; i < dihedralVector.size(); i++){
                dihedralArray[i] = dihedralVector[i];
            }

            moleculePattern.push_back(createMolecule(-1, atomArray, angleArray, bondArray, dihedralArray, 
                 atomVector.size(), angleVector.size(), bondVector.size(), dihedralVector.size()));

            atomVector.clear();
            bondVector.clear();
            angleVector.clear();
            dihedralVector.clear();

            startNewMolecule = false;
        } 
    }

    zmatrixScanner.close();
	 
    }

    return 0;
}

void Zmatrix_Scan::parseLine(string line, int numOfLines){

    string atomID, atomType, oplsA, oplsB, bondWith, bondDistance, angleWith, angleMeasure, dihedralWith, dihedralMeasure;

    stringstream ss;

    //check if line contains correct format
    int format = checkFormat(line);

    if(format == 1){
        //read in strings in columns and store the data in temporary variables
        ss << line;    	
        ss >> atomID >> atomType >> oplsA >> oplsB >> bondWith >> bondDistance >> angleWith >> angleMeasure >> dihedralWith >> dihedralMeasure;

        //setup structures for permanent encapsulation
        Atom lineAtom;
        Bond lineBond;
        Angle lineAngle;
        Dihedral lineDihedral;
		  
        if (oplsA.compare("-1") != 0)
        {
            lineAtom = oplsScanner->getAtom(oplsA);
            lineAtom.id = atoi(atomID.c_str());
            lineAtom.x = 0;
            lineAtom.y = 0;
            lineAtom.z = 0;
        }
        else//dummy atom
        {
        	char dummy = 'X';
            lineAtom = createAtom(atoi(atomID.c_str()), -1, -1, -1, -1, -1, -1, dummy);
        }
		  atomVector.push_back(lineAtom);

        if (bondWith.compare("0") != 0){
            lineBond.atom1 = lineAtom.id;
            lineBond.atom2 = atoi(bondWith.c_str());
            lineBond.distance = atof(bondDistance.c_str());
            lineBond.variable = false;
            bondVector.push_back(lineBond);
        }

        if (angleWith.compare("0") != 0){
            lineAngle.atom1 = lineAtom.id;
            lineAngle.atom2 = atoi(angleWith.c_str());
            lineAngle.value = atof(angleMeasure.c_str());
            lineAngle.variable = false;
            angleVector.push_back(lineAngle);
        }

        if (dihedralWith.compare("0") != 0){
            lineDihedral.atom1 = lineAtom.id;
            lineDihedral.atom2 = atoi(dihedralWith.c_str());
            lineDihedral.value = atof(dihedralMeasure.c_str());
            lineDihedral.variable = false;
            dihedralVector.push_back(lineDihedral);
        }
    } //end if format == 1

    else if(format == 2){
        startNewMolecule = true;
	 }
    else if(format == 3){
        startNewMolecule = true;
    }
    if (previousFormat >= 3 && format == -1)
        handleZAdditions(line, previousFormat);

    previousFormat = format;
}
    
int Zmatrix_Scan::checkFormat(string line){
    int format =-1; 
    stringstream iss(line);
    stringstream iss2(line);
    string atomType, someLine;
    int atomID, oplsA, oplsB, bondWith, angleWith,dihedralWith,extra;
    double bondDistance, angleMeasure, dihedralMeasure;	 

    // check if it is the normal 11 line format
    if( iss >> atomID >> atomType >> 
        oplsA >> oplsB >> 
        bondWith >> bondDistance >> 
        angleWith >> angleMeasure >> 
        dihedralWith >> dihedralMeasure >> extra)
        format = 1;
    else{
        someLine = line;
        if(someLine.find("TERZ")!=string::npos){
            format = 2;
		  }
        else if(someLine.find("Geometry Variations")!=string::npos)
            format = 3;
        else if(someLine.find("Variable Bonds")!=string::npos)
            format = 4;
        else if(someLine.find("Additional Bonds")!=string::npos)
            format = 5;
        else if(someLine.find("Harmonic Constraints")!=string::npos)
            format = 6;
        else if(someLine.find("Variable Bond Angles")!=string::npos)
            format = 7;
        else if(someLine.find("Additional Bond Angles")!=string::npos)
            format = 8;
        else if(someLine.find("Variable Dihedrals")!=string::npos)
            format = 9;
        else if(someLine.find("Additional Dihedrals")!=string::npos)
            format = 10;
        else if(someLine.find("Domain Definitions")!=string::npos)
            format = 11;
        else if(someLine.find("Final blank line")!=string::npos)
            format = -2;  
    }

    return format;
}

void Zmatrix_Scan::handleZAdditions(string line, int cmdFormat){
    vector<int> atomIds;
    int id;
    stringstream tss(line.substr(0,15) );

    if(line.find("AUTO")!=string::npos){
	     //Do stuff for AUTO
    }
    else{
        while(tss >> id){
            atomIds.push_back(id);
            if(tss.peek()=='-'||tss.peek()==','||tss.peek()==' ')
                tss.ignore();
        }

        int start, end=0;
        if( atomIds.size()== 1){
            start = atomIds[0];
            end = atomIds[0]; 
        }
        else if(atomIds.size() == 2){
            start = atomIds[0]; end = atomIds[1];
        }

        switch(cmdFormat){
            case 3:
                // Geometry Variations follow 
            break;
            case 4:
                // Variable Bonds follow
                for(int i=0; i< moleculePattern[0].numOfBonds; i++){
                    if(  moleculePattern[0].bonds[i].atom1 >= start &&  moleculePattern[0].bonds[i].atom1 <= end){
                        //cout << "Bond Atom1: "<<  moleculePattern[0].bonds[i].atom1 << " : " <<  moleculePattern[0].bonds[i].variable<<endl;//DEBUG
                        moleculePattern[0].bonds[i].variable = true;
                    }
                }
                break;
            case 5:
                //  Additional Bonds follow 
                break;
            case 6:
                // Harmonic Constraints follow 
                break;
            case 7:
                //  Variable Bond Angles follow
                for(int i=0; i<  moleculePattern[0].numOfAngles; i++){
                    if(  moleculePattern[0].angles[i].atom1 >= start && moleculePattern[0].angles[i].atom1 <= end){
                    //cout << "Angle Atom1: "<<  moleculePattern[0].angles[i].atom1 << " : " << moleculePattern[0].angles[i].variable << endl;//DEBUG
                    moleculePattern[0].angles[i].variable = true;
                    }
                }
                break;
            case 8:
                // Additional Bond Angles follow
                break;
            case 9:
                // Variable Dihedrals follow
                for(int i=0; i< moleculePattern[0].numOfDihedrals; i++){
                    if(  moleculePattern[0].dihedrals[i].atom1 >= start &&  moleculePattern[0].dihedrals[i].atom1 <= end ) {
                    //cout << "Dihedral Atom1: "<<  moleculePattern[0].dihedrals[i].atom1 << " : " <<   moleculePattern[0].dihedrals[i].variable << endl;//DEBUG
                    moleculePattern[0].dihedrals[i].variable = true;
                    }
                }
                break;
            case 10:
                //  Domain Definitions follow
                break;
            default:
                //Do nothing
                break;
        }
    }
}

vector<Hop> Zmatrix_Scan::calculateHops(Molecule molec){
    vector<Hop> newHops;
    int **graph;
    int size = molec.numOfAtoms;
	 int startId = molec.atoms[0].id;

    buildAdjacencyMatrix(graph,molec);

    for(int atom1=0; atom1<size; atom1++){
        for(int atom2=atom1+1; atom2<size; atom2++){
            int distance = findHopDistance(atom1,atom2,size,graph);
            if(distance >=3){
					 Hop tempHop = createHop(atom1+startId,atom2+startId,distance); //+startId because atoms may not start at 1
                newHops.push_back(tempHop);					
            }  		      
        }
    }

    return newHops; 
}

bool Zmatrix_Scan::contains(vector<int> &vect, int item){
    for(int i=0; i<vect.size(); i++){
        if(vect[i]==item)
            return true;
    }

    return false;
}


int Zmatrix_Scan::findHopDistance(int atom1,int atom2,int size, int **graph){
    map<int,int> distance;
    queue<int> Queue;
    vector<int> checked;
    vector<int> bonds;


    Queue.push(atom1);
    checked.push_back(atom1);
    distance.insert( pair<int,int>(atom1,0) );	

    while(!Queue.empty()){
        int target = Queue.front();
        Queue.pop();
        if(target == atom2)
            return distance[target];

        bonds.clear();
        for(int col=0;col<size;col++){
            if( graph[target][col]==1 )
                bonds.push_back(col);
        }

        for(int x=0; x<bonds.size();x++){
            int currentBond = bonds[x];
            if(!contains(checked,currentBond) ){
                checked.push_back(currentBond);
                int newDistance = distance[target]+1;
                distance.insert(pair<int,int>(currentBond, newDistance));
                Queue.push(currentBond);
            }
        }
    }
}

void Zmatrix_Scan::buildAdjacencyMatrix(int **&graph, Molecule molec){
    int size = molec.numOfAtoms;
	 int startId = molec.atoms[0].id; //the first atom ID in the molecule
	 int lastId = startId + molec.numOfAtoms -1; //the last atom ID in the molecule
    graph =  new int*[size]; //create colums
    for(int i=0; i<size; i++) //create rows
        graph[i]=new int[size];	

    //fill with zero
    for(int c=0; c<size; c++)
        for(int r=0; r<size; r++)
            graph[c][r]=0;

    //fill with adjacent array with bonds
    for(int x=0; x<molec.numOfBonds; x++){
        Bond bond = molec.bonds[x];
		  //make sure the bond is intermolecular
		  if( (bond.atom1 >= startId && bond.atom1 <= lastId) &&
		     (bond.atom2 >= startId && bond.atom2 <= lastId) ){
			       graph[bond.atom1-startId][bond.atom2-startId]=1;
                graph[bond.atom2-startId][bond.atom1-startId]=1;
			  }       
    }
}

vector<Molecule> Zmatrix_Scan::buildMolecule(int startingID){
	 int numOfMolec = moleculePattern.size();
	 Molecule newMolecules[numOfMolec];
	 
    //need a deep copy of molecule pattern incase it is modified.
    for (int i = 0; i < moleculePattern.size(); i++){
        Atom *atomCopy = new Atom[ moleculePattern[i].numOfAtoms] ;
        for(int a=0; a <  moleculePattern[i].numOfAtoms ; a++)
            atomCopy[a]=  moleculePattern[i].atoms[a];

        Bond *bondCopy = new Bond[ moleculePattern[i].numOfBonds] ;
        for(int a=0; a <  moleculePattern[i].numOfBonds ; a++)
            bondCopy[a]=  moleculePattern[i].bonds[a];

        Angle *angleCopy = new Angle[ moleculePattern[i].numOfAngles] ;
        for(int a=0; a <  moleculePattern[i].numOfAngles ; a++)
            angleCopy[a]=  moleculePattern[i].angles[a];

        Dihedral *dihedCopy = new Dihedral[ moleculePattern[i].numOfDihedrals];
        for(int a=0; a <  moleculePattern[i].numOfDihedrals ; a++)
            dihedCopy[a]=  moleculePattern[i].dihedrals[a];

        //calculate and add array of Hops to the molecule
        vector<Hop> calculatedHops;
        calculatedHops = calculateHops(moleculePattern[i]);
        int numOfHops = calculatedHops.size();
        Hop *hopCopy = new Hop[numOfHops];
        for(int a=0; a < numOfHops; a++)
            hopCopy[a] = calculatedHops[a];


        Molecule molecCopy = createMolecule(-1,atomCopy, angleCopy, bondCopy, dihedCopy, hopCopy, 
                                    moleculePattern[i].numOfAtoms, 
                                    moleculePattern[i].numOfAngles,
                                    moleculePattern[i].numOfBonds,
                                    moleculePattern[i].numOfDihedrals,
                                    numOfHops);	
		  newMolecules[i] = molecCopy; 
             
    }
				
	 //Assign/calculate the appropiate x,y,z positions to the molecules. 									
	 buildMoleculeInSpace(newMolecules, numOfMolec);
	 

    for (int i = 0; i < numOfMolec; i++)
    {
        if(i == 0){
            newMolecules[i].id = startingID;
        }
        else
            newMolecules[i].id = newMolecules[i-1].id + newMolecules[i-1].numOfAtoms; 
    }
	 
    for (int j = 0; j < numOfMolec; j++)
    {
        Molecule newMolecule = newMolecules[j];
        //map unique IDs to atoms within structs based on startingID
        for(int i = 0; i < newMolecules[j].numOfAtoms; i++){
            int atomID = newMolecule.atoms[i].id - 1;
            //newMolecule.atoms[i].id = atomID + newMolecule.id;
				newMolecule.atoms[i].id = atomID + startingID;

        }
        for (int i = 0; i < newMolecule.numOfBonds; i++){
            int atom1ID = newMolecule.bonds[i].atom1 - 1;
            int atom2ID = newMolecule.bonds[i].atom2 - 1;

            //newMolecule.bonds[i].atom1 = atom1ID + newMolecule.id;
            //newMolecule.bonds[i].atom2 = atom2ID + newMolecule.id;
				newMolecule.bonds[i].atom1 = atom1ID + startingID;
            newMolecule.bonds[i].atom2 = atom2ID + startingID;
        }
        for (int i = 0; i < newMolecule.numOfAngles; i++){
            int atom1ID = newMolecule.angles[i].atom1 - 1;
            int atom2ID = newMolecule.angles[i].atom2 - 1;

            //newMolecule.angles[i].atom1 = atom1ID + newMolecule.id;
            //newMolecule.angles[i].atom2 = atom2ID + newMolecule.id;
				newMolecule.angles[i].atom1 = atom1ID + startingID;
            newMolecule.angles[i].atom2 = atom2ID + startingID;
        }
        for (int i = 0; i < newMolecule.numOfDihedrals; i++){
            int atom1ID = newMolecule.dihedrals[i].atom1 - 1;
            int atom2ID = newMolecule.dihedrals[i].atom2 - 1;

            //newMolecule.dihedrals[i].atom1 = atom1ID + newMolecule.id;
            //newMolecule.dihedrals[i].atom2 = atom2ID + newMolecule.id;
				newMolecule.dihedrals[i].atom1 = atom1ID + startingID;
            newMolecule.dihedrals[i].atom2 = atom2ID + startingID;
        }
        for (int i = 0; i < newMolecule.numOfHops; i++){
            int atom1ID = newMolecule.hops[i].atom1 - 1;
            int atom2ID = newMolecule.hops[i].atom2 - 1;

            //newMolecule.hops[i].atom1 = atom1ID + newMolecule.id;
            //newMolecule.hops[i].atom2 = atom2ID + newMolecule.id;
				newMolecule.hops[i].atom1 = atom1ID + startingID;
            newMolecule.hops[i].atom2 = atom2ID + startingID;
        }
    }

    return vector<Molecule>(newMolecules,newMolecules+sizeof(newMolecules)/sizeof(Molecule) );
}
