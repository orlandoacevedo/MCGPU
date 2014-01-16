#include "zMatrixTest.h"

Molecule createMeshZMolecules(Opls_Scan *scanner){

    //1 S    200    0    0    0.000000   0    0.000000   0    0.000000        0  
    Atom atom1=scanner->getAtom("200");
    atom1.id=1;

    // 2 DUM   -1    0    1    0.500000   0    0.000000   0    0.000000        0
    Atom atom2=createAtom(2,-1,-1,-1,-1,-1,-1,NULL);
    Bond bond2=createBond(2,1,0.5,false);

    //3 DUM   -1    0    2    0.500000   1   90.000000   0    0.000000        0 
    Atom atom3=createAtom(3,-1,-1,-1,-1,-1,-1,NULL);
    Bond bond3=createBond(3,2,0.5,false);
    Angle angle3=createAngle(3,1,90,false);

    //4 hH    204    0    1    1.336532   2   90.000000   3  180.000000        0    
    Atom atom4=scanner->getAtom("204");
    atom4.id=4;
    Bond bond4=createBond(4,1,1.336532,true);
    Angle angle4=createAngle(4,2,90,false);
    Dihedral dihed4=createDihedral(4,3,180,false);

    //5 C    217    0    1    1.811119   4   96.401770   2  180.000000        0
    Atom atom5=scanner->getAtom("217");
    atom5.id=5;
    Bond bond5=createBond(5,1,1.811119,true);
    Angle angle5=createAngle(5,4,96.401770,true);
    Dihedral dihed5=createDihedral(5,2,180,false);

    //6 HC   140    0    5    1.090187   1  110.255589   4  179.999947        0
    Atom atom6=scanner->getAtom("140");
    atom6.id=6;
    Bond bond6=createBond(6,5,1.090187,true);
    Angle angle6=createAngle(6,1,110.255589,true);
    Dihedral dihed6=createDihedral(6,4,179.999947,true);    

    //7 HC   140    0    5    1.090135   6  108.527646   1  121.053891        0
    Atom atom7=scanner->getAtom("140");
    atom7.id=7;
    Bond bond7=createBond(7,5,1.090135,true);
    Angle angle7=createAngle(7,6,108.527646,true);
    Dihedral dihed7=createDihedral(7,1,121.053891,true);    

    //8 HC   140    0    5    1.090135   6  108.527646   1  238.946114        0
    Atom atom8=scanner->getAtom("140");
    atom8.id=8;
    Bond bond8=createBond(8,5,1.090135,true);
    Angle angle8=createAngle(8,6,108.527646,true);
    Dihedral dihed8=createDihedral(8,1,238.946114,true);

    /* HOPS and BONDS Diagram of Mesh.z
    //      1--2--3
    //      |\
    //      | \
    //      4  5
    //        /|\
    //       / | \
    //      6  7  8
    */
    //All hops that have a hop distance >= 3

    Hop hop1= createHop(2,6,3);
    Hop hop2= createHop(2,7,3);
    Hop hop3= createHop(2,8,3);
    Hop hop4= createHop(3,4,3);
    Hop hop5= createHop(3,5,3);
    Hop hop6= createHop(3,6,4);
    Hop hop7= createHop(3,7,4);
    Hop hop8= createHop(3,8,4);
    Hop hop9= createHop(4,6,3);
    Hop hop10= createHop(4,7,3);
    Hop hop11= createHop(4,8,3);

    Atom *atomPtr = new Atom[8];
    Bond *bondPtr = new Bond[7];
    Angle *anglePtr = new Angle[6];
    Dihedral *dihedPtr = new Dihedral[5];
    Hop *hopPtr = new Hop[11];


    atomPtr[0]=atom1;  atomPtr[1]=atom2; atomPtr[2]=atom3; atomPtr[3]=atom4; atomPtr[4]=atom5; atomPtr[5]=atom6; atomPtr[6]=atom7; atomPtr[7]=atom8;
    bondPtr[0]=bond2; bondPtr[1]=bond3; bondPtr[2]=bond4; bondPtr[3]=bond5; bondPtr[4]=bond6; bondPtr[5]=bond7; bondPtr[6]=bond8;
    anglePtr[0]=angle3; anglePtr[1]=angle4; anglePtr[2]=angle5; anglePtr[3]=angle6; anglePtr[4]=angle7; anglePtr[5]=angle8;
    dihedPtr[0]=dihed4; dihedPtr[1]=dihed5; dihedPtr[2]=dihed6; dihedPtr[3]=dihed7; dihedPtr[4]=dihed8;
    hopPtr[0]=hop1; hopPtr[1]=hop2; hopPtr[2]=hop3; hopPtr[3]=hop4; hopPtr[4]=hop5; hopPtr[5]=hop6;  hopPtr[6]=hop7; hopPtr[7]=hop8; hopPtr[8]=hop9; 
    hopPtr[9]=hop10; hopPtr[10]=hop11;


    return createMolecule(1,atomPtr,anglePtr,bondPtr,dihedPtr,hopPtr,8,6,7,5,11);
}

void addPositionsToMeshZ(Molecule *meshMolec){
    //the hardcoded positions(x y z) of the mezh.z molecule
    // Atom#   X    Y    Z
	 // 0       0    0    0
	 meshMolec->atoms[0].x = 0;
	 meshMolec->atoms[0].y = 0;
    meshMolec->atoms[0].z = 0;
	 // 1       0    0.5    0
	 meshMolec->atoms[1].x = 0;
	 meshMolec->atoms[1].y = 0.5;
    meshMolec->atoms[1].z = 0;
	 // 2     -0.5    0.5   -0.5
	 meshMolec->atoms[2].x = -0.5;
	 meshMolec->atoms[2].y = 0.5;
    meshMolec->atoms[2].z = -0.5;
	 // 3    1.33653   2.39894e-9   1.33653
	 meshMolec->atoms[3].x = 1.33653;
	 meshMolec->atoms[3].y = 0.00000000239894;
    meshMolec->atoms[3].z = 1.33653;
	 // 4   -0.857136    -1.36313   -0.823996
	 meshMolec->atoms[4].x = -0.857136 ;
	 meshMolec->atoms[4].y = -1.36313;
    meshMolec->atoms[4].z = -0.823996;
	 // 5   -1.92671    -1.51275   -0.405374
	 meshMolec->atoms[5].x = -1.92671;
	 meshMolec->atoms[5].y = -1.51275;
    meshMolec->atoms[5].z = -0.405374;
	 // 6   -0.878138    -1.16034   -1.94557
	 meshMolec->atoms[6].x = -0.878138;
	 meshMolec->atoms[6].y = -1.16034;
    meshMolec->atoms[6].z = -1.94557;
	 // 7    -0.244983    -2.314   -0.710251
	 meshMolec->atoms[7].x = -0.244983;
	 meshMolec->atoms[7].y = -2.314;
    meshMolec->atoms[7].z = -0.710251;
}

vector<Molecule> createT3pdimMolecules(Opls_Scan *scanner){
   // TIP3P Water Dimer                                      Tot. E =     -6.5396
   //1  O   111  111    0     .000000   0     .000000   0     .000000        0
	Atom atom1 = scanner->getAtom("111");
	atom1.id=1;    
	
   //2  DU   -1   -1    1    1.000000   0     .000000   0     .000000        0 
	Atom atom2=createAtom(2,-1,-1,-1,-1,-1,-1,NULL);
   Bond bond2=createBond(2,1,1.0,false);   
	
   //3  DU   -1   -1    2    1.000000   1   90.000000   0     .000000        0
	Atom atom3=createAtom(3,-1,-1,-1,-1,-1,-1,NULL);
   Bond bond3=createBond(3,2,1.0,false);
	Angle angle3=createAngle(3,1,90.0,false);    
	
   //4  HO  112  112    1     .957200   2  127.740000   3   90.000000        0
	Atom atom4=scanner->getAtom("112");
	atom4.id=4;
	Bond bond4=createBond(4,1,.957200,false);
	Angle angle4=createAngle(4,2,127.740000,false);
	Dihedral dihed4=createDihedral(4,3,90.00,false);
	    
   //5  HO  112  112    1     .957200   4  104.520000   2  180.000000        0 
	Atom atom5=scanner->getAtom("112");
	atom5.id=5;
	Bond bond5=createBond(5,1,.957200,false);
	Angle angle5=createAngle(5,4,104.520000,false);
	Dihedral dihed5=createDihedral(5,2,180.00,false); 
	     
   //6  X    -1   -1    1     .150000   4   52.260000   5     .000000        0
	Atom atom6=createAtom(6,-1,-1,-1,-1,-1,-1,NULL);
	Bond bond6=createBond(6,1,.150000,false);
	Angle angle6=createAngle(6,4,52.260000,false);
	Dihedral dihed6=createDihedral(6,5,.00,false); 
	
	//TERZ 2nd Bonded Molecule                                                                           
   //7  O   111  111    1    2.751259   2  131.716186   3  269.946364        0
	Atom atom7=scanner->getAtom("111");
	atom7.id=7;
	Bond bond7=createBond(7,1,2.751259,false);
	Angle angle7=createAngle(7,2,131.716186,false);
	Dihedral dihed7=createDihedral(7,3,269.946364,false);
	    
   //8  DU   -1   -1    7    1.000000   1   21.472391   2  179.887530        0
	Atom atom8=createAtom(8,-1,-1,-1,-1,-1,-1,NULL);
	Bond bond8=createBond(8,7,1.0,false);
	Angle angle8=createAngle(8,1,21.472391,false);
	Dihedral dihed8=createDihedral(8,2,179.887530,false); 
	    
   //9  DU   -1   -1    8    1.000000   7   90.000000   1  179.624613        0
	Atom atom9=createAtom(9,-1,-1,-1,-1,-1,-1,NULL);
	Bond bond9=createBond(9,8,1.0,false);
	Angle angle9=createAngle(9,7,90.000000,false);
	Dihedral dihed9=createDihedral(9,1,179.624613,false); 
	    
   //10  HO  112  112    7     .957200   8  127.740000   9   90.000000        0
	Atom atom10=scanner->getAtom("112");
	atom10.id=10;
	Bond bond10=createBond(10,7,0.957200,false);
	Angle angle10=createAngle(10,8,127.74,false);
	Dihedral dihed10=createDihedral(10,9,90.0,false);

    
   //11  HO  112  112    7     .957200  10  104.520000   8  180.000000        0
	Atom atom11=scanner->getAtom("112");
	atom11.id=11;
	Bond bond11=createBond(11,7,0.957200,false);
	Angle angle11=createAngle(11,10,104.52,false);
	Dihedral dihed11=createDihedral(11,8,180.0,false);
	    
   //12  X    -1   -1    7     .150000  10   52.260000  11     .000000        0 
	Atom atom12=createAtom(12,-1,-1,-1,-1,-1,-1,NULL);
	Bond bond12=createBond(12,7,0.15,false);
	Angle angle12=createAngle(12,10,52.26,false);
	Dihedral dihed12=createDihedral(12,11,0.0,false);
	
	/* HOPS and BONDS Diagram of Mesh.z
   //      
	//
	//      10  11 12
	//        \ | /
	//         \|/
	//          7--8--9
	//          |
	//          |
	//          1--2--3
	//         /|\
	//        / | \
	//       4  5  6
	*/
	//All intermolecular hops that have a hop distance >= 3 
	
	//bottom Molecule 
	Hop hop1=createHop(3,4,3);
	Hop hop2=createHop(3,5,3);
	Hop hop3=createHop(3,6,3);
	//top Molecule 
	Hop hop4=createHop(9,10,3);
	Hop hop5=createHop(9,11,3);
	Hop hop6=createHop(9,12,3);
	
	Atom *atomPtr = new Atom[6];
   Bond *bondPtr = new Bond[5];
   Angle *anglePtr = new Angle[4];
   Dihedral *dihedPtr = new Dihedral[3];
   Hop *hopPtr = new Hop[3];
	vector<Molecule> retVect;
	
	atomPtr[0]=atom1; atomPtr[1]=atom2; atomPtr[2]=atom3; atomPtr[3]=atom4; atomPtr[4]=atom5; atomPtr[5]=atom6;
	bondPtr[0]=bond2; bondPtr[1]=bond3; bondPtr[2]=bond4; bondPtr[3]=bond5; bondPtr[4]=bond6;
	anglePtr[0]=angle3; anglePtr[1]=angle4; anglePtr[2]=angle5; anglePtr[3]=angle6;
	dihedPtr[0]=dihed4; dihedPtr[1]=dihed5; dihedPtr[2]=dihed6;
	hopPtr[0]=hop1; hopPtr[1]=hop2; hopPtr[2]=hop3;  
	retVect.push_back( createMolecule(1,atomPtr,anglePtr,bondPtr,dihedPtr,hopPtr,6,4,5,3,3) );
	
	Atom *atomPt = new Atom[6];
   Bond *bondPt = new Bond[6];
   Angle *anglePt = new Angle[6];
   Dihedral *dihedPt = new Dihedral[6];
   Hop *hopPt = new Hop[3];
	
	atomPt[0]=atom7; atomPt[1]=atom8; atomPt[2]=atom9; atomPt[3]=atom10; atomPt[4]=atom11; atomPt[5]=atom12;
	bondPt[0]=bond7; bondPt[1]=bond8; bondPt[2]=bond9; bondPt[3]=bond10; bondPt[4]=bond11;  bondPt[5]=bond12;
	anglePt[0]=angle7; anglePt[1]=angle8; anglePt[2]=angle9; anglePt[3]=angle10; anglePt[4]=angle11; anglePt[5]=angle12;
	dihedPt[0]=dihed7; dihedPt[1]=dihed8; dihedPt[2]=dihed9; dihedPt[3]=dihed10; dihedPt[4]=dihed11; dihedPt[5]=dihed12;
	hopPt[0]=hop4; hopPt[1]=hop5; hopPt[2]=hop6;  
	retVect.push_back( createMolecule(7,atomPt,anglePt,bondPt,dihedPt,hopPt,6,6,6,6,3) );
	
	return retVect;
}

void addPositionsToT3pdim(vector<Molecule> t3pdimMolec){
   //the hardcoded positions(x y z) of the t3pdim.z molecule
   // Atom#   X    Y    Z
   //1, 0, 0, 0
	t3pdimMolec[0].atoms[0].x=0;
	t3pdimMolec[0].atoms[0].y=0;
	t3pdimMolec[0].atoms[0].z=0;
	
	//2, 0, 1, 0
	t3pdimMolec[0].atoms[1].x=0;
	t3pdimMolec[0].atoms[1].y=1;
	t3pdimMolec[0].atoms[1].z=0;
	
	//3, -1, 1, -1
	t3pdimMolec[0].atoms[2].x=-1;
	t3pdimMolec[0].atoms[2].y=1;
	t3pdimMolec[0].atoms[2].z=-1;
	
	//4, 0.862904, -0.585882, -0.862904
	t3pdimMolec[0].atoms[3].x=0.862904;
	t3pdimMolec[0].atoms[3].y=-0.585882;
	t3pdimMolec[0].atoms[3].z=-0.862904;
	
	//5, 2.92389, -0.142182, -0.0567901
	t3pdimMolec[0].atoms[4].x=2.92389;
	t3pdimMolec[0].atoms[4].y=-0.1412182;
	t3pdimMolec[0].atoms[4].z=-0.0567901;
	
	//6, 0.100926, 0.0671973, -0.00196027
	t3pdimMolec[0].atoms[5].x=0.100926;
	t3pdimMolec[0].atoms[5].y=0.0671973;
	t3pdimMolec[0].atoms[5].z=-0.00196027;
	
	//7, -2.42538, -1.8308, 2.42993
	t3pdimMolec[1].atoms[0].x=-2.42538;
	t3pdimMolec[1].atoms[0].y=-1.8308;
	t3pdimMolec[1].atoms[0].z=2.42993;
	
	//8, 1.23403, 0.714883, -0.186453
	t3pdimMolec[1].atoms[1].x=1.23403;
	t3pdimMolec[1].atoms[1].y=0.714883;
	t3pdimMolec[1].atoms[1].z=-0.186453;

	//9, 5.08911, 0.528186, -3.5686
	t3pdimMolec[1].atoms[2].x=5.08911;
	t3pdimMolec[1].atoms[2].y=0.528186;
	t3pdimMolec[1].atoms[2].z=-3.5686;

	//10, 0.859175, -5.50301, 4.27592
	t3pdimMolec[1].atoms[3].x=0.859175;
	t3pdimMolec[1].atoms[3].y=-5.50301;
	t3pdimMolec[1].atoms[3].z=4.27592;

	//11, 0.621781, 2.83868, 0.176906
	t3pdimMolec[1].atoms[4].x=0.621781;
	t3pdimMolec[1].atoms[4].y=2.83868;
	t3pdimMolec[1].atoms[4].z=0.176906;

	//12, 0.210078, -1.31891, -0.204393
	t3pdimMolec[1].atoms[5].x=0.210078;
	t3pdimMolec[1].atoms[5].y=-1.31891;
	t3pdimMolec[1].atoms[5].z=-0.204393;

}  

void compareTestMolecules(Molecule molec1, Molecule molec2){
    // check if id's are equal
    assert(molec1.id == molec2.id);	 

    //check atoms
    assert(molec1.numOfAtoms == molec2.numOfAtoms);
    for(int i=0; i< molec1.numOfAtoms; i++){
	     assert(molec1.atoms[i].id==molec2.atoms[i].id);
        assert(percentDifference(molec1.atoms[i].charge, molec2.atoms[i].charge));
        assert(percentDifference(molec1.atoms[i].sigma, molec2.atoms[i].sigma));
        assert(percentDifference(molec1.atoms[i].epsilon, molec2.atoms[i].epsilon));
        assert(percentDifference(molec1.atoms[i].x, molec2.atoms[i].x));
        assert(percentDifference(molec1.atoms[i].y, molec2.atoms[i].y));
        assert(percentDifference(molec1.atoms[i].z, molec2.atoms[i].z));			
    }

    //check bond
    assert(molec1.numOfBonds == molec2.numOfBonds);
    for(int i=0; i< molec1.numOfBonds; i++){
        assert(molec1.bonds[i].atom1 == molec2.bonds[i].atom1);
        assert(molec1.bonds[i].atom2 == molec2.bonds[i].atom2);
        assert(percentDifference(molec1.bonds[i].distance, molec2.bonds[i].distance));
        assert(asserTwoBool( molec1.bonds[i].variable,molec2.bonds[i].variable));
    }

    //check angles
    assert(molec1.numOfAngles == molec2.numOfAngles);
    for(int i=0; i< molec1.numOfAngles; i++){
        assert(molec1.angles[i].atom1 == molec2.angles[i].atom1);
        assert(molec1.angles[i].atom2 == molec2.angles[i].atom2);
        assert(percentDifference(molec1.angles[i].value,molec2.angles[i].value));
        assert(asserTwoBool(molec1.angles[i].variable, molec2.angles[i].variable));
    }

    //check dihederals
    assert(molec1.numOfDihedrals == molec2.numOfDihedrals);
    for(int i=0; i< molec1.numOfDihedrals; i++){
        assert(molec1.dihedrals[i].atom1 == molec2.dihedrals[i].atom1);
        assert(molec1.dihedrals[i].atom2 == molec2.dihedrals[i].atom2);
        assert(percentDifference(molec1.dihedrals[i].value,molec2.dihedrals[i].value));
        assert(asserTwoBool(molec1.dihedrals[i].variable,molec2.dihedrals[i].variable));
    }

    //check hops
    assert(molec1.numOfHops == molec2.numOfHops);
    for(int i=0; i< molec1.numOfHops; i++){
        assert(molec1.hops[i].atom1 == molec2.hops[i].atom1);
        assert(molec1.hops[i].atom2 == molec2.hops[i].atom2);
        assert(molec1.hops[i].hop == molec2.hops[i].hop);
    }

}


void testZmatrixScanner(Opls_Scan *opls){
    cout << "Testing Z-matrix scanner"<< endl;
    string zMatrixFile1 = "Utilities/bossFiles/mesh.z";
    Molecule meshZ= createMeshZMolecules(opls);	 
	 addPositionsToMeshZ(&meshZ);
	 
    Zmatrix_Scan zScan (zMatrixFile1,opls);
    vector<Molecule> scannedInMolecules;

    int open = zScan.scanInZmatrix();
    if(open == -1)
        cout << "Zmatrix file: " << zMatrixFile1 << "Failed to Open" << endl;
    else{
        scannedInMolecules = zScan.buildMolecule(1);
    }

    assert(scannedInMolecules.size()==1);
    compareTestMolecules(scannedInMolecules[0],meshZ);  
    cout << "Testing Z-matrix scanner Completed\n"<< endl;  
}

void testZmatrixScanner_multipleSingle(Opls_Scan *opls){
    cout << "Testing Z-matrix scanner with reuse 100 molec : 1 atom each"<< endl;
    string zMatrixFile1 = "Utilities/bossFiles/testZ.z";
    Zmatrix_Scan zScan (zMatrixFile1,opls);
    vector<Molecule> scannedInMolecules;
    
    cout << "Zmatrix_Scan created." << endl;

    int open = zScan.scanInZmatrix();
    if(open == -1){
        cout << "Zmatrix file: " << zMatrixFile1 << "Failed to Open" << endl;
        assert(open >= 0 );
    }
    
    cout << "ZmatrixFile opened and read." << endl;

    Molecule molec;
    for(int i=0; i<100; i++){
        scannedInMolecules = zScan.buildMolecule(i);
        assert(scannedInMolecules.size()==1);
        molec = scannedInMolecules[0];
        assert(molec.id == i);
        assert(molec.numOfAtoms==1);
        assert(molec.atoms[0].id == i);
    }
     
    cout << "Testing Z-matrix scanner with reuse 100 molec : 1 atom each Complete\n"<< endl;	 
}

void testZmatrixScanner_BondedMolecules(Opls_Scan *opls){
    cout <<"Testing Z-matrix scanner with Bonded molecules: t3pdim"<< endl;
	 string zMatrixFile1 =  "Utilities/bossFiles/t3pdim.z";
    Zmatrix_Scan zScan (zMatrixFile1,opls);
    cout << "Zmatrix_Scan created" << endl;
    vector<Molecule> t3pdim= createT3pdimMolecules(opls) ;
    cout << "t3pdim created." << endl;	
    addPositionsToT3pdim(t3pdim);
    int open = zScan.scanInZmatrix();
    if(open == -1){
        cout << "Zmatrix file: " << zMatrixFile1 << "Failed to Open" << endl;
        assert(open >= 0 );
    }
	 vector<Molecule> scannedInMolecules;
	 scannedInMolecules = zScan.buildMolecule(1);
	 assert(scannedInMolecules.size()==t3pdim.size());
	 for(int i=0; i<t3pdim.size(); i++){
	     cout<<"Comparing Molecule: "<<i+1<<endl;
	     compareTestMolecules(scannedInMolecules[i],t3pdim[i]);
		  }	
    cout <<"Testing Z-matrix scanner with Bonded molecules: t3pdim Complete\n"<< endl;
}

void testZmatrix(Opls_Scan *scanner){
    //Test Zmatrix
    testZmatrixScanner(scanner); 
	testZmatrixScanner_multipleSingle(scanner);
	testZmatrixScanner_BondedMolecules(scanner);
}
