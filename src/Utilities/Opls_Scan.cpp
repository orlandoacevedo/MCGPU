/*!\file
  \brief Structures and functions used to read and write OPLS files.
  \author Alexander Luchs, Riley Spahn, Seth Wooten
 
 */

#include "Opls_Scan.h"
#include <exception>
#include <stdexcept>

Opls_Scan::Opls_Scan(string filename)
{
   fileName = filename;
}

Opls_Scan::~Opls_Scan()
{
    oplsTable.clear();
}

int Opls_Scan::scanInOpls(string filename)
{
    int numOfLines=0;
    ifstream oplsScanner(filename.c_str());
    if( !oplsScanner.is_open() )
        return -1;
    else {
        string line; 
        while( oplsScanner.good() )
        {
            numOfLines++;
            getline(oplsScanner,line);

            //check if it is a commented line,
            //or if it is a title line
            try{
                if(line.at(0) != '#' && numOfLines > 1)
                    addLineToTable(line,numOfLines);
            }
            catch (std::out_of_range& e){}
        }
        oplsScanner.close();
        logErrors();
    }
    return 0;
}

void Opls_Scan::addLineToTable(string line, int numOfLines)
{
    string hashNum;
    int secCol;
    double charge,sigma,epsilon;
    string name, extra;
    stringstream ss(line);

    //check to see what format it is opls, V value, or neither
    int format = checkFormat(line);
              
    if(format == 1)
    {      	
        ss >> hashNum >> secCol >> name >> charge >> sigma >> epsilon;
        char *atomtype = (char*)name.c_str(); 
         
        Atom temp = createAtom(0, -1, -1, -1, sigma, epsilon, charge, *atomtype);
        pair<map<string,Atom>::iterator,bool> ret;
        ret = oplsTable.insert( pair<string,Atom>(hashNum,temp) );

        if (ret.second==false)
        {
            errHashes.push_back(hashNum);
        }
    }
    else if(format == 2)
    {
        double v0,v1,v2,v3;
        ss >> hashNum >> v0 >> v1 >> v2 >> v3 ;
        Fourier vValues = {v0,v1,v2,v3};
        pair<map<string,Fourier>::iterator,bool> ret2;
        ret2 = fourierTable.insert( pair<string,Fourier>(hashNum,vValues) );

        if (ret2.second==false)
        {
            errHashesFourier.push_back(hashNum);	
        }	  
    }
    else
    {
	     errLines.push_back(numOfLines);
    }
}

int Opls_Scan::checkFormat(string line)
{   	 
    int hashNum, secCol;
    double charge,sigma,epsilon;
    string name, extra;
    stringstream iss(line);

    double v1,v2,v3,v4;
    stringstream issw(line);

    //see if format is the V values for the diherdral format
    if((issw >> hashNum >> v1 >> v2 >> v3 >> v4) )
    {
        return 2;
    }
    //else see if format is normal opls format
    else if((iss >> hashNum >> secCol >> name >> charge >> sigma >> epsilon ))
    {
        return 1;
    }
    //if neither return -1
    else
    {
        return -1;
    }
}

void Opls_Scan::logErrors()
{
    stringstream output;
    // See if there were any errors
    if(errLines.empty() || errHashes.empty()|| errHashesFourier.empty())
    {
	     //Errors in the format
		  output<<"Errors found in the OPLS file: "<< fileName<<endl;
        if(!errLines.empty())
        {
		      output << "Found Errors in the Format of the following Lines: " << endl;
				for(int a=0; a<errLines.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << endl;
                    }
				    output << errLines[a]<< " ";
				}
				output << endl<< endl;
		  }
		  if(!errHashes.empty())
          {
		      output << "Error - The following OPLS values existed more than once: " << endl;
				for(int a=0; a<errHashes.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << endl;
                    }
				    output << errHashes[a]<< " ";
				}
				output << endl<< endl;
		  }
		  if(!errHashesFourier.empty())
          {
		      output << "Error - The following Fourier Coefficent values existed more than once: " << endl;
				for(int a=0; a<errHashesFourier.size(); a++)
                {
				    if(a%10==0 && a!=0) //ten per line
                    {
					     output << endl;
                    }
				    output << errHashesFourier[a]<< " ";
				}
				output << endl<< endl;
		  }
		  writeToLog(output,OPLS);
	}
}

Atom Opls_Scan::getAtom(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        return oplsTable[hashNum];
	}
	else
    {
	    cerr << "Index does not exist: "<< hashNum <<endl;
		return createAtom(0, -1, -1, -1, -1, -1, -1, NULL);
	}
}

double Opls_Scan::getSigma(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.sigma;
    }
    else
    {
        cerr << "Index does not exist: "<< hashNum <<endl;
        return -1;
    }
}

double Opls_Scan::getEpsilon(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.epsilon;
    }
    else
    {
        cerr << "Index does not exist: "<< hashNum <<endl;
        return -1;
    }
}

double Opls_Scan::getCharge(string hashNum)
{
    if(oplsTable.count(hashNum)>0 )
    {
        Atom temp = oplsTable[hashNum];
        return temp.charge;
    }
    else
    {
        cerr << "Index does not exist: "<< hashNum <<endl;
        return -1;
    }
}

Fourier Opls_Scan::getFourier(string hashNum)
{
    if(fourierTable.count(hashNum)>0 )
    {
        Fourier temp = fourierTable[hashNum];
        return temp;
    }
    else
    {	    
        cerr << "Index does not exist: "<< hashNum <<endl;
        Fourier temp ={-1,-1,-1,-1};
        return temp;
    }
}
