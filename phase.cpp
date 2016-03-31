#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <float.h>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <vector>
#include "randomc.h"
#include "ranrotw.cpp"

#include <mpi.h>

using namespace std;

double sign(double v)
{
double factor;
if (v<= 0) { factor=-1;}
else { factor=1; }
return factor;
}

int rownd (double a) {
	return(int(a+ 0.5));
}

long int rowndlong (double a) {
	return((long int)(a+0.5));
}

double distance(int ii, int jj, int x, int y)
{
double dist;
dist = sqrt(pow(double(ii)-double(x)-double(((jj+y)%2)*(0.5*double(jj%2)-0.5*double(y%2))),2)+pow(0.8660254*(double(jj)-double(y)),2));
return dist;      
}


double prob1 ( int a, int b, double selin ) {

	double p1result;
	double p1num;
	double p1den;
	if (a==1){ p1num=1;} else {
	p1num = pow(1-selin,double(a-1));
				  }
	if (b==1){ p1den=1;} else {
	p1den = pow(1-selin,double(b-1));
				  }
	p1result = p1num/(p1num+p1den);

	return p1result;
}



int main(int argc, char *argv[])
{

	MPI::Init(argc,argv);

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	int root = 0;

    
    const double EPSILON=0.0000000001;
 
    string statsfilename;  
    ofstream outstats;
    ifstream testoutstats;
	int filecounter=0;
	string tempstr;
   
    int32 seed = time(0)+rownd(double(rank)*double(time(0))/double(size));
    
    int leftedge;
    int rightedge;
    int numruns;
    double numrunsd;
    int it;
    double sectorwidth;
    
    double rannum, rannum1, rannum2;
    double selection, mutf, mutb;

    int latticesizex;
    double latticesizexd;
    int startseed;
    int heteroyes;
    int numcolors;
    double basegrowth=0;
    
    long int latticesizet;
    int nstep, nslices, ncounter;
    int nslices2;

	ostringstream tempstring;	


if (rank == 0) {

    cin >> latticesizex >> latticesizet >> nslices >> numruns >> startseed >> numcolors;

    tempstring << numcolors;
    statsfilename += "nc";
    statsfilename += tempstring.str();
    tempstring.str("");
    tempstring.clear();
    tempstring << startseed;
    statsfilename += "ss";
    statsfilename += tempstring.str();

	tempstring.str("");
	tempstring.clear();
	tempstring << filecounter;	
	
	statsfilename += "run";
	statsfilename += tempstring.str();

	testoutstats.open(statsfilename.c_str());
	testoutstats.close();
	
	while (!testoutstats.fail())
{
	tempstr = tempstring.str();
	statsfilename.erase(statsfilename.end()-tempstr.size(),statsfilename.end());
	filecounter++;
	tempstring.str("");
	tempstring.clear();
	tempstring << filecounter;
	statsfilename += tempstring.str();
	testoutstats.open(statsfilename.c_str());
	testoutstats.close();

}
	testoutstats.clear(ios::failbit);
	outstats.open(statsfilename.c_str());


	cout << statsfilename.c_str() << endl;
	cout << "# Lx: " << latticesizex << " Lt: " << latticesizet << endl;
	cout << "# Nruns: " << numruns << " Nslices: " << nslices << " StartS: " << startseed << " ";
	cout << "# Colors: " << numcolors << endl;
	


    outstats << "# Lx: " << latticesizex << " Lt: " << latticesizet << endl;
    outstats << "# Nruns: " << numruns << " Nslices: " << nslices << " StartS: " << startseed << " ";
    outstats << "# Colors: " << numcolors << endl;


}

MPI::COMM_WORLD.Barrier();

	MPI::COMM_WORLD.Bcast(&latticesizex, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&latticesizet, 1, MPI::LONG, root);
	MPI::COMM_WORLD.Bcast(&numruns, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&nslices, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&startseed, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&numcolors, 1, MPI::INT, root);


MPI::COMM_WORLD.Barrier();


	nslices2=nslices*nslices;

    int centerx = rownd(latticesizex/2);
    int centery = rownd(latticesizet/2);
    numrunsd = double(numruns);
    latticesizexd=double(latticesizex);
    

    TRanrotWGenerator rg(seed);
    
    vector< int > latteven(latticesizex);
    vector< int > lattodd(latticesizex);
 


    vector< double > sdencount;
    vector< double > scolorprobperc;
    vector< double > scolorprobperctemp(numcolors);

    vector< double > dencount;
    vector< double > colorprobperc;
    vector< double > colorprobperctemp(numcolors);


    

            scolorprobperc.clear();
            sdencount.clear();

            colorprobperc.clear();
            dencount.clear();


           for (int j=0; j<numcolors; j++) {
               scolorprobperc.push_back(0);
               sdencount.push_back(0);
                }    
                           
           for (int j=0; j<numcolors; j++) {
               colorprobperc.push_back(0);
               dencount.push_back(0);
                }    
     
mutf=0;
ncounter=0;
for (int t1=0; t1<nslices; t1++) {

selection=0;
for (int t2=0; t2<nslices; t2++) {
      
           for (int j=0; j<numcolors; j++) {
               scolorprobperc[j]=0;
               sdencount[j]=0;
                }    
                           
           for (int j=0; j<numcolors; j++) {
               colorprobperc[j]=0;
               dencount[j]=0;
                }    
 
      
for (int run=1; run<=numruns; run++) {
    
  
  
  if(startseed>0) {
          leftedge=rownd(latticesizex/2)-rownd(startseed/2);
          rightedge=rownd(latticesizex/2)-rownd(startseed/2)+startseed;
          for (int i1=0; i1 < leftedge; i1++) latteven[i1]=numcolors;        
          for (int i2=leftedge; i2 < rightedge; i2++) latteven[i2]=1;         
          for (int i3=rightedge; i3 < latticesizex; i3++) latteven[i3]=numcolors;
                  }
  else {
 if (mutf>EPSILON){     for (int i=0; i < latticesizex; i++) latteven[i]=1; }
 else { for (int i=0; i <latticesizex; i++) latteven[i]=rg.IRandom(1,2); }
       }
       
   for (int t=1; t<=latticesizet; t++) {        
       
       // Odd steps:
              rannum1=rg.Random();
              rannum2=rg.Random();
       if (latteven[0]!=latteven[latticesizex-1]) {
	if (rannum1 < prob1(latteven[0],latteven[latticesizex-1],selection))
                 {
                      lattodd[0]=latteven[0];      
                            } else {
                                   lattodd[0]=latteven[latticesizex-1];
                                   }
                                                  }      else {lattodd[0]=latteven[0];} 
        if(rannum2<mutf && lattodd[0]<numcolors) lattodd[0]++;

                   
       for (int i=1; i < latticesizex; i++)
       {
              rannum1=rg.Random();
              rannum2=rg.Random();    
             
           if (latteven[i]!=latteven[i-1]) {
                    if(rannum1 < prob1(latteven[i],latteven[i-1],selection))
                    {
                         lattodd[i]=latteven[i];      
                               } else { 
                                      lattodd[i]=latteven[i-1];
                                       }
                                                         }
           else {lattodd[i]=latteven[i];}
           
           if(rannum2<mutf && lattodd[i]<numcolors) lattodd[i]++;
           
           }

       t++;
       
       // Even steps:
      for (int i=0; i < latticesizex; i++) 
       {
              rannum1=rg.Random();
              rannum2=rg.Random();
                             
          if (lattodd[i]!=lattodd[(i+1)%latticesizex]) {
                if (rannum1 < prob1(lattodd[i],lattodd[(i+1)%latticesizex],selection))
                              {
                      latteven[i]=lattodd[i];                                                                                                                     
                            } else {
                                   latteven[i]=lattodd[(i+1)%latticesizex];
                                   }                          
                                          }  else 
                            {
                               latteven[i]=lattodd[i];           
                                          } 
          if (rannum2 < mutf && latteven[i]<numcolors) latteven[i]++;                                   
       }

			}
                   it=0;
                   sectorwidth=0; 
                   for(int j=0; j < numcolors; j++) colorprobperctemp[j]=0;            
                   for(int i=0; i < latticesizex; i++) {            
                         dencount[lattodd[(i+it)%latticesizex]-1]+=1/latticesizexd/numrunsd; 
                        colorprobperctemp[lattodd[(i+it)%latticesizex]-1]=1;  
                               }  
                   for(int j=0; j < numcolors; j++) colorprobperc[j]+=colorprobperctemp[j]/numrunsd;            
                                                   }


		for (int j=0; j < numcolors; j++) {
			dencount[j] = dencount[j]/double(size);
		        colorprobperc[j] = colorprobperc[j]/double(size);
		}


MPI::COMM_WORLD.Barrier();


	MPI::COMM_WORLD.Reduce(&dencount.front(),&sdencount.front() , numcolors , MPI::DOUBLE, MPI::SUM, root);
	MPI::COMM_WORLD.Reduce(&colorprobperc.front(),&scolorprobperc.front() , numcolors , MPI::DOUBLE, MPI::SUM, root);

MPI::COMM_WORLD.Barrier();

            if (rank==0) { 
                   outstats << selection << " " << mutf;
                 for (int j=0; j<numcolors; j++) {
                          if (startseed==0) {          
   outstats << " " <<  sdencount[j]; 
					   }
                          if (startseed>0) {
  outstats << " " << scolorprobperc[j];
						 }
                     }
                   outstats << endl;   
              }      


        ncounter++;
	selection=selection+1/double(nslices-1);
	}
	mutf=mutf+0.5/double(nslices-1);
	}



MPI::COMM_WORLD.Barrier();


if (rank==0)  outstats.close();
                   


MPI::COMM_WORLD.Barrier();             

	MPI::Finalize();
  
    return( 0 );
}

