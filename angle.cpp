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



double distance(int ii, int jj, int x, int y)
{
double dist;
dist = sqrt(pow(double(ii)-double(x)-double(((jj+y)%2)*(0.5*double(jj%2)-0.5*double(y%2))),2)+pow(0.8660254*(double(jj)-double(y)),2));
return dist;      
}


int main(int argc, char *argv[])
{

	MPI::Init(argc,argv);

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	int root = 0;

    
    const double EPSILON=0.0000001;
   ofstream outrun;
 
    outrun.open("clonalintlin");
 
    string statsfilename;  
    ofstream outstats;
    ifstream testoutstats;
	int filecounter=0;
	string tempstr;
   
    int32 seed = time(0)+rownd(double(rank)*double(time(0))/double(size));
    
    int leftedge, seedloc;
	bool survived;
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
    int numcolors;
    double basegrowth = -1;
    
    long int latticesizet;
    int nstep, nslices, ncounter;

	ostringstream tempstring;	



if (rank == 0) {


    cin >> latticesizex >> latticesizet >> selection >> mutf >> nslices >> numruns >> startseed >> numcolors;

    statsfilename = "s";
    tempstring << selection;	
    statsfilename += tempstring.str();
    tempstring.str("");
    tempstring.clear();
    tempstring << mutf;	
    statsfilename += "mf";
    statsfilename += tempstring.str();
    tempstring.str("");
    tempstring.clear();
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



    outstats << "# Lx: " << latticesizex << " Lt: " << latticesizet << " sel: " << selection << " mutf: " << mutf << endl;
    outstats << "# Nruns: " << numruns << " Nslices: " << nslices << " StartS: " << startseed << " Colors: " << numcolors << endl;


}

MPI::COMM_WORLD.Barrier();

	MPI::COMM_WORLD.Bcast(&latticesizex, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&latticesizet, 1, MPI::LONG, root);
	MPI::COMM_WORLD.Bcast(&selection, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&mutf, 1, MPI::DOUBLE, root);
	MPI::COMM_WORLD.Bcast(&numruns, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&nslices, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&startseed, 1, MPI::INT, root);
	MPI::COMM_WORLD.Bcast(&numcolors, 1, MPI::INT, root);


MPI::COMM_WORLD.Barrier();



    int centerx = rownd(latticesizex/2);
    int centery = rownd(latticesizet/2);
    numrunsd = double(numruns);
    latticesizexd=double(latticesizex);
    
    nstep=rownd(double(latticesizet)/double(nslices));
    if (nstep%2==1) {
                    nstep++;
                    nslices=rownd(double(latticesizet)/double(nstep));
                    }
                    

    TRanrotWGenerator rg(seed);
    
    vector< int > latteven(latticesizex);
    vector< int > lattodd(latticesizex);
 


    vector< vector< double > > sdencount(nslices+5);
    vector< vector< double > > scolorprobperc(nslices+5);
    vector< double > scolorprobperctemp(numcolors);
    vector< vector< double > > ssectorcount(nslices+5);
    vector< vector< double > > ssectorwidths(nslices+5);

	
    vector< vector< double > > dencount(nslices+5);
    vector< vector< double > > colorprobperc(nslices+5);
    vector< double > colorprobperctemp(numcolors);
    vector< vector< double > > sectorcount(nslices+5);
    vector< vector< double > > sectorwidths(nslices+5);
    vector< double > probperc(nslices+5);
    vector< double > sprobperc(nslices+5);

	vector< double > spread(nslices+5);
	vector< double > sspread(nslices+5);



    
    for(int i=0; i < nslices+5; i++) {

		spread[i]=0;
		sspread[i]=0;
		probperc[i]=0;
		sprobperc[i]=0;
            scolorprobperc[i].clear();
            sdencount[i].clear();
if (startseed>0){ 
            ssectorcount[i].clear();
            ssectorwidths[i].clear(); }

            colorprobperc[i].clear();
            dencount[i].clear();
if (startseed>0){ 
            sectorcount[i].clear();
            sectorwidths[i].clear(); }

            } 

    for (int i=0; i<nslices+5; i++)
    {


           for (int j=0; j<numcolors; j++) {
               scolorprobperc[i].push_back(0);
            if (startseed>0) {   ssectorwidths[i].push_back(0);
               			ssectorcount[i].push_back(0);   }
               sdencount[i].push_back(0);
                }    
                        

           for (int j=0; j<numcolors; j++) {
               colorprobperc[i].push_back(0);
            if (startseed>0) {   sectorwidths[i].push_back(0);
               			sectorcount[i].push_back(0);   }
               dencount[i].push_back(0);
                }    
                       


        }
    
// start the simulation runs
             
for (int run=1; run<=numruns; run++) {
 
// set up initial conditions 
  
  for (int i=0; i < latticesizex; i++) {latteven[i]=-1; lattodd[i]=-1;}
  
  if(startseed>0) {
          leftedge=rownd(latticesizex/2)-rownd(startseed/2);
          rightedge=rownd(latticesizex/2)-rownd(startseed/2)+startseed;
          for (int i1=0; i1 < leftedge; i1++) latteven[i1]=numcolors;        
          for (int i2=leftedge; i2 < rightedge; i2++) latteven[i2]=1;

		seedloc = leftedge;
         
          for (int i3=rightedge; i3 < latticesizex; i3++) latteven[i3]=numcolors;
                  }
  else {
   for (int i=0; i < latticesizex; i++) latteven[i]=1; 
       }
       
       ncounter=0;
        

                   it=0;
                   while ( (latteven[it]==latteven[(it+1)%latticesizex]) && it<=latticesizex) {
                         it++;
                         }             
                    it=(it+1)%latticesizex;           
                   sectorwidth=0; 
                   for(int j=0; j < numcolors; j++) colorprobperctemp[j]=0; 
			
				leftedge=latticesizex;
				rightedge=0;
				survived=0;
           
                   for(int i=0; i < latticesizex; i++) {            
                         dencount[ncounter][latteven[(i+it)%latticesizex]-1]+=1/numrunsd/latticesizexd;     
                         colorprobperctemp[latteven[(i+it)%latticesizex]-1]=1;  
                         if (startseed>0) {
                         if (latteven[(i+it)%latticesizex]!=latteven[(i+it+1)%latticesizex])
                         {
                                sectorwidths[ncounter][latteven[(i+it)%latticesizex]-1]+=sectorwidth/numrunsd;
                                sectorcount[ncounter][latteven[(i+it)%latticesizex]-1]+=1/numrunsd;
                                sectorwidth=0;                                                            
                                                        } else { sectorwidth=sectorwidth+1; } }
                         if (latteven[i]==1) {
					survived = 1;
					if (i<leftedge) leftedge=i;
					if (i>rightedge) rightedge=i;
					}                        
                               }


			if (survived) {
					spread[ncounter]+=rightedge-leftedge;
				 	probperc[ncounter]++;	
				}

                   for(int j=0; j < numcolors; j++) colorprobperc[ncounter][j]+=colorprobperctemp[j]/numrunsd;            
                   
                   if (sectorwidth==latticesizex && startseed>0) {
                           sectorcount[ncounter][latteven[0]-1]+=1/numrunsd; 
                           sectorwidths[ncounter][latteven[0]-1]+=sectorwidth/numrunsd; 
                                                  }            
                                                         
       ncounter++;
	   
  //  Begin simulation run:
   for (int t=1; t<=latticesizet; t++) {        
       
       // Odd steps:
              rannum1=rg.Random();
              rannum2=rg.Random();
       if (latteven[0]!=latteven[latticesizex-1]) {
                 if(rannum1 < (pow(1-selection,double(latteven[0]-1)))/(pow(1-selection,double(latteven[0]-1))+pow(1-selection,double(latteven[latticesizex-1]-1))))
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
                    if(rannum1 < (pow(1-selection,double(latteven[i]-1)))/(pow(1-selection,double(latteven[i]-1))+pow(1-selection,double(latteven[i-1]-1))))
                    {
                         lattodd[i]=latteven[i];      
                               } else { 
                                      lattodd[i]=latteven[i-1];
                                       }
                                                         }
           else {lattodd[i]=latteven[i];}
           if(rannum2<mutf && lattodd[i]<numcolors) lattodd[i]++;
           }

// If computing just one run, output the lattice configuration:

                    if (numruns==1) {
                           
                       outrun << latteven[latticesizex-1] << " "; 
                       for (int i=0; i < latticesizex-1; i++)      {     outrun << latteven[i] << " " << latteven[i] << " ";   }
                       outrun << latteven[latticesizex-1] << endl;
                       
                       outrun << latteven[latticesizex-1] << " "; 
                       for (int i=0; i < latticesizex-1; i++)      {     outrun << latteven[i] << " " << latteven[i] << " ";   }
                       outrun << latteven[latticesizex-1] << endl;
                       
                               
                       for (int i=0; i < latticesizex; i++)      {     outrun << lattodd[i] << " " << lattodd[i] << " ";    }
                           outrun << endl;  
                       
                       for (int i=0; i < latticesizex; i++)      {     outrun << lattodd[i] << " " << lattodd[i] << " ";    }
                           outrun << endl;           
                          
                             }

       t++;
       
       // Even steps:
      for (int i=0; i < latticesizex; i++) 
       {
              rannum1=rg.Random();
              rannum2=rg.Random();
                             
          if (lattodd[i]!=lattodd[(i+1)%latticesizex]) {
                if (rannum1 < (pow(1-selection,double(lattodd[i]-1)))/(pow(1-selection,double(lattodd[i]-1))+pow(1-selection,double(lattodd[i+1]-1))))
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
       
      // if time step is multiple of nstep, calculate statistics for the run: 
       if (t%nstep==0 && t<=latticesizet) {
 
                   it=0;
                   while ( (latteven[it]==latteven[(it+1)%latticesizex]) && it<=latticesizex) {
                         it++;
                         }             
                    it=(it+1)%latticesizex;           
                   sectorwidth=0; 
                   for(int j=0; j < numcolors; j++) colorprobperctemp[j]=0;  


			survived=0;
			leftedge=latticesizex-1;
			rightedge=0;

          
                   for(int i=0; i < latticesizex; i++) {            
                         dencount[ncounter][latteven[(i+it)%latticesizex]-1]+=1/latticesizexd/numrunsd; 
                        colorprobperctemp[latteven[(i+it)%latticesizex]-1]=1;  
                         
                         if (startseed>0) {
                         if (latteven[(i+it)%latticesizex]!=latteven[(i+it+1)%latticesizex])
                         {
                                sectorwidths[ncounter][latteven[(i+it)%latticesizex]-1]+=sectorwidth/numrunsd;
                                sectorcount[ncounter][latteven[(i+it)%latticesizex]-1]+=1/numrunsd;
                                sectorwidth=0;                                                            
                                                        } else { sectorwidth=sectorwidth+1; }
                               
			    if (latteven[i]==1) {
					survived=1;
					if (i<leftedge) leftedge=i;
					if (i>rightedge) rightedge=i;
					}
				          }               
                               }
		      
			if (survived) {
				spread[ncounter]+=rightedge-leftedge;
				probperc[ncounter]++;
					}	
                   for(int j=0; j < numcolors; j++) colorprobperc[ncounter][j]+=colorprobperctemp[j]/numrunsd;            
                   if (sectorwidth==latticesizex && startseed>0) {
                           sectorcount[ncounter][latteven[0]-1]+=1/numrunsd; 
                           sectorwidths[ncounter][latteven[0]-1]+=sectorwidth/numrunsd; 
                                                  }           
              ncounter++;   
                       } 
                                      } // end of time-step loop
                                                   } // end of simulation run loop
	for (int i=0; i < nslices; i++) {
		for (int j=0; j < numcolors; j++) {
			dencount[i][j] = dencount[i][j]/double(size);
		        colorprobperc[i][j]=colorprobperc[i][j]/double(size);
		if (startseed>0) {
			 sectorcount[i][j] = sectorcount[i][j]/double(size);			
			 sectorwidths[i][j] = sectorwidths[i][j]/double(size);
			}
		}
		if (probperc[i]>EPSILON) {
			spread[i]=spread[i]/probperc[i];	
				}
	}

if (numruns==1) outrun.close();


// compile the statistics from all processes:

MPI::COMM_WORLD.Barrier();

	MPI::COMM_WORLD.Reduce(&spread.front(),&sspread.front() , spread.size() , MPI::DOUBLE, MPI::SUM, root);
	MPI::COMM_WORLD.Reduce(&probperc.front(),&sprobperc.front() , probperc.size() , MPI::DOUBLE, MPI::SUM, root);
	

	for(int i=0; i < nslices; i++)
{
	MPI::COMM_WORLD.Reduce(&dencount[i].front(),&sdencount[i].front() , numcolors , MPI::DOUBLE, MPI::SUM, root);
	MPI::COMM_WORLD.Reduce(&colorprobperc[i].front(),&scolorprobperc[i].front() , numcolors , MPI::DOUBLE, MPI::SUM, root);
	if (startseed>0) {
	MPI::COMM_WORLD.Reduce(&sectorcount[i].front(),&ssectorcount[i].front() ,sectorcount[i].size(), MPI::DOUBLE, MPI::SUM, root);
	MPI::COMM_WORLD.Reduce(&sectorwidths[i].front(),&ssectorwidths[i].front() ,sectorwidths[i].size(), MPI::DOUBLE, MPI::SUM, root);
                         }
}

// output the compiled statistics:

MPI::COMM_WORLD.Barrier();

if (rank==0) {  
             for (int i=0; i < nslices; i++)
             {
                   outstats << double(i*nstep) << " " << sspread[i]/double(size) << " " << sprobperc[i];
     for (int j=0; j < 1; j++) {
          if (startseed==0) {          
				outstats << " " <<  sdencount[i][j]; 
					   }
                          if (startseed>0) {
  outstats << " " << scolorprobperc[i][j]  << " " << ssectorcount[i][j] << " " << ssectorwidths[i][j];
						 }
                     }
                   outstats << endl;   
             }              
		outstats.close();
             }
				  
MPI::COMM_WORLD.Barrier();             
MPI::Finalize();


  
   
    return( 0 );
}

