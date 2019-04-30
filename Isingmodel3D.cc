/* 
Ising Model for Computational Physics
*/

#include <iostream> //cout
#include <fstream>  //instream, ofstream
#include <cstdlib>  //rand, srand
#include <ctime>    //time_t, time
#include <cmath>    //abs, exp, pow
using namespace std;

const int SIZE=40;
const int NSTEPS=20;
const int INTERACTION=1;
const int EXTERNAL_FIELD=0;
const double START_TEMPERATURE=3.0;
const double STOP_TEMPERATURE=6.0;
const double TEMPERATURE_STEP_SIZE=0.05;
const bool USE_WOLFF=true;

class Lattice{
  public:
    Lattice();
    void configurePositive();
    void configureNegative();
    void configureRandomly();
    double computeHamiltonian();
    int computeMagnetisation();
    double getHamiltonian();
    int getMagnetisation();
    double getTemperature();
    void setTemperature(double temperature);
    double computeHamiltonianDifference(int i, int j, int k);
    void doTimestep();
    void doWolffRecursion(int i, int j, int k);
    void doWolffTimestep();
    void thermalize();
  private:
    int grid[SIZE][SIZE][SIZE];
    int J;
    int B;
    double beta;
    double hamiltonian;
    double magnetisation;
    double probabilityWolff;
};//Lattice

Lattice::Lattice(){
  J=INTERACTION;
  B=EXTERNAL_FIELD;
  beta=1/START_TEMPERATURE;
  probabilityWolff=1-exp(-2*beta*J);
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      for(int k=0;k<SIZE;k++){
        grid[i][j][k]=1;
      }
    }
  }
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::constructor

void Lattice::configurePositive(){
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      for(int k=0;k<SIZE;k++){
        grid[i][j][k]=1;
      }
    }
  }
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::configurePositive

void Lattice::configureNegative(){
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      for(int k=0;k<SIZE;k++){
        grid[i][j][k]=-1;
      }
    }
  }
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::configureNegative

void Lattice::configureRandomly(){
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      for(int k=0;k<SIZE;k++){
        grid[i][j][k]=2*(rand()%2)-1;
      }
    }
  }
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::configureRandomly

double Lattice::computeHamiltonian(){
  double H=0;
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      for(int k=0;k<SIZE;k++){
        H=H-J*grid[i][j][k]*(grid[(i+1)%SIZE][j][k]+grid[i][(j+1)%SIZE][k]+grid[i][j][(k+1)%SIZE])-B*grid[i][j][k];
      }
    }
  }
  return H;
}//Lattice::computeHamiltonian

int Lattice::computeMagnetisation(){
  int M=0;
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      for(int k=0;k<SIZE;k++){
        M=M+grid[i][j][k];
      }
    }
  }
  return M;
}//Lattice::computeMagnetisation

double Lattice::getHamiltonian(){
  return hamiltonian;
}//Lattice::getHamiltonian

int Lattice::getMagnetisation(){
  return magnetisation;
}//Lattice::getMagnetisation

double Lattice::getTemperature(){
  return 1/beta;
}//Lattice::getTemperature

void Lattice::setTemperature(double temperature){
  beta=1/temperature;
  probabilityWolff=1-exp(-2*beta*J);
}//Lattice::setTemperature

double Lattice::computeHamiltonianDifference(int i, int j, int k){
  return(2*J*grid[i][j][k]*(grid[(i+1)%SIZE][j][k]+grid[i][(j+1)%SIZE][k]+grid[i][j][(k+1)%SIZE]+
         grid[(SIZE+i-1)%SIZE][j][k]+grid[i][(SIZE+j-1)%SIZE][k]+grid[i][j][(SIZE+k-1)%SIZE])+2*B*grid[i][j][k]);
}//Lattice::computeHamiltonianDifference

void Lattice::doTimestep(){
  int i, j, k;
  i=rand()%SIZE;
  j=rand()%SIZE;
  k=rand()%SIZE;
  double randomUniform=rand();
  double difference=computeHamiltonianDifference(i, j, k);
  if((difference<0) or (randomUniform/RAND_MAX < exp(-beta*(difference)))){
    hamiltonian+=difference;
    magnetisation-=2*grid[i][j][k];
    grid[i][j][k]=-grid[i][j][k];
  }
}//Lattice::doTimestep

void Lattice::doWolffRecursion(int i, int j, int k){
  double difference=computeHamiltonianDifference(i, j, k);
  double random1, random2, random3, random4, random5, random6;
  hamiltonian+=difference;
  magnetisation-=2*grid[i][j][k];
  grid[i][j][k]=-grid[i][j][k];
  random1=rand();
  random2=rand();
  random3=rand();
  random4=rand();
  random5=rand();
  random6=rand();
  if((grid[(i+1)%SIZE][j][k]==-grid[i][j][k]) and (random1/RAND_MAX<probabilityWolff)){
    doWolffRecursion((i+1)%SIZE, j, k);
  }
  if((grid[i][(j+1)%SIZE][k]==-grid[i][j][k]) and (random2/RAND_MAX<probabilityWolff)){
    doWolffRecursion(i, (j+1)%SIZE, k);
  }
  if((grid[i][j][(k+1)%SIZE]==-grid[i][j][k]) and (random3/RAND_MAX<probabilityWolff)){
    doWolffRecursion(i, j, (k+1)%SIZE);
  }
  if((grid[(SIZE+i-1)%SIZE][j][k]==-grid[i][j][k]) and (random4/RAND_MAX<probabilityWolff)){
    doWolffRecursion((SIZE+i-1)%SIZE, j, k);    
  }
  if((grid[i][(SIZE+j-1)%SIZE][k]==-grid[i][j][k]) and (random5/RAND_MAX<probabilityWolff)){
    doWolffRecursion(i, (SIZE+j-1)%SIZE, k);    
  }
  if((grid[i][j][(SIZE+k-1)%SIZE]==-grid[i][j][k]) and (random6/RAND_MAX<probabilityWolff)){
    doWolffRecursion(i, j, (SIZE+k-1)%SIZE);    
  }
}//Lattice::doWolfRecursion

void Lattice::doWolffTimestep(){
  int i, j, k;
  i=rand()%SIZE;
  j=rand()%SIZE;
  k=rand()%SIZE;
  doWolffRecursion(i, j, k);
}//Lattice::doWolffTimestep

void Lattice::thermalize(){
  int nSwitchPositiveHamiltonian=0, nSwitchNegativeHamiltonian=0, 
      nSwitchPositiveMagnetisation=0, nSwitchNegativeMagnetisation=0, i=0;
  Lattice positiveModel, negativeModel;
  negativeModel.configureNegative();
  positiveModel.setTemperature(getTemperature());
  negativeModel.setTemperature(getTemperature());
  configureRandomly();
   
  //The model is thermalized if all three models reach similar energy levels.	 
  while(nSwitchPositiveHamiltonian<500 or nSwitchNegativeHamiltonian<500 or 
        nSwitchPositiveMagnetisation<500 or nSwitchNegativeMagnetisation<500){
    i++;
    if(USE_WOLFF){
      doWolffTimestep();
      positiveModel.doWolffTimestep();
      negativeModel.doWolffTimestep();
    }
    else{
      doTimestep();
      positiveModel.doTimestep();
      negativeModel.doTimestep();
    }
    if(pow(-1.0,nSwitchPositiveHamiltonian)*(getHamiltonian()-positiveModel.getHamiltonian())>=0){
      nSwitchPositiveHamiltonian++;
    }
    if(pow(-1.0,nSwitchNegativeHamiltonian)*(getHamiltonian()-negativeModel.getHamiltonian())>=0){
      nSwitchNegativeHamiltonian++;
    }
    if(pow(-1.0,nSwitchPositiveMagnetisation)*(abs(getMagnetisation())-abs(positiveModel.getMagnetisation()))>=0){
      nSwitchPositiveMagnetisation++;
    }
    if(pow(-1.0,nSwitchNegativeHamiltonian)*(abs(getMagnetisation())-abs(negativeModel.getMagnetisation()))>=0){
      nSwitchNegativeMagnetisation++;
    }
  }
}//Lattice::thermalize

int main(){
  srand(time(NULL));
  int stepSize=pow(SIZE,3);
  int nTemperatureSteps=(STOP_TEMPERATURE-START_TEMPERATURE)/TEMPERATURE_STEP_SIZE;
  Lattice model;
  model.configureRandomly();
  model.thermalize();
  
  ofstream output("data.txt");
  output <<"Temperature, Energy, Magnetisation\n";
  for(int k=1;k<=nTemperatureSteps;k++){
    for(int i=1;i<=NSTEPS;i++){
      output <<model.getTemperature()<<", "
             <<model.getHamiltonian()<<", "<<model.getMagnetisation()<<"\n";
      if(USE_WOLFF){
        model.doWolffTimestep();
      }
      else{
        for(int j=1;j<=stepSize;j++){
          model.doTimestep();
        }
      }
    }
    model.setTemperature(model.getTemperature()+TEMPERATURE_STEP_SIZE);
    model.thermalize();
  }
  
  output.close();
  return 0;
}//main

