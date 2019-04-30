/* 
Ising Model for Computational Physics
*/

#include <iostream> //cout
#include <fstream>  //instream, ofstream
#include <cstdlib>  //rand, srand
#include <ctime>    //time_t, time
#include <cmath>    //abs, exp, pow
using namespace std;

const int SIZE=10;
const int NSTEPS=20000;
const int INTERACTION=1;
const int EXTERNAL_FIELD=0;
const double START_TEMPERATURE=0.0;
const double STOP_TEMPERATURE=5.0;
const double TEMPERATURE_STEP_SIZE=0.01;
const bool USE_WOLFF=true;

class Lattice{
  public:
    Lattice();
    void configurePositive();
    void configureNegative();
    void configureRandomly();
    void configureManually();
    void print();
    double computeHamiltonian();
    int computeMagnetisation();
    double getHamiltonian();
    int getMagnetisation();
    double getTemperature();
    void setTemperature(double temperature);
    double computeHamiltonianDifference(int i, int j);
    void doTimestep();
    void doWolffRecursion(int i, int j);
    void doWolffTimestep();
    void thermalize();
//  private:
    int grid[SIZE][SIZE];
    int J;
    int B;
    double beta;
    double hamiltonian;
    double magnetisation;
};//Lattice

Lattice::Lattice(){
  J=INTERACTION;
  B=EXTERNAL_FIELD;
  beta=1/START_TEMPERATURE;
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      grid[i][j]=1;
    }
  }
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::constructor

void Lattice::configurePositive(){
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      grid[i][j]=1;
    }
  }
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::configurePositive

void Lattice::configureNegative(){
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      grid[i][j]=-1;
    }
  }
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::configureNegative

void Lattice::configureRandomly(){
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      grid[i][j]=2*(rand()%2)-1;
    }
  }
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::configureRandomly

void Lattice::configureManually(){
  char input_char;
  ifstream input("startingconfiguration.txt");
  if(input.is_open()){
    for(int i=0;i<SIZE;i++){
      for(int j=0;j<SIZE;j++){
        input_char=input.get();
        while(!(input_char=='0' or input_char=='1') and input.good()){
          input_char=input.get();
        }
        grid[i][j]=2*(input_char-'0')-1;
      }
    }
  }
  else{
    cout <<"Could not open file./n";
  }
  input.close();
  hamiltonian=computeHamiltonian();
  magnetisation=computeMagnetisation();
}//Lattice::configureManually

void Lattice::print(){
  ofstream output("lattice.txt");
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      output <<(grid[i][j]+1)/2<<" ";
    }
    output <<"\n";
  }
  output.close();
}//Lattice::print

double Lattice::computeHamiltonian(){
  double H=0;
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      H=H-J*grid[i][j]*(grid[(i+1)%SIZE][j]+grid[i][(j+1)%SIZE])-B*grid[i][j];
    }
  }
  return H;
}//Lattice::computeHamiltonian

int Lattice::computeMagnetisation(){
  int M=0;
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      M=M+grid[i][j];
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
}//Lattice::setTemperature

double Lattice::computeHamiltonianDifference(int i, int j){
  return(2*J*grid[i][j]*(grid[(i+1)%SIZE][j]+grid[i][(j+1)%SIZE]+
         grid[(SIZE+i-1)%SIZE][j]+grid[i][(SIZE+j-1)%SIZE])+2*B*grid[i][j]);
}//Lattice::computeHamiltonianDifference

void Lattice::doTimestep(){
  int i, j;
  i=rand()%SIZE;
  j=rand()%SIZE;
  double randomUniform=rand();
  double difference=computeHamiltonianDifference(i, j);
  if((difference<0) or (randomUniform/RAND_MAX < exp(-beta*(difference)))){
    hamiltonian+=difference;
    magnetisation-=2*grid[i][j];
    grid[i][j]=-grid[i][j];
  }
}//Lattice::doTimestep

void Lattice::doWolffRecursion(int i, int j){
  double difference=computeHamiltonianDifference(i, j);
  double random1, random2, random3, random4;
  hamiltonian+=difference;
  magnetisation-=2*grid[i][j];
  grid[i][j]=-grid[i][j];
  random1=rand();
  random2=rand();
  random3=rand();
  random4=rand();
  if((grid[(i+1)%SIZE][j]==-grid[i][j]) and (random1/RAND_MAX<(1-exp(-2*beta*J)))){
    doWolffRecursion((i+1)%SIZE, j);
  }
  if((grid[i][(j+1)%SIZE]==-grid[i][j]) and (random2/RAND_MAX<(1-exp(-2*beta*J)))){
    doWolffRecursion(i, (j+1)%SIZE);
  }
  if((grid[(SIZE+i-1)%SIZE][j]==-grid[i][j]) and (random3/RAND_MAX<(1-exp(-2*beta*J)))){
    doWolffRecursion((SIZE+i-1)%SIZE, j);
  }
  if((grid[i][(SIZE+j-1)%SIZE]==-grid[i][j]) and (random4/RAND_MAX<(1-exp(-2*beta*J)))){
    doWolffRecursion(i, (SIZE+j-1)%SIZE);    
  }
}//Lattice::doWolfRecursion

void Lattice::doWolffTimestep(){
  int i, j;
  i=rand()%SIZE;
  j=rand()%SIZE;
  doWolffRecursion(i, j);
}//Lattice::doWolffTimestep

void Lattice::thermalize(){
  int nSwitchPositiveHamiltonian=0, nSwitchNegativeHamiltonian=0, 
      nSwitchPositiveMagnetisation=0, nSwitchNegativeMagnetisation=0, i=0;
  Lattice positiveModel, negativeModel;
  negativeModel.configureNegative();
  positiveModel.setTemperature(getTemperature());
  negativeModel.setTemperature(getTemperature());
  configureRandomly();
   
  //The model is thermalized if all three models reach similar energy and magnetisation levels.	 
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
  int stepSize=pow(SIZE,2);
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

