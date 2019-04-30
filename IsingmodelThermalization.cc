/* 
Ising Model for Computational Physics
*/

#include <iostream> //cout
#include <fstream>  //instream, ofstream
#include <cstdlib>  //rand, srand
#include <ctime>    //time_t, time
#include <cmath>    //abs, exp, pow
using namespace std;

const int SIZE=20;
const int NSTEPS=5;
const int INTERACTION=1;
const int EXTERNAL_FIELD=0;
const double START_TEMPERATURE=2;
const double STOP_TEMPERATURE=2;
const double TEMPERATURE_STEP_SIZE=0.1;
const bool USE_WOLFF=false;

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
  ofstream output("Lattice.txt");
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      output <<(grid[i][j]+1)/2<<" ";
    }
    output <<"\n";
  }
  output.close();
}//Lattice::print

double Lattice::computeHamiltonian(){
  double H=0.0;
  for(int i=0;i<SIZE;i++){
    for(int j=0;j<SIZE;j++){
      H=H-J*(grid[i][j]*grid[(i+1)%SIZE][j]+grid[i][j]*grid[i][(j+1)%SIZE])-B*grid[i][j];
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
  return(2*J*(grid[i][j]*grid[(i+1)%SIZE][j]+grid[i][j]*grid[i][(j+1)%SIZE]+
    grid[i][j]*grid[(SIZE+i-1)%SIZE][j]+grid[i][j]*grid[i][(SIZE+j-1)%SIZE])+2*B*grid[i][j]);
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
  hamiltonian+=difference;
  magnetisation-=2*grid[i][j];
  grid[i][j]=-grid[i][j];
  if((grid[(i+1)%SIZE][j]==-grid[i][j]) and (rand()/RAND_MAX < (1-exp(-2*beta*J)))){
    doWolffRecursion((i+1)%SIZE, j);
  }
  if((grid[i][(j+1)%SIZE]==-grid[i][j]) and (rand()/RAND_MAX < (1-exp(-2*beta*J)))){
    doWolffRecursion(i, (j+1)%SIZE);
  }
  if((grid[(SIZE+i-1)%SIZE][j]==-grid[i][j]) and (rand()/RAND_MAX < (1-exp(-2*beta*J)))){
    doWolffRecursion((SIZE+i-1)%SIZE, j);
  }
  if((grid[i][(SIZE+j-1)%SIZE]==-grid[i][j]) and (rand()/RAND_MAX < (1-exp(-2*beta*J)))){
    doWolffRecursion(i, (SIZE+j-1)%SIZE);    
  }
}//Lattice::doWolfRecursion

void Lattice::doWolffTimestep(){
  int i, j;
  i=rand()%SIZE;
  j=rand()%SIZE;
  doWolffRecursion(i, j);
}//Lattice::doWolffTimestep

int main(){
  srand(time(NULL));
  int stepSize=pow(SIZE,2);
  double nTemperatureSteps=(STOP_TEMPERATURE-START_TEMPERATURE)/TEMPERATURE_STEP_SIZE;
  int nSwitchPositiveHamiltonian=0, nSwitchNegativeHamiltonian=0, 
      nSwitchPositiveMagnetisation=0, nSwitchNegativeMagnetisation=0, i=0;
  Lattice model, positiveModel, negativeModel;
  
  ofstream output("thermalizationdata.txt");
  output <<"Temperature, Iteration, " 
         <<"ModelEnergy, ModelMagnetisation, " 
         <<"PositiveModelEnergy, PositiveModelMagnetisation, "
         <<"NegativeModelEnergy, NegativeModelMagnetisation\n";
  for(int j=0;j<=nTemperatureSteps;j++){
    for(int k=1;k<=NSTEPS;k++){
      positiveModel.configurePositive();
      negativeModel.configureNegative();
      positiveModel.setTemperature(model.getTemperature());
      negativeModel.setTemperature(model.getTemperature());
      model.configureRandomly();
      while(nSwitchPositiveHamiltonian<500 or nSwitchNegativeHamiltonian<500 or 
            nSwitchPositiveMagnetisation<500 or nSwitchNegativeMagnetisation<500){
        if(pow(-1.0,nSwitchPositiveHamiltonian)*(model.getHamiltonian()-positiveModel.getHamiltonian())>=0){
          nSwitchPositiveHamiltonian++;
        }
        if(pow(-1.0,nSwitchNegativeHamiltonian)*(model.getHamiltonian()-negativeModel.getHamiltonian())>=0){
          nSwitchNegativeHamiltonian++;
        }
        if(pow(-1.0,nSwitchPositiveMagnetisation)*(abs(model.getMagnetisation())-abs(positiveModel.getMagnetisation()))>=0){
          nSwitchPositiveMagnetisation++;
        }
        if(pow(-1.0,nSwitchNegativeHamiltonian)*(abs(model.getMagnetisation())-abs(negativeModel.getMagnetisation()))>=0){
          nSwitchNegativeMagnetisation++;
        }
        if(i%stepSize==0){
          output <<model.getTemperature()<<", "<<k<<", "
                 <<model.getHamiltonian()<<", "<<model.getMagnetisation()<<", "
                 <<positiveModel.getHamiltonian()<<", "<<positiveModel.getMagnetisation()<<", "
                 <<negativeModel.getHamiltonian()<<", "<<negativeModel.getMagnetisation()<<"\n";
        }
        i++;
        model.doTimestep();
        positiveModel.doTimestep();
        negativeModel.doTimestep();
      }
      nSwitchPositiveHamiltonian=0;
      nSwitchNegativeHamiltonian=0;
      nSwitchPositiveMagnetisation=0;
      nSwitchNegativeMagnetisation=0;
      i=0;
    }
    model.setTemperature(model.getTemperature()+TEMPERATURE_STEP_SIZE);
  }
  output.close();
  return 0;
}//main

