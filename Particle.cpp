#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>

  
    
    

     ParticleType* Particle::fParticleTypes[fMaxNumParticleType];
     int Particle::fNParticleType{0};
     
Particle::Particle() : fIndex{}, fPx{}, fPy{}, fPz{} {};

    Particle::Particle(const char* Name, double Px, double Py, double Pz):
    fPx{Px}, fPy{Py}, fPz{Pz}
    {
        const int Index = FindParticle(Name);
        if (Index==-1){
            std::cout<<"Error, no particle called" << Name << "found."<<'\n';
        }
        fIndex = Index;
    }
    void Particle::Print() const{
        std::cout<<"Particle type Index: " << fIndex <<'\n';
        std::cout<<"Particle name: " << fParticleTypes[fIndex]->GetName() <<'\n';
        std::cout<<"Particle momentum: P=("<<fPx <<","<<fPy <<","<<fPz<<")"<<'\n';
    }
    double Particle::GetMass() const{
        return fParticleTypes[fIndex]->GetMass();
    }
    double Particle::Momentum() const{
        return sqrt(fPx*fPx+ fPy*fPy + fPz*fPz);
    }
    double Particle::Momentum2() const{
        return fPx*fPx + fPy*fPy + fPz*fPz;
    }
    
    double Particle::GetEnergy() const{
        return sqrt((pow(GetMass(),2)+Momentum2()));
    }

    int Particle::GetCharge() const { return fParticleTypes[fIndex]->GetCharge(); }

    double Particle::InvariantMass(Particle& Particle) const {
        double TotalEnergySum = GetEnergy()+Particle.GetEnergy();
        double MomentumSum = Momentum()+Particle.Momentum();
        return sqrt(pow(TotalEnergySum, 2)- pow(MomentumSum, 2));
    }
    bool Particle::AddParticleType(const char* Name, double Mass, int Charge, double Width){
        if(FindParticle(Name) != -1){
            std::cout<< "Particle name already in the Particles Array." <<'\n';
            return false;
        }
        if (fNParticleType == 10){
            std::cout <<"The Particles Array is full." <<'\n';
            return false;
        }
        if (Width == 0){
        fParticleTypes[fNParticleType ++] = new ParticleType(Name, Mass, Charge);
        return true;
        }
        fParticleTypes[fNParticleType ++] = new ResonanceType(Name, Mass, Charge, Width);
        return true;
    }
       
       int Particle::FindParticle(const char* Name){
           for(int i=0; i< fNParticleType; ++i){
              if(std::strcmp(fParticleTypes[i]->GetName(), Name)==0){
                  return i;
              }
           }
           return -1;
       }

       void Particle::PrintParticleTypes(){
           for(int i=0; i< fNParticleType; ++i){
               fParticleTypes[i]->PrintParticleType();
           }

       }
       int Particle::GetIndex() const{
           return fIndex;
       }
       double Particle::GetPx() const {
           return fPx;
       }
       double Particle::GetPy() const {
           return fPy;
       }
       double Particle::GetPz() const {
           return fPz;
       }


        bool Particle::SetIndex(const char* Name){
            return SetIndex(FindParticle(Name));
        }

        void Particle::SetP(double Px, double Py, double Pz){
            fPx=Px;
            fPy=Py;
            fPz=Pz;
        }       
       
       
       int Particle::Decay2body(Particle &dau1,Particle &dau2) const {
  if(GetMass() == 0.0){
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }
  
  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if(fIndex > -1){ // add width effect

    // gaussian random numbers

    float x1, x2, w, y1, y2;
    
    double invnum = 1./RAND_MAX;
    do {
      x1 = 2.0 * rand()*invnum - 1.0;
      x2 = 2.0 * rand()*invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += fParticleTypes[fIndex]->GetWidth() * y1;

  }

  if(massMot < massDau1 + massDau2){
    printf("Decayment cannot be preformed because mass is too low in this channel\n");
    return 2;
  }
  
  double pout = sqrt((massMot*massMot - (massDau1+massDau2)*(massDau1+massDau2))*(massMot*massMot - (massDau1-massDau2)*(massDau1-massDau2)))/massMot*0.5;

  double norm = 2*M_PI/RAND_MAX;

  double phi = rand()*norm;
  double theta = rand()*norm*0.5 - M_PI/2.;
  dau1.SetP(pout*sin(theta)*cos(phi),pout*sin(theta)*sin(phi),pout*cos(theta));
  dau2.SetP(-pout*sin(theta)*cos(phi),-pout*sin(theta)*sin(phi),-pout*cos(theta));

  double energy = sqrt(fPx*fPx + fPy*fPy + fPz*fPz + massMot*massMot);

  double bx = fPx/energy;
  double by = fPy/energy;
  double bz = fPz/energy;

  dau1.Boost(bx,by,bz);
  dau2.Boost(bx,by,bz);

  return 0;
}
void Particle::Boost(double bx, double by, double bz)
{

  const double energy = GetEnergy();

  //Boost this Lorentz vector
  double b2 = bx*bx + by*by + bz*bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx*fPx + by*fPy + bz*fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0)/b2 : 0.0;

  fPx += gamma2*bp*bx + gamma*bx*energy;
  fPy += gamma2*bp*by + gamma*by*energy;
  fPz += gamma2*bp*bz + gamma*bz*energy;
}


       
    








