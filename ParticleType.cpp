#include "ParticleType.hpp"
#include <iostream>

    ParticleType::ParticleType(const char* Name, const double Mass, const int Charge):
    fName{Name}, fMass{Mass}, fCharge{Charge}
    {}

    const char* ParticleType::GetName() const{
        return fName;
    }
    double ParticleType::GetMass() const{
        return fMass;
    }
    int ParticleType::GetCharge() const{
        return fCharge;
    }
    void ParticleType::PrintParticleType() const{
        std::cout<<"Particle Name:" << fName <<'\n';
        std::cout<<"Particle Mass:" << fMass <<'\n';
        std::cout<<"Particle Charge:"<< fCharge<<'\n';
    }
    double ParticleType::GetWidth() {
        return 0;
    }





