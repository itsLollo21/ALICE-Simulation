#include "ResonanceType.hpp"
#include <iostream>

    ResonanceType::ResonanceType(const char* Name, const double Mass, const int Charge, const double Width):
    ParticleType(Name, Mass, Charge),
    fWidth{Width}
    {}
    double ResonanceType::GetWidth() const{
        return fWidth;
    }
    void ResonanceType::PrintParticle() const{
        ParticleType::PrintParticleType();
        std::cout<<"Resonance Width:"<< fWidth <<'\n';
    }



