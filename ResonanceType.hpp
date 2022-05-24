#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP
#include "ParticleType.hpp"


    class ResonanceType: public ParticleType{
        public:
        ResonanceType(const char* Name, const double Mass, const int Charge, const double Width);
double GetWidth() const;
void PrintParticle() const;

        private:
      const double fWidth;
    };

#endif