#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP

 
    class ParticleType{
        public:
        ParticleType(const char* Name,const double Mass,const int Charge);
        char const* GetName() const;
        double GetMass() const;
        int GetCharge() const;
        void PrintParticleType() const;
        virtual double GetWidth();
      
        private:
         const char* fName;
         const double fMass;
         const int fCharge;
    };


#endif