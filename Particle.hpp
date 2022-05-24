
#ifndef PARTICLE_HPP
#define PARTICLE_HPP

class ParticleType;

class Particle {
 public:
  
  explicit Particle(const char* Name, double Px = 0, double Py = 0, double Pz = 0);
  static const int fMaxNumParticleType = 10;

  static void PrintParticleTypes();

Particle();

  void Print() const;
  double GetMass() const;
  double GetEnergy() const;
  int GetCharge() const;

  double Momentum2() const;
  double Momentum() const;

  

  double InvariantMass(Particle& Particle) const;

  static bool AddParticleType(const char* Name, double Mass, int Charge, double Width = 0);

  int GetIndex() const;
  double GetPx() const;
  double GetPy() const;
  double GetPz() const;

  bool SetIndex(int Index);
  
  bool SetIndex(const char* Name);

 
  void SetP(double Px, double Py, double Pz);

  int Decay2body(Particle &dau1,Particle &dau2) const;
  void Boost(double bx, double by, double bz);
  

 private:
  // Array containing all active particle types
  static ParticleType* fParticleTypes[fMaxNumParticleType];

  // Number of types in ParticleTypes_ array
  static int fNParticleType;

  // Index of this particle in ParticleTypes_ array
  int fIndex;

  // Momentum components
  double fPx{0};
  double fPy{0};
  double fPz{0};

  // Methods ///////////////////////////////////////////////////////////////////
  // Find a particle type in fParticleTypes array by name, return its index.
  // Returns -1 if not found.
  static int FindParticle(const char* Name);

  
};

#endif