#ifndef METROPOLIS_ENSEMBLE_H
#define METROPOLIS_ENSEMBLE_H

#include <stdbool.h>

#include "SimulationStep.h"
#include "GPUCopy.h"

class Ensemble {

public:
  virtual bool accept(Real deltaE) = 0;

  virtual void update() = 0;
};

class NVTEnsemble: public Ensemble {

public:
  void update();

  bool accept(Real deltaE);
};

class NPTEnsemble: public Ensemble {

public:
  Real oldVolume;
  Real pressure;
  int volumeInterval;
  int step;

  NPTEnsemble(Real pressure_in, int volumeInterval_in);

  void update();

  bool accept(Real deltaE);
};

#endif
