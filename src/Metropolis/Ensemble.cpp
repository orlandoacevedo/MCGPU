#include "Ensemble.h"

void NVTEnsemble::update() {}

bool NVTEnsemble::accept(Real deltaE) {
  return SimCalcs::acceptMove(deltaE);
}

void NPTEnsemble::update() {
  step++;
  if (step % volumeInterval != 0) return;

  Real* bSize = GPUCopy::sizePtr();
  oldVolume =  1.0;
  for (int i = 0; i < NUM_DIMENSIONS; i++) {
    oldVolume *= bSize[i];
  }
  SimCalcs::resizeBox(1.1);
}

bool NPTEnsemble::accept(Real deltaE) {
  return SimCalcs::acceptVolumeMove(deltaE, oldVolume, pressure);
}

NPTEnsemble::NPTEnsemble(Real pressure_in, int volumeInterval_in) {
  pressure = pressure_in;
  volumeInterval = volumeInterval_in;
}

