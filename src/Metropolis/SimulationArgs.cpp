#include "SimulationArgs.h"

/** Convert a string strategy to a SimulationStrategy type */
SimulationStrategy Strategy::fromString(std::string type) {
  if (type == "brute" || type == "brute-force") {
    return Strategy::BruteForce;
  } else if (type == "prox" || type == "proximity-matrix") {
    return Strategy::ProximityMatrix;
  } else if (type == "nlc" || type == "neighbor-linked-cell") {
    return Strategy::NLC;
  } else if (type == "verlet" || type == "verlet-list") {
    return Strategy::Verlet;
  } else {
    return Strategy::Unknown;
  }
}
