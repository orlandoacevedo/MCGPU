#include "SimulationArgs.h"

/** Convert a string strategy to a SimulationStrategy type */
SimulationStrategy Strategy::fromString(std::string type) {
  if (type == "brute" || type == "brute-force") {
    return Strategy::BruteForce;
  } else if (type == "prox" || type == "proximity-matrix") {
    return Strategy::ProximityMatrix;
  } else {
    return Strategy::Unknown;
  }
}
