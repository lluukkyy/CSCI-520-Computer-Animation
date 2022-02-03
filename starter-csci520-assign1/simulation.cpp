#include "jello.h"
#include "simulation.h"
#include "physics.h"

void step(struct world * jello){
    Euler(jello);
}