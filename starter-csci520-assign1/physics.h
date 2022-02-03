/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);

struct point structuralForce(struct world * jello, int i, int j, int k);
struct point sheerForce(struct world * jello, int i, int j, int k);
struct point bendForce(struct world * jello, int i, int j, int k);
struct point fieldForce(struct world * jello, int i, int j, int k);
struct point wallColForce(struct world * jello, int i, int j, int k);

struct point hookForce(struct world * jello, int i1, int j1, int k1, int i2, int j2, int k2);
struct point dampingForce(struct world * jello, int i1, int j1, int k1, int i2, int j2, int k2);
struct point collisionForce(struct world * jello, int i, int j, int k, double a, double b, double c, double d);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);

#endif

