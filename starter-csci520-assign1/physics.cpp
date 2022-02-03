/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <algorithm>

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  /* for you to implement ... */
  // struct point x = {1,1,1};
  // struct point nx = x.normalize();
  // printf("\n\n\n??? %f (%f %f %f) \n\n\n", x.norm(), nx.x,nx.y,nx.z);
  struct point stForce, shForce, bdForce, fdForce, wcForce;
  for (int i=0; i<8; i++)
    for (int j=0; j<8; j++)
      for (int k=0; k<8; k++){
        stForce = structuralForce(jello,i,j,k);
        shForce = sheerForce(jello,i,j,k) *1.5;
        bdForce = bendForce(jello,i,j,k) *2;
        fdForce = fieldForce(jello,i,j,k);
        wcForce = wallColForce(jello,i,j,k);
        a[i][j][k] = ( stForce + shForce + bdForce + fdForce + wcForce )/jello->mass;
        if (i==3 && j == 3 && k == 3){
          printf("\n\n a %f %f %f p %f %f %f\n field %f %f %f \n collision %f %f %f \n stforce %f %f %f\n shforce %f %f %f\n bdForce %f %f %f",
            a[i][j][k].x,a[i][j][k].y,a[i][j][k].z,jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z,
            fdForce.x,fdForce.y,fdForce.z,
            wcForce.x,wcForce.y,wcForce.z,
            stForce.x,stForce.y,stForce.z,
            shForce.x,shForce.y,shForce.z,
            bdForce.x,bdForce.y,bdForce.z);
        }
      }
}

struct point wallColForce(struct world * jello, int i, int j, int k){
  double wall[6][4] = {
    {0,1,0,2}, {0,-1,0,2},
    {1,0,0,2}, {-1,0,0,2},
    {0,0,1,2}, {0,0,-1,2}
  };

  struct point rc = {0,0,0};
  for(int index=0; index<6; index++){
    struct point cf = collisionForce(jello, i,j,k, wall[index][0],wall[index][1],wall[index][2],wall[index][3]);
    rc = rc + cf;
  }
  
  return rc;
}

struct point bendForce(struct world * jello, int i, int j, int k){
  int others[6][3] = {
    {i+2,j,k},{i-2,j,k},
    {i,j+2,k},{i,j-2,k},
    {i,j,k+2},{i,j,k-2}
  };

  struct point rc = {0,0,0};
  struct point hf,df;
  for(int index=0; index < 6; index++){
    hf = hookForce(jello, i,j,k, others[index][0],others[index][1],others[index][2]);
    df = dampingForce(jello, i,j,k, others[index][0],others[index][1],others[index][2]);
    rc = hf + df + rc ;
  }
  
  return rc;
}

struct point sheerForce(struct world * jello, int i, int j, int k){
  int others[20][3] = {
    {i+1,j+1,k}, {i+1,j-1,k}, {i-1,j+1,k}, {i-1,j-1,k},
    {i+1,j,k+1}, {i+1,j,k-1}, {i-1,j,k+1}, {i-1,j,k-1},
    {i,j+1,k+1}, {i,j-1,k+1}, {i,j+1,k-1}, {i,j-1,k-1},
    {i+1,j+1,k+1}, {i+1,j+1,k-1}, {i+1,j-1,k+1}, {i+1,j-1,k-1},
    {i-1,j+1,k+1}, {i-1,j+1,k-1}, {i-1,j-1,k+1}, {i-1,j-1,k-1}
  };

  struct point rc = {0,0,0};
  struct point hf,df; 
  for(int index=0; index < 20; index++){
    hf = hookForce(jello, i,j,k, others[index][0],others[index][1],others[index][2]);
    df = dampingForce(jello, i,j,k, others[index][0],others[index][1],others[index][2]);
    rc = hf + df + rc ;
  }
  return rc;
}

struct point structuralForce(struct world * jello, int i, int j, int k){
  //set up 6 directional springs
  int others[6][3] = {
    {i+1,j,k},{i-1,j,k},
    {i,j+1,k},{i,j-1,k},
    {i,j,k+1},{i,j,k-1}
  };

  struct point rc = {0,0,0};
  struct point hf,df;
  for(int index=0; index < 6; index++){
    hf = hookForce(jello, i,j,k, others[index][0],others[index][1],others[index][2]);
    df = dampingForce(jello, i,j,k,others[index][0],others[index][1],others[index][2]);
    rc = hf + df + rc ;
  }
  
  return rc;
}

struct point fieldForce(struct world * jello, int i, int j, int k){
  
  //if the point is out of the bounding box
  if(
    jello->p[i][j][k].x > 2 ||  jello->p[i][j][k].x < -2 ||
    jello->p[i][j][k].y > 2 ||  jello->p[i][j][k].y < -2 ||
    jello->p[i][j][k].z > 2 ||  jello->p[i][j][k].z < -2 
  ) {
    struct point rc = {0,0,0};
    return rc;
  }

  //get interpolation index in the force field
  double d = 4.0/(jello->resolution-1);
  int imin = ceil((jello->p[i][j][k].x+2)/d);
  int jmin = ceil((jello->p[i][j][k].y+2)/d);
  int kmin = ceil((jello->p[i][j][k].z+2)/d);

  //perform trilinear interpolation
  struct point c000 = jello->forceField[imin*jello->resolution*jello->resolution + jmin*jello->resolution + kmin];
  struct point c001 = jello->forceField[imin*jello->resolution*jello->resolution + jmin*jello->resolution + kmin+1];
  struct point c010 = jello->forceField[imin*jello->resolution*jello->resolution + (jmin+1)*jello->resolution + kmin];
  struct point c011 = jello->forceField[imin*jello->resolution*jello->resolution + (jmin+1)*jello->resolution + kmin+1];
  struct point c100 = jello->forceField[(imin+1)*jello->resolution*jello->resolution + jmin*jello->resolution + kmin];
  struct point c101 = jello->forceField[(imin+1)*jello->resolution*jello->resolution + jmin*jello->resolution + kmin+1];
  struct point c110 = jello->forceField[(imin+1)*jello->resolution*jello->resolution + (jmin+1)*jello->resolution + kmin];
  struct point c111 = jello->forceField[(imin+1)*jello->resolution*jello->resolution + (jmin+1)*jello->resolution + kmin+1];

  double xd = (jello->p[i][j][k].x+2 - imin * d)/d;
  double yd = (jello->p[i][j][k].y+2 - jmin * d)/d;
  double zd = (jello->p[i][j][k].z+2 - kmin * d)/d;

  struct point c00 = c000*(1-xd); c00 = c100*xd + c00;
  struct point c01 = c001*(1-xd); c01 = c101*xd + c01;
  struct point c10 = c010*(1-xd); c10 = c110*xd + c10;
  struct point c11 = c011*(1-xd); c11 = c111*xd + c11;

  struct point c0 = c00*(1-yd); c0 = c10*yd + c0;
  struct point c1 = c01*(1-yd); c1 = c11*yd + c1;
  
  struct point rc = c0*(1-zd); rc = c1*zd + rc;

  return rc;
}

struct point collisionForce(struct world * jello, int i, int j, int k, double a, double b, double c, double d){
  struct point n = {a,b,c};
  double distance = (a*jello->p[i][j][k].x+b*jello->p[i][j][k].y+c*jello->p[i][j][k].z+d) / jello->p[i][j][k].norm();
  if(distance >= 0){
    struct point rc = {0,0,0};
    return rc;
  }
  
  struct point hk = n.normalize() * (-jello->kCollision * distance);
  struct point dp = jello->v[i][j][k] * (- jello->dCollision);
  return hk + dp;
}

struct point hookForce(struct world * jello, int i1, int j1, int k1, int i2, int j2, int k2){
  if (
    i1 < 0 || i1 > 7 ||
    i2 < 0 || i2 > 7 ||
    j1 < 0 || j1 > 7 ||
    j2 < 0 || j2 > 7 ||
    k1 < 0 || k1 > 7 ||
    k2 < 0 || k2 > 7 ) {
      struct point rc = {0,0,0};
      return rc;
    }

  struct point direction = jello->p[i1][j1][k1] - jello->p[i2][j2][k2];
  double restLength = pow((double) (pow(i1-i2,2) + pow(j1-j2,2) + pow(k1-k2,2)), 0.5)/7;
  struct point rc;
  rc =  direction.normalize() * ( -jello->kElastic * (direction.norm() - restLength)  );
  return rc;
}

struct point dampingForce(struct world * jello, int i1, int j1, int k1, int i2, int j2, int k2){
  if ( 
    i1 < 0 || i1 > 7 ||
    i2 < 0 || i2 > 7 ||
    j1 < 0 || j1 > 7 ||
    j2 < 0 || j2 > 7 ||
    k1 < 0 || k1 > 7 ||
    k2 < 0 || k2 > 7 ) {
      struct point rc = {0,0,0};
      return rc;
    }

  struct point directionN = (jello->p[i1][j1][k1] - jello->p[i2][j2][k2]).normalize();
  struct point dv = jello->v[i1][j1][k1] - jello->v[i2][j2][k2];
  struct point rc = (dv.cross(directionN)).cross(directionN) * (jello->dElastic);
  return rc;
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
