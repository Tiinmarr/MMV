#pragma once

#include "mathematics.h"
#include "node.h" 

//Sphere tracing :
Vector intersection(Vector origine, Vector direction, Noeud* surface, int step) {
  Vector p = origine;
    for (int i = 0; i <step; i++) {
      if (surface->Value(p) < 0 ) {
        return p;
      }
      p.operator+=(direction);
    }
    return p;
  }