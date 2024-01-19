#pragma once
#include "mathematics.h"
#include "node.h" 

class Scale : public Noeud {
  public :
  Noeud* enfant;
  float point;

  Scale(Noeud* g, const float c) {
    point= c;
    enfant =g;
  }

  float Value(const Vector& p) override {
    Vector v=p;
    return enfant->Value(v/point)*point;
  }
};