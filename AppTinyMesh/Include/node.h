#pragma once

#include "mathematics.h"

class Noeud {
public:
    virtual float Value(const Vector& p) = 0;
};

class OperateurBinaire : public Noeud {
public:
    Noeud* gauche;
    Noeud* droit;
};

class Arbre : public Noeud {
public:
    Noeud* racine;
    float Value(const Vector& point) override {
        return racine->Value(point);
    }
};