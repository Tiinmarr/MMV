// Box

#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"

class Box2 {
public:
    Vec2 min, max;

    Box2(const Vec2& min, const Vec2& max) : min(min), max(max) {}

    bool Inside(const Vec2& point) const {
        if (operator<=(min,point) + operator>=(max,point) == 2 ){
          return true;
        }
        else {
          return false;
        } 
    }

    bool Intersect(const Box2& other) const {
        if (operator>=(max,other.min) + operator<=(min,other.min)  == 2 ){
          return true;
        }
        else {
          return false;
        }
    }

};


class Grid : public Box2 {
public:
    int rows, cols;

    Grid(const Vec2& min, const Vec2& max, int rows, int cols) : Box2(min, max), rows(rows), cols(cols) {}

    int Index(int i, int j) const {
        // Calcul de l'index dans le tableau en fonction des coordonnées (i, j)
        return i * cols + j;
    }

    bool Inside(int i, int j) const {
        // Vérifier si les coordonnées (i, j) sont à l'intérieur des bornes de la grille
        return (i >= 0 && i < rows && j >= 0 && j < cols);
    }
};