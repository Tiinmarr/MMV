// HeightField.h
#pragma once

#include "ScalarField.h"
#include "mathematics.h"
#include <vector>

class HeightField : public ScalarField {
public:
    HeightField(const Vec2& min, const Vec2& max, int rows, int cols);

    double Height(int i, int j) const;
    double Slope(int i, int j) const;
    double AverageSlope(int i, int j) const;

    Vector Vertex(int i, int j) const;
    Vector Normal(int i, int j) const;
    QImage Shade(const Vector& lightDirection) const;
    void Export(const QString& filename,const Vector& lightDirection) const;
    void Export_Slope(const QString& filename) const;
    void Export_Laplacian(const QString& filename) const;
    void Export_Accesibility(const QString& filename) const;
};
