// Box

#pragma once

#include <vector>
#include <iostream>
#include <QImage>
#include <QColor>

#include "mathematics.h"
#include "box2.h"

class ScalarField : public Grid {
public:
    std::vector<double> values;

    ScalarField(const Vec2& min, const Vec2& max, int rows, int cols);

    void saveScalarFieldImage(const QString& filename) const;

    Vec2 Gradient(int i, int j) const;

    ScalarField GradientNorm() const;

    ScalarField Laplacian() const;

    void Smooth();
};