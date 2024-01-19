// ScalarField 

// Self include
#include "ScalarField.h"
#include <iostream>

ScalarField::ScalarField(const Vec2& min, const Vec2& max, int rows, int cols)
    : Grid(min, max, rows, cols) {
    values.resize(rows * cols, 0.0);
}

void ScalarField::saveScalarFieldImage(const QString& filename) const {
    QImage image(cols, rows, QImage::Format_RGB32);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double value = values[Index(i, j)];
            int color = static_cast<int>(value * 255.0);
            image.setPixel(j, i, qRgb(color, color, color));
        }
    }

    image.save(filename);
}

Vec2 ScalarField::Gradient(int i, int j) const {
    double dx = (i > 0 && i < rows - 1) ? (values[Index(i + 1, j)] - values[Index(i - 1, j)]) * 0.5 : 0.0;
    double dy = (j > 0 && j < cols - 1) ? (values[Index(i, j + 1)] - values[Index(i, j - 1)]) * 0.5 : 0.0;
    return Vec2(dx, dy);
}

ScalarField ScalarField::GradientNorm() const {
        ScalarField result(min, max, rows, cols);

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                Vec2 gradient = Gradient(i, j);
                result.values[Index(i, j)] = std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]);
            }
        }
  return result;
}

ScalarField ScalarField::Laplacian() const {
    ScalarField result(min, max, rows, cols);

    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            double laplacian = values[Index(i + 1, j)] + values[Index(i - 1, j)] +
                                values[Index(i, j + 1)] + values[Index(i, j - 1)] -
                                4.0 * values[Index(i, j)];
            result.values[Index(i, j)] = laplacian;
        }
    }

    return result;
}

  void ScalarField::Smooth() {
    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            double smoothedValue = (values[Index(i - 1, j - 1)] + 2.0 * values[Index(i - 1, j)] + values[Index(i - 1, j + 1)] +
                                    2.0 * values[Index(i, j - 1)] + 4.0 * values[Index(i, j)] + 2.0 * values[Index(i, j + 1)] +
                                    values[Index(i + 1, j - 1)] + 2.0 * values[Index(i + 1, j)] + values[Index(i + 1, j + 1)]) / 16.0;
            values[Index(i, j)] = smoothedValue;
        }
    }
}
