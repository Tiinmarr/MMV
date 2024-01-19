// HeightField.cpp
#include "HeightField.h"
#include <fstream>
#include <random>
#include <cmath>

Vector RGB(float t) {
    if (t < 0.5) return Vector(255*2*t,255*2*t,255);
    else return Vector(255 ,255 * 2 * (1 - t) ,255 * 2 * (1 - t));
}

Vector RGB2(float t) {
    return Vector(255*t,255*t,255);
}

double access(const HeightField& H, int i, int j) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> al(0,360);
    std::uniform_real_distribution<double> th(0,90);
    int result = 0;
    for (int k=0; k<10; k++) {
        int iter = 0;
        double alpha = al(gen);
        double theta = th(gen);
        while(iter <5){
            int x = i + cos(alpha) * sin(theta);
            int y = j + sin(alpha) * sin(theta);
            double z = H.Height(i,j) + cos(theta);
            if (x < 0 || x >= H.rows || y < 0 || y >= H.cols) {
                result += 1;
                break;
            }
        
            if (z < H.Height(x,y)) {
                break;
            }
            result += 1;
            iter += 1;
        }
    }
    return result / 10.0;
}

HeightField::HeightField(const Vec2& min, const Vec2& max, int rows, int cols)
    : ScalarField(min, max, rows, cols) {
    // Initialisation spécifique à HeightField si nécessaire
}

double HeightField::Height(int i, int j) const {
    return values[Index(i, j)];
}

double HeightField::Slope(int i, int j) const {
    Vec2 gradient = Gradient(i, j);
    // return std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]);
    return std::atan(std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]));
}

double HeightField::AverageSlope(int i, int j) const {
        double totalSlope = 0.0;

    for (int di = -1; di <= 1; ++di) {
        for (int dj = -1; dj <= 1; ++dj) {
            if (di == 0 && dj == 0) continue;

            totalSlope += Slope(i + di, j + dj);
        }
    }

    return totalSlope / 8.0;
}

Vector HeightField::Vertex(int i, int j) const {
    double x = i;
    double y = j;
    double z = Height(i, j);

    return Vector(x, y, z);
}

Vector HeightField::Normal(int i, int j) const {
    Vec2 gradient = Gradient(i, j);
    return Normalized(Vector(-gradient[0], -gradient[1], 1.0));
}

QImage HeightField::Shade(const Vector& lightDirection) const {
    QImage image(cols, rows, QImage::Format_RGB32);
    double minElevation = std::numeric_limits<double>::max();
    double maxElevation = std::numeric_limits<double>::min();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double elevation = values[Index(i, j)];
            minElevation = std::min(minElevation, elevation);
            maxElevation = std::max(maxElevation, elevation);
        }
    }

    double elevationRange = maxElevation - minElevation;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Vector normal = Normal(i, j);
            // double intensity = Math::Max(0.0, SquaredNorm2(normal,lightDirection));
            double normalizedElevation = (values[Index(i, j)] - minElevation) / elevationRange;
            Vector color = RGB(normalizedElevation);
            image.setPixel(j, i, qRgb(color[0], color[1], color[2]));
        }
    }
    return image;
}

void HeightField::Export(const QString& filename, const Vector& lightDirection) const {
    // std::ofstream objFile(filename);
    QImage shadedImage = Shade(lightDirection);
    shadedImage.save(filename);
}

void HeightField::Export_Slope(const QString& filename) const {
    QImage image(cols, rows, QImage::Format_RGB32);
    double minSlope = std::numeric_limits<double>::max();
    double maxSlope = std::numeric_limits<double>::min();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double slope = Slope(i, j);
            minSlope = std::min(minSlope, slope);
            maxSlope = std::max(maxSlope, slope);
        }
    }

    double SlopeRange = maxSlope - minSlope;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double normalizedSlope = (Slope(i,j) - minSlope) / SlopeRange;
            Vector color = RGB2(normalizedSlope);
            image.setPixel(j, i, qRgb(color[0], color[1], color[2]));
        }
    }
    image.save(filename);
}

void HeightField::Export_Laplacian(const QString& filename) const {
    QImage image(cols, rows, QImage::Format_RGB32);
    double minL = std::numeric_limits<double>::max();
    double maxL = std::numeric_limits<double>::min();

    ScalarField L = Laplacian();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double laplacian = L.values[Index(i, j)];
            minL = std::min(minL, laplacian);
            maxL = std::max(maxL, laplacian);
        }
    }

    double LRange = maxL - minL;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double normalizedL = (L.values[Index(i, j)] - minL) / LRange;
            // Vector color = RGB(normalizedL);
            // image.setPixel(j, i, qRgb(color[0], color[1], color[2]));
            double value = normalizedL * 255.0;
            int color = static_cast<int>(value);
            image.setPixel(j, i, qRgb(color, 0, 128));
        }
    }
    image.save(filename);
}

void HeightField::Export_Accesibility(const QString& filename) const {
// Determine the accesibility in each point of the Heightfield :
    ScalarField A(min,max,rows,cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            A.values[Index(i, j)] = access(*this,i,j);
        }
    }

    QImage image(cols, rows, QImage::Format_RGB32);
    double minA = std::numeric_limits<double>::max();
    double maxA = std::numeric_limits<double>::min();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double acc = A.values[Index(i, j)];
            minA = std::min(minA, acc);
            maxA = std::max(maxA, acc);
        }
    }

    double ARange = maxA - minA;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double normalizedA = (A.values[Index(i, j)] - minA) / ARange;
            Vector color = RGB2(normalizedA);
            image.setPixel(j, i, qRgb(color[0], color[1], color[2]));
        }
    }
    image.save(filename);

}

