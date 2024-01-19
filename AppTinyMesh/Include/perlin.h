// Fonction de Weierstrass

#include <cmath>
#include <random>
#include "SimplexNoise.h"

class FBM {
public:
    FBM(int alpha, double lacunarity, double persistence)
        : alpha(alpha), lacunarity(lacunarity), persistence(persistence) {}

    double weierstrass(double x,double y, double a, double b) {
        float z= 0;
        for (int i = 0; i < 10; i++) {
            z += a * std::cos(b * std::pow(x, 2 * i) * std::pow(y, 2 * i));
        }
        return z;
    }

    double generate(double x, double y) {
        double value = 0.0;

        for (int i = 0; i < 10; i++) {
        double lambda = lacunarity * std::pow(persistence, i);
        double dx = lambda * x;
        double dy = lambda * y;
        double t = weierstrass(dx,dy, alpha, 1);
        value += alpha * t / lambda;
        }
    return value;
    }
private:
    int alpha;
    double lacunarity;
    double persistence;
};


// class RidgeNoise 
// {
// public : 
//     RidgeNoise(int octaves, float lacunarity, float gain) : octaves(octaves), lacunarity(lacunarity), gain(gain)
//     {
//         FastNoiseLite noise;
//         noise.SetFractalOctaves(octaves);
//         noise.SetFractalLacunarity(lacunarity);
//         noise.SetFractalGain(gain);
//     }

//     float ridgenoise(float x, float y) {
//         return 2 * (0.5 - std::abs(0.5 - noise.GetNoise(x, y)));
//     }

// private :
// int octaves;
// float lacunarity;
// float gain;
// FastNoiseLite noise;
// };

class RidgeNoise 
{
public : 
    RidgeNoise(int o)
    {
        octaves = o;
    }

    float ridgenoise(float x, float y) {
        SimplexNoise gen;
        return 2*(1 - std::abs(gen.noise(x, y))) -0.2;
    }

    float sharpenednoise(float x, float y) {
        SimplexNoise gen;
        return 2*std::pow((1 - std::abs(gen.noise(x, y))),2) -0.2;
    }

private :
int octaves;
};