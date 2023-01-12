#ifndef PULSE_H_
#define PULSE_H_

#include "au.hpp"
#include <cmath>

namespace pulse {
void ToAU(double IoW, double weV, double &IoAU, double &wAU);

double Sine_T(double w, int cycles);

void SineA_Setup(double Io, double w, double cepd, int cycles, double &Ao,
                 double &cepds, double &Wenv);

double SineA_A(double Ao, double w, double cepds, double Wenv, double t);

double SineA_E(double Ao, double w, double cepds, double Wenv, double t);

void SineE_Setup(double Io, double w, double cepd, int cycles, double &Eo,
                 double &cepds, double &Wenv);

double SineE_A(double Eo, double w, double cepds, double Wenv, double t);

double SineE_E(double Eo, double w, double cepds, double Wenv, double t);
} // namespace pulse

#endif // PULSE_H_