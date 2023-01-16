#include "pulse.hpp"

void pulse::toAU(double IoW, double weV, double &IoAU, double &wAU) {
  IoAU = IoW * conv::I_W_cm2_au_;
  wAU = weV * conv::En_ev_au_;
}

double pulse::sineT(double w, int cycles) { return cycles * 2.0 * M_PI / w; }

void pulse::sineASetup(double Io, double w, double cepd, int cycles, double &Ao,
                       double &cepds, double &Wenv) {
  Ao = sqrt(Io) / w;
  cepds = cepd * (2.0 * M_PI / w);
  Wenv = w / (2.0 * cycles);
}

// double Ao = sqrt(Io)/w;
// double cepds = cepd*(2.0*M_PI/w);
// double Wenv = w/(2.0*cycles);
double pulse::sineAA(double Ao, double w, double cepds, double Wenv, double t) {
  double sinWt = sin(Wenv * t);
  return (double)(t >= 0.0 && t <= M_PI / Wenv) * Ao * sinWt * sinWt *
         sin(w * t + cepds);
}

// double Ao = sqrt(Io)/w;
// double cepds = cepd*(2.0*M_PI/w);
// double Wenv = w/(2.0*cycles);
double pulse::sineAE(double Ao, double w, double cepds, double Wenv, double t) {
  double Wt = Wenv * t;
  double sinWt = sin(Wt);
  double wtc = w * t + cepds;
  return (double)(t >= 0.0 && t <= M_PI / Wenv) * Ao *
         (w * sinWt * sinWt * cos(wtc) + Wenv * sin(2 * Wt) * sin(wtc));
}

void pulse::sineESetup(double Io, double w, double cepd, int cycles, double &Eo,
                       double &cepds, double &Wenv) {
  Eo = sqrt(Io);
  cepds = cepd * (2 * M_PI / w);
  Wenv = w / (2.0 * cycles);
}

// double Eo = sqrt(Io);
// double cepd = cepd * (2*M_PI/w);
// double Wenv = w/(2.0*cycles);
// double WpP = w + 2.0*Wenv;
// double WpN = w - 2.0*Wenv;
double pulse::sineEA(double Eo, double w, double cepds, double Wenv, double t) {
  double WpP = w + 2.0 * Wenv;
  double WpN = w - 2.0 * Wenv;
  return (double)(t >= 0.0 && t <= M_PI / Wenv) * 0.5 * Eo *
         (((1.0 - cos(w * t + cepds)) / w) -
          0.5 * ((1.0 - cos(WpP * t)) / WpP + (1.0 - cos(WpN * t)) / WpN));
}

// double Eo = sqrt(Io);
// double cepds = cepd*(2.0*M_PI/w);
// double Wenv = w/(2.0*cycles);
double pulse::sineEE(double Eo, double w, double cepds, double Wenv, double t) {
  double sinWt = sin(Wenv * t);
  return (double)(t >= 0.0 && t <= M_PI / Wenv) * Eo * sinWt * sinWt *
         sin(w * t + cepds);
}