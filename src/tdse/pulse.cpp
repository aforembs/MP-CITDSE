#include "pulse.hpp"

void pulse::toAU(double IoW, double weV, double &IoAU, double &wAU) {
  IoAU = IoW * conv::I_W_cm2_au_;
  wAU = weV * conv::En_ev_au_;
}

double pulse::period(double w, int cycles) { return cycles * 2.0 * M_PI / w; }

void pulse::gaussESetup(double Io, double w, double tau, int cycles,
                        pulse::params &pars) {
  constexpr double width = 0.58870501125773734549;
  pars.Eo = sqrt(Io);
  pars.w = w;
  double To = 2 * M_PI / w;
  int cycles_fwhm = 2 + int(cycles / 10);
  pars.cepd = 0.5 * To * cycles_fwhm / width;
  pars.Wenv = tau; // pass T to pulse::gauss while fitting function template
}

void pulse::gaussASetup(double Io, double w, double tau, int cycles,
                        pulse::params &pars) {
  constexpr double width = 0.58870501125773734549;
  pars.Eo = sqrt(Io) / w;
  pars.w = w;
  double To = 2 * M_PI / w;
  int cycles_fwhm = 2 + int(cycles / 10);
  pars.cepd = 0.5 * To * cycles_fwhm / width;
  pars.Wenv = tau; // pass T to pulse::gauss while fitting function template
}

double pulse::gauss(pulse::params &pars, double t) {
  double t_diff = t - 0.5 * pars.Wenv;
  return (double)(t >= 0.0 && t <= pars.Wenv) * pars.Eo *
         exp(-t_diff * t_diff / (pars.cepd * pars.cepd)) * sin(pars.w * t);
}

void pulse::sineASetup(double Io, double w, double cepd, int cycles,
                       pulse::params &pars) {
  pars.Eo = sqrt(Io) / w;
  pars.w = w;
  pars.cepd = cepd * (2.0 * M_PI / w);
  pars.Wenv = w / (2.0 * cycles);
}

// double Ao = sqrt(Io)/w;
// double cepds = cepd*(2.0*M_PI/w);
// double Wenv = w/(2.0*cycles);
double pulse::sineAA(pulse::params &pars, double t) {
  double sinWt = sin(pars.Wenv * t);
  return (double)(t >= 0.0 && t <= M_PI / pars.Wenv) * pars.Eo * sinWt * sinWt *
         sin(pars.w * t + pars.cepd);
}

// double Ao = sqrt(Io)/w;
// double cepds = cepd*(2.0*M_PI/w);
// double Wenv = w/(2.0*cycles);
double pulse::sineAE(pulse::params &pars, double t) {
  double Wt = pars.Wenv * t;
  double sinWt = sin(Wt);
  double wtc = pars.w * t + pars.cepd;
  return (double)(t >= 0.0 && t <= M_PI / pars.Wenv) * pars.Eo *
         (pars.w * sinWt * sinWt * cos(wtc) +
          pars.Wenv * sin(2 * Wt) * sin(wtc));
}

void pulse::sineESetup(double Io, double w, double cepd, int cycles,
                       pulse::params &pars) {
  pars.Eo = sqrt(Io);
  pars.w = w;
  pars.cepd = cepd * (2 * M_PI / w);
  pars.Wenv = w / (2.0 * cycles);
}

// double Eo = sqrt(Io);
// double cepd = cepd * (2*M_PI/w);
// double Wenv = w/(2.0*cycles);
// double WpP = w + 2.0*Wenv;
// double WpN = w - 2.0*Wenv;
double pulse::sineEA(pulse::params &pars, double t) {
  double WpP = pars.w + 2.0 * pars.Wenv;
  double WpN = pars.w - 2.0 * pars.Wenv;
  return (double)(t >= 0.0 && t <= M_PI / pars.Wenv) * 0.5 * pars.Eo *
         (((1.0 - cos(pars.w * t + pars.cepd)) / pars.w) -
          0.5 * ((1.0 - cos(WpP * t)) / WpP + (1.0 - cos(WpN * t)) / WpN));
}

// double Eo = sqrt(Io);
// double cepds = cepd*(2.0*M_PI/w);
// double Wenv = w/(2.0*cycles);
double pulse::sineEE(pulse::params &pars, double t) {
  double sinWt = sin(pars.Wenv * t);
  return (double)(t >= 0.0 && t <= M_PI / pars.Wenv) * pars.Eo * sinWt * sinWt *
         sin(pars.w * t + pars.cepd);
}
