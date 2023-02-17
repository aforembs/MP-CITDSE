#ifndef PULSE_H_
#define PULSE_H_

/**
 * @file pulse.hpp
 * @author Andrew Forembski (andrew.forembski2@mail.dcu.ie)
 * @brief Functions used for modelling the laser pulse
 * @version 1.0
 */

#include "au.hpp"
#include <cmath>

namespace pulse {

struct params {
  double Eo;
  double w;
  double cepd;
  double Wenv;
};

/**
 * @brief Function for converting eV and W/cm^2 to atomic units
 *
 * @param IoW peak intensity in W/cm^2
 * @param weV photon energy in eV
 * @param IoAU peak intensity in atomic units
 * @param wAU photon energy in atomic units
 */
void toAU(double IoW, double weV, double &IoAU, double &wAU);

/**
 * @brief Convert cycles to pulse duration
 *
 * @param w photon energy
 * @param cycles number of cycles
 * @return double pulse duration
 */
double period(double w, int cycles);

/**
 *
 *
 */
void gaussESetup(double Io, double w, double tau, int cycles,
                 pulse::params &pars);

/**
 * @brief
 *
 * @param Io
 * @param w
 * @param tau
 * @param cycles
 * @param pars
 */
void gaussASetup(double Io, double w, double tau, int cycles,
                 pulse::params &pars);

/**
 * @brief
 *
 * @param pars
 * @param t
 * @return double
 */
double gauss(pulse::params &pars, double t);

/**
 * @brief Function setting up the parameters for a pulse defined by a vector
 * potential
 *
 * @param Io peak intensity in a.u.
 * @param w photon energy in a.u.
 * @param cepd phase offset
 * @param cycles number of cycles
 * @param Ao peak vector potential in a.u.
 * @param cepds adjusted pulse frequency
 * @param Wenv pulse envelope frequency
 */
void sineASetup(double Io, double w, double cepd, int cycles,
                pulse::params &pars);

/**
 * @brief Function returning the vector potential of the pulse at time 't' since
 * it's start.
 *
 * @param Ao peak vector potential in a.u.
 * @param w photon energy in a.u.
 * @param cepds adjusted pulse frequency
 * @param Wenv pulse envelope frequency
 * @param t time in a.u. since the start of the pulse
 * @return double the value of the vector potential at time t
 */
double sineAA(pulse::params &pars, double t);

/**
 * @brief Function returning the electric field of a pulse defined using a
 * vector potential at time 't' since it's start.
 *
 * @param Ao peak vector potential in a.u.
 * @param w photon energy in a.u.
 * @param cepds adjusted pulse frequency
 * @param Wenv pulse envelope frequency
 * @param t time in a.u. since the start of the pulse
 * @return double the value of the electric field at time t
 */
double sineAE(pulse::params &pars, double t);

/**
 * @brief Setup of a pulse defined by an electric field
 *
 * @param Io peak intensity in a.u.
 * @param w photon energy in a.u.
 * @param cepd phase offset
 * @param cycles number of cycles
 * @param Eo peak electric field in a.u.
 * @param cepds adjusted pulse frequency
 * @param Wenv pulse envelope frequency
 */
void sineESetup(double Io, double w, double cepd, int cycles,
                pulse::params &pars);

/**
 * @brief Function returning the vector potential of a pulse defined using an
 * electric field at time 't' since it's start.
 *
 * @param Eo peak electric field in a.u.
 * @param w photon energy in a.u.
 * @param cepds adjusted pulse frequency
 * @param Wenv pulse envelope frequency
 * @param t time in a.u. since the start of the pulse
 * @return double the value of the vector potential at time t
 */
double sineEA(pulse::params &pars, double t);

/**
 * @brief Function returning electric field of the pulse at time 't' since it's
 * start.
 *
 * @param Ao peak electric field in a.u.
 * @param w photon energy in a.u.
 * @param cepds adjusted pulse frequency
 * @param Wenv pulse envelope frequency
 * @param t time in a.u. since the start of the pulse
 * @return double the value of the electric field at time t
 */
double sineEE(pulse::params &pars, double t);
} // namespace pulse

#endif // PULSE_H_
