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

/**
 * @brief Structure for storing the static parameters of the pulse
 *
 */
struct params {
  double Eo;
  double w;
  double cepd;
  double a;
  double b;
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
 * @brief Function setting up the parameters for a gaussian pulse defined by an
 * electric field
 *
 * @param Io peak intensity in a.u.
 * @param w photon energy in a.u.
 * @param tau the period of the pulse envelope
 * @param cycles the total number of cycles in the pulse
 * @param cepd phase offset
 * @param pars structure containing the remaining static pulse parameters
 */
void gaussESetup(double Io, double w, double tau, int cycles, double cepd,
                 pulse::params &pars);

/**
 * @brief Function setting up the parameters for a gaussian pulse defined by a
 * vector potential
 *
 * @param Io peak intensity in a.u.
 * @param w photon energy in a.u.
 * @param tau the period of the pulse envelope
 * @param cycles the total number of cycles in the pulse
 * @param cepd phase offset
 * @param pars structure containing the remaining static pulse parameters
 */
void gaussASetup(double Io, double w, double tau, int cycles, double cepd,
                 pulse::params &pars);

/**
 * @brief
 *
 * @param pars structure containing static pulse parameters
 * @param t time in a.u. since the start of the pulse
 * @return double the value of the pulse E(t) or A(t)
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
 * @param pars structure containing the remaining static pulse parameters
 */
void sineASetup(double Io, double w, double cepd, int cycles,
                pulse::params &pars);

/**
 * @brief Function returning the vector potential of the pulse at time 't' since
 * it's start.
 *
 * @param pars structure containing static pulse parameters
 * @param t time in a.u. since the start of the pulse
 * @return double the value of the vector potential at time t
 */
double sineAA(pulse::params &pars, double t);

/**
 * @brief Function returning the electric field of a pulse defined using a
 * vector potential at time 't' since it's start.
 *
 * @param pars structure containing static pulse parameters
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
 * @param pars structure containing the remaining static pulse parameters
 */
void sineESetup(double Io, double w, double cepd, int cycles,
                pulse::params &pars);

/**
 * @brief Function returning the vector potential of a pulse defined using an
 * electric field at time 't' since it's start.
 *
 * @param pars structure containing static pulse parameters
 * @param t time in a.u. since the start of the pulse
 * @return double the value of the vector potential at time t
 */
double sineEA(pulse::params &pars, double t);

/**
 * @brief Function returning electric field of the pulse at time 't' since it's
 * start.
 *
 * @param pars structure containing static pulse parameters
 * @param t time in a.u. since the start of the pulse
 * @return double the value of the electric field at time t
 */
double sineEE(pulse::params &pars, double t);
} // namespace pulse

#endif // PULSE_H_
