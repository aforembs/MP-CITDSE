#ifndef PROP_INC_AU_H_
#define PROP_INC_AU_H_

/**
 * @brief Namespace containing constants in atomic units
 */
namespace au {
static constexpr double m_pi_ = 3.14159265358979323846;
static constexpr double m_e_ = 1.0;           // mass of electron
static constexpr double m_p_ = 1836.15270137; // proton  mass
static constexpr double m_n_ = 1838.68366240; // neutron mass
static constexpr double e_ = 1.0;             // electron charge
static constexpr double h_ = 1.0;             // Planck's constant (h/2pi)
static constexpr double c_ = 137.035989561;   // speed of light in vacuum
static constexpr double a_ = 1.0 / c_;        // fine stucture constant
static constexpr double to_ = 1.0;            // atomic time unit
static constexpr double ao_ = 1.0;            // atomic length unit(Bohr radius)
static constexpr double eo_ = 1. / 4 * m_pi_; // permittivity of vacuum
static constexpr double mo_ = 4 * m_pi_ / (c_ * c_); // permeability of vacuum
static constexpr double Ryd_ = 0.5;
static constexpr double GEn_H_ = -0.5;
// note that eo_*mo_*c_*c_=k*k where k depends on the unit system
} // namespace au

/**
 * @brief Namespace containing conversion factors to and from atomic units
 */
namespace conv {
static constexpr double En_Ryd_au_ = 0.5;
static constexpr double En_ev_au_ =
    1. / 27.211396181; // AMO HANDBOOK, G.DRAGE API 1996
static constexpr double I_W_cm2_au_ = (1. / 3.509338) * 1.e-16;
static constexpr double t_s_au_ = (1. / 2.41889) * 1.e17;
static constexpr double t_au_fs_ = 100. / 2.41889;
static constexpr double ev_J_ = 1.6021773349 * 1.e-19;
static constexpr double enau = 27.211396181; // AMO HANDBOOK, G.DRAGE API 1996
} // namespace conv

#endif // PROP_INC_AU_H_