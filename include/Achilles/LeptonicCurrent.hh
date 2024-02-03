#pragma once

#include "Achilles/Current.hh"

#include <complex>
#include <map>
#include <vector>

namespace achilles {

class FourVector;
class Particle;
class Nucleus;
class Event;
class NuclearModel;
class PID;
struct FormFactorInfo;
class ProcessInfo;

using Particles = std::vector<Particle>;
using Current = std::vector<VCurrent>;
using Currents = std::map<int, Current>;
using FFDictionary = std::map<std::pair<PID, PID>, std::vector<FormFactorInfo>>;




// general 2->2 scattering current
// TODO: is there a way to support 2->3? 2->N?
class LeptonicCurrent {
  public:
    LeptonicCurrentGeneral() = default;
    virtual void Initialize(const ProcessInfo);
    virtual FFDictionary GetFormFactor();
    virtual Currents CalcCurrents(const FourVector &, const FourVector &) const;
    virtual size_t NSpins() const {};
};




class NeutrinoCurrent : public LeptonicCurrent {
  public:
    NeutrinoCurrent() = default;
    void override Initialize(const ProcessInfo &);
    FFDictionary override GetFormFactor();
    Currents override CalcCurrents(const FourVector &, const FourVector &) const;
    size_t NSpins() const override { return 4; }

  private:
    bool NeutralCurrent(PID, PID) const;
    bool ChargedCurrent(bool, PID, PID) const;
    std::complex<double> coupl_left{}, coupl_right{};
    double mass{}, width{};
    int pid{};
    bool anti{};
};




class PrimakoffCurrent : public LeptonicCurrent {
  public:
    PrimakoffCurrent() = default;
    void override Initialize(const ProcessInfo &, const double &, const double &, const double &, const double &);
    FFDictionary override GetFormFactor();
    Currents override CalcCurrents(const FourVector &, const FourVector &) const;
    size_t NSpins() const override { return 2; }   // TODO(AT): what function?

  private:
    double coupling_dim5{}, alp_mass{}, vector_med_mass{}, vector_med_width{};
    int pid{};

};

} // namespace achilles
