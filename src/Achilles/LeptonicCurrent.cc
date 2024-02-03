#include "Achilles/LeptonicCurrent.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/Spinor.hh"

using achilles::LeptonicCurrent;

void NeutrinoCurrent::Initialize(const ProcessInfo &process) {
    using namespace achilles::Constant;
    const std::complex<double> i(0, 1);
    // Determine process
    bool init_neutrino = ParticleInfo(process.m_leptonic.first).IsNeutrino();
    bool neutral_current = NeutralCurrent(process.m_leptonic.first, process.m_leptonic.second[0]);
    bool charged_current =
        ChargedCurrent(init_neutrino, process.m_leptonic.first, process.m_leptonic.second[0]);
    if(!neutral_current && !charged_current)
        throw std::runtime_error("HardScattering: Invalid process");

    // TODO: Define couplings correctly
    if(charged_current) {
        pid = init_neutrino ? (process.m_leptonic.first.AsInt() < 0 ? 24 : -24)
                            : (process.m_leptonic.first.AsInt() < 0 ? -24 : 24);
        coupl_right = 0;
        coupl_left = ee * i / (sw * sqrt(2));
        mass = Constant::MW;
        width = Constant::GAMW;
    } else if(neutral_current) {
        if(init_neutrino) {
            coupl_left = (cw * ee * i) / (2 * sw) + (ee * i * sw) / (2 * cw);
            coupl_right = 0;
            pid = 23;
            mass = Constant::MZ;
            width = Constant::GAMZ;
        } else {
            coupl_right = -ee * i;
            coupl_left = coupl_right;
            pid = 22;
        }
    }
    anti = process.m_leptonic.first.AsInt() < 0;
}

bool NeutrinoCurrent::NeutralCurrent(achilles::PID initial, achilles::PID final) const {
    return initial == final;
}

bool NeutrinoCurrent::ChargedCurrent(bool neutrino, achilles::PID initial,
                                     achilles::PID final) const {
    int sign = std::signbit(initial.AsInt()) ? -1 : 1;
    return initial.AsInt() - sign * (2 * neutrino - 1) == final.AsInt();
}

achilles::FFDictionary NeutrinoCurrent::GetFormFactor() {
    FFDictionary results;
    static constexpr std::complex<double> i(0, 1);
    using namespace achilles::Constant;
    // TODO: Double check form factors
    if(pid == 24) {
        const std::complex<double> coupl = Vud * ee * i / (sw * sqrt(2) * 2);
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                         {FormFactorInfo::Type::F1n, -coupl},
                                         {FormFactorInfo::Type::F2p, coupl},
                                         {FormFactorInfo::Type::F2n, -coupl},
                                         {FormFactorInfo::Type::FA, coupl}};
        results[{PID::neutron(), pid}] = {};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == -24) {
        const std::complex<double> coupl = Vud * ee * i / (sw * sqrt(2) * 2);
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                          {FormFactorInfo::Type::F1n, -coupl},
                                          {FormFactorInfo::Type::F2p, coupl},
                                          {FormFactorInfo::Type::F2n, -coupl},
                                          {FormFactorInfo::Type::FA, coupl}};
        results[{PID::proton(), pid}] = {};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == 23) {
        const std::complex<double> coupl1 = (ee * i / (4 * sin2w * cw)) * (0.5 - 2 * sin2w);
        const std::complex<double> coupl2 = (ee * i / (4 * sw * cw));
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl1},
                                         {FormFactorInfo::Type::F1n, -coupl2},
                                         {FormFactorInfo::Type::F2p, coupl1},
                                         {FormFactorInfo::Type::F2n, -coupl2},
                                         {FormFactorInfo::Type::FA, coupl2}};
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1n, coupl1},
                                          {FormFactorInfo::Type::F1p, -coupl2},
                                          {FormFactorInfo::Type::F2n, coupl1},
                                          {FormFactorInfo::Type::F2p, -coupl2},
                                          {FormFactorInfo::Type::FA, -coupl2}};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == 22) {
        const std::complex<double> coupl = i * ee;
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                         {FormFactorInfo::Type::F2p, coupl}};
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1n, coupl},
                                          {FormFactorInfo::Type::F2n, coupl}};
        results[{PID::carbon(), pid}] = {{FormFactorInfo::Type::FCoh, 6.0 * coupl}};
    } else {
        throw std::runtime_error("NeutrinoCurrent: Invalid probe");
    }

    return results;
}

achilles::Currents NeutrinoCurrent::CalcCurrents(const FourVector &p_in,
                                                 const FourVector &p_out) const {
    Currents currents;

    // Setup spinors
    FourVector pU, pUBar;
    if(anti) {
        pUBar = -p_in;
        pU = p_out;
    } else {
        pU = -p_in;
        pUBar = p_out;
    }
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pUBar);
    ubar[1] = UBarSpinor(1, pUBar);
    u[0] = USpinor(-1, pU);
    u[1] = USpinor(1, pU);

    // Calculate currents
    Current result;
    double q2 = (p_in - p_out).M2();
    std::complex<double> prop =
        std::complex<double>(0, 1) / (q2 - mass * mass - std::complex<double>(0, 1) * mass * width);
    spdlog::trace("Calculating Current for {}", pid);
    for(size_t i = 0; i < 2; ++i) {
        for(size_t j = 0; j < 2; ++j) {
            VCurrent subcur;
            for(size_t mu = 0; mu < 4; ++mu) {
                subcur[mu] = ubar[i] *
                             (coupl_left * SpinMatrix::GammaMu(mu) * SpinMatrix::PL() +
                              coupl_right * SpinMatrix::GammaMu(mu) * SpinMatrix::PR()) *
                             u[j] * prop;
                spdlog::trace("Current[{}][{}] = {}", 2 * i + j, mu, subcur[mu]);
            }
            result.push_back(subcur);
        }
    }
    currents[pid] = result;

    return currents;
}




// Initialize PrimakoffCurrent
void PrimakoffCurrent::Initialize(const ProcessInfo &process,
                                  const double &coupling
                                  const double &mass,
                                  const double &mediator_mass,
                                  const double &mediator_width) {
    
    alp_mass = mass;
    coupling_dim5 = coupling;
    vector_med_mass = mediator_mass;
    vector_med_width = mediator_width;
    pid = 81;  // TODO: use different PID? should ALP PID added to yaml file?

}

achilles::FFDictionary PrimakoffCurrent::GetFormFactor() {
    FFDictionary results;
    static constexpr std::complex<double> i(0, 1);
    using namespace achilles::Constant;

    results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, 0.0},
                                        {FormFactorInfo::Type::F1n, 0.0},
                                        {FormFactorInfo::Type::F2p, 0.0},
                                        {FormFactorInfo::Type::F2n, 0.0},
                                        {FormFactorInfo::Type::FA, coupling_dim5}};
    results[{PID::neutron(), pid}] = {};
    results[{PID::carbon(), pid}] = {};

    return results;
}



achilles::Currents PrimakoffCurrent::CalcCurrents(const FourVector &p_in,
                                                 const FourVector &p_out) const {
    Currents currents;

    // Calculate currents: incoming scalar particle, outgoing photon, vector mediator
    Current result;
    FourVector q = p_in - p_out;
    double q2 = q.M2();
    std::complex<double> prop =
        std::complex<double>(0, 1) / (q2 - vector_med_mass * vector_med_mass - std::complex<double>(0, 1) * vector_med_mass * vector_med_width);

    FourVector eps;  // placeholder <--- need PolarizationVector class to replace this for the outgoing polarization vector
    double eps_dot_q = eps * q;
    double p_out_dot_q = p_out * q;

    // TODO: initialize eps in the for loop over polarizations

    spdlog::trace("Calculating Current for {}", pid);
    for(size_t i = 0; i < 2; ++i) {  // sum over final photon polarization
        VCurrent subcur;
        for(size_t mu = 0; mu < 4; ++mu) {
            // for Primakoff, we have a massive incoming scalar/pseudoscalar with a massive vector mediator
            subcur[mu] = coupling*(eps_dot_q*p_out[mu] - p_out_dot_q*[mu]) * prop;
            spdlog::trace("Current[{}][{}] = {}", 2 * i + j, mu, subcur[mu]);
        result.push_back(subcur);
        }
    }
    currents[pid] = result;

    return currents;
}

