#ifndef FORMFACTORY_HH
#define FORMFACTORY_HH

#include "Achilles/Factory.hh"

#include <array>
#include <complex>
#include <memory>

#include "Achilles/Achilles.hh"

#include "fmt/format.h"

#include "yaml-cpp/node/node.h"

namespace achilles {

class FormFactorImpl;
class FormFactorBuilder;

struct FormFactorInfo {
    enum class Type {
        F1,
        F2,
        F1p,
        F1n,
        F2p,
        F2n,
        FA,
        FCoh,
    };

    Type form_factor;
    std::complex<double> coupling{};

    bool operator==(const FormFactorInfo &other) const {
        return form_factor == other.form_factor && coupling == other.coupling;
    }
};

class FormFactor {
  public:
    struct Values {
        double Gep{}, Gen{}, Gmp{}, Gmn{};
        double F1p{}, F2p{}, F1n{}, F2n{};
        double FA{}, FAs{};
        double Fcoh{};
    };

    FormFactor() = default;
    MOCK ~FormFactor() = default;

    MOCK Values operator()(double Q2) const;

    friend FormFactorBuilder;

  private:
    std::shared_ptr<FormFactorImpl> vector = nullptr;
    std::shared_ptr<FormFactorImpl> axial = nullptr;
    std::shared_ptr<FormFactorImpl> coherent = nullptr;
};

enum class FFType { vector, axial, coherent };

inline std::string FFTypeToString(FFType type) {
    std::string result;
    switch(type) {
    case FFType::vector:
        result = "vector";
        break;
    case FFType::axial:
        result = "axial";
        break;
    case FFType::coherent:
        result = "coherent";
        break;
    }
    return result;
}

template <typename Derived>
using RegistrableFormFactor = Registrable<FormFactorImpl, Derived, FFType, const YAML::Node &>;
using FormFactorFactory = Factory<FormFactorImpl, FFType, const YAML::Node &>;

class FormFactorBuilder {
  public:
    FormFactorBuilder() { form_factor = std::make_unique<FormFactor>(); }
    static FormFactorBuilder &Instance() {
        static FormFactorBuilder instance;
        return instance;
    }
    MOCK ~FormFactorBuilder() = default;
    MOCK void Reset() { form_factor = std::make_unique<FormFactor>(); }
    MOCK FormFactorBuilder &Vector(const std::string &, const YAML::Node &);
    MOCK FormFactorBuilder &AxialVector(const std::string &, const YAML::Node &);
    MOCK FormFactorBuilder &Coherent(const std::string &, const YAML::Node &);

    MOCK std::unique_ptr<FormFactor> build() { return std::move(form_factor); }

  private:
    std::unique_ptr<FormFactor> form_factor = nullptr;
};

class FormFactorImpl {
  public:
    FormFactorImpl() = default;
    FormFactorImpl(const FormFactorImpl &) = default;
    FormFactorImpl(FormFactorImpl &&) = default;
    FormFactorImpl &operator=(const FormFactorImpl &) = default;
    FormFactorImpl &operator=(FormFactorImpl &&) = default;
    virtual ~FormFactorImpl() = default;
    virtual void Evaluate(double, FormFactor::Values &) const = 0;

    static std::unique_ptr<FormFactorImpl> Build(const YAML::Node &);
    static std::string Name() { return "Form Factor"; }

  protected:
    void Fill(double, FormFactor::Values &) const;
    template <typename Derived> static void Validate(FFType other) {
        if(Derived::Type() != other)
            throw std::runtime_error(fmt::format("FormFactor: Expected type {}, got type {}",
                                                 FFTypeToString(Derived::Type()),
                                                 FFTypeToString(other)));
    }
};

class VectorDipole : public FormFactorImpl, RegistrableFormFactor<VectorDipole> {
  public:
    VectorDipole(const YAML::Node &);
    void Evaluate(double, FormFactor::Values &) const override;

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType, const YAML::Node &);
    static std::string Name() { return "VectorDipole"; }
    static FFType Type() { return FFType::vector; }

  private:
    double lambda, muP, muN;
};

class VectorDummy : public FormFactorImpl, RegistrableFormFactor<VectorDummy> {
  public:
    VectorDummy(const YAML::Node &) {}
    void Evaluate(double, FormFactor::Values &vals) const override {
        vals.F1p = 1;
        vals.F1n = 0;
        vals.F2p = 0;
        vals.F2n = 0;
    }

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType type, const YAML::Node &node) {
        Validate<VectorDummy>(type);
        return std::make_unique<VectorDummy>(node);
    }
    static std::string Name() { return "VectorDummy"; }
    static FFType Type() { return FFType::vector; }
};

class AxialDummy : public FormFactorImpl, RegistrableFormFactor<AxialDummy> {
  public:
    AxialDummy(const YAML::Node &) {}
    void Evaluate(double, FormFactor::Values &vals) const override {
        vals.FA = 0;
        vals.FAs = 0;
    }

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType type, const YAML::Node &node) {
        Validate<AxialDummy>(type);
        return std::make_unique<AxialDummy>(node);
    }
    static std::string Name() { return "AxialDummy"; }
    static FFType Type() { return FFType::axial; }
};

class AxialDipole : public FormFactorImpl, RegistrableFormFactor<AxialDipole> {
  public:
    AxialDipole(const YAML::Node &);
    void Evaluate(double, FormFactor::Values &) const override;

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType, const YAML::Node &);
    static std::string Name() { return "AxialDipole"; }
    static FFType Type() { return FFType::axial; }

  private:
    double MA, gan1, gans;
};

class AxialZExpansion : public FormFactorImpl, RegistrableFormFactor<AxialZExpansion> {
  public:
    AxialZExpansion(const YAML::Node &);
    void Evaluate(double, FormFactor::Values &) const override;

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType, const YAML::Node &);
    static std::string Name() { return "AxialZExpansion"; }
    static FFType Type() { return FFType::axial; }

  private:
    double tcut, t0;
    std::vector<double> cc_params{}, strange_params{};
    double ZExpand(std::vector<double> A, double z) const {
        // Compute the z-expansion, \sum{i=0}^{i_max} a_i z^i
        double result = 0;
        for(size_t i = A.size() - 1; i > 0; --i) {
            result += A[i];
            result *= z;
        }
        return result + A[0];
    }
};

class Kelly : public FormFactorImpl, RegistrableFormFactor<Kelly> {
  public:
    Kelly(const YAML::Node &);
    void Evaluate(double, FormFactor::Values &) const override;

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType, const YAML::Node &);
    static std::string Name() { return "Kelly"; }
    static FFType Type() { return FFType::vector; }

  private:
    double lambdasq, muP, muN;
    std::array<double, 4> termsEp{};
    std::array<double, 2> termsEn{};
    std::array<double, 4> termsMp{};
    std::array<double, 4> termsMn{};

    double Parameterization(std::array<double, 4> A, double x) const {
        return (1 + A[0] * x) / (1 + A[1] * x + A[2] * x * x + A[3] * x * x * x);
    }
};

class BBBA : public FormFactorImpl, RegistrableFormFactor<BBBA> {
  public:
    BBBA(const YAML::Node &);
    void Evaluate(double, FormFactor::Values &) const override;

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType, const YAML::Node &);
    static std::string Name() { return "BBBA"; }
    static FFType Type() { return FFType::vector; }

  private:
    double muP, muN;
    std::array<double, 4> numEp{}, denEp{};
    std::array<double, 4> numEn{}, denEn{};
    std::array<double, 4> numMp{}, denMp{};
    std::array<double, 4> numMn{}, denMn{};

    double Numerator(std::array<double, 4> A, double x) const {
        double result = 0;
        for(size_t i = A.size() - 1; i > 0; --i) {
            result += A[i];
            result *= x;
        }

        return result + A[0];
    }

    double Denominator(std::array<double, 4> B, double x) const {
        double result = 0;
        for(size_t i = B.size(); i > 0; --i) {
            result += B[i - 1];
            result *= x;
        }

        return result + 1;
    }
};

class ArringtonHill : public FormFactorImpl, RegistrableFormFactor<ArringtonHill> {
  public:
    ArringtonHill(const YAML::Node &);
    void Evaluate(double, FormFactor::Values &) const override;

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType, const YAML::Node &);
    static std::string Name() { return "ArringtonHill"; }
    static FFType Type() { return FFType::vector; }

  private:
    double muP, muN, tcut, t0;
    std::array<double, 13> epParams{}, enParams{}, mpParams{}, mnParams{};

    double ZExpand(std::array<double, 13> A, double z) const {
        double result = 0;
        for(size_t i = A.size() - 1; i > 0; --i) {
            result += A[i];
            result *= z;
        }
        return result + A[0];
    }
};

class HelmFormFactor : public FormFactorImpl, RegistrableFormFactor<HelmFormFactor> {
  public:
    HelmFormFactor(const YAML::Node &);
    void Evaluate(double, FormFactor::Values &) const override;

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType, const YAML::Node &);
    static std::string Name() { return "Helm"; }
    static FFType Type() { return FFType::coherent; }

  private:
    double s, r;
};

class LovatoFormFactor : public FormFactorImpl, RegistrableFormFactor<LovatoFormFactor> {
  public:
    LovatoFormFactor(const YAML::Node &);
    void Evaluate(double, FormFactor::Values &) const override;

    // Required factory methods
    static std::unique_ptr<FormFactorImpl> Construct(FFType, const YAML::Node &);
    static std::string Name() { return "Lovato"; }
    static FFType Type() { return FFType::coherent; }

  private:
    double b;
    std::array<double, 5> c{};
};

} // namespace achilles

namespace fmt {

template <> struct formatter<achilles::FormFactorInfo::Type> {
    template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const achilles::FormFactorInfo::Type &ffit, FormatContext &ctx) const {
        switch(ffit) {
        case achilles::FormFactorInfo::Type::F1:
            return format_to(ctx.out(), "F1({})", static_cast<int>(ffit));
        case achilles::FormFactorInfo::Type::F2:
            return format_to(ctx.out(), "F2({})", static_cast<int>(ffit));
        case achilles::FormFactorInfo::Type::F1p:
            return format_to(ctx.out(), "F1p({})", static_cast<int>(ffit));
        case achilles::FormFactorInfo::Type::F1n:
            return format_to(ctx.out(), "F1n({})", static_cast<int>(ffit));
        case achilles::FormFactorInfo::Type::F2p:
            return format_to(ctx.out(), "F2p({})", static_cast<int>(ffit));
        case achilles::FormFactorInfo::Type::F2n:
            return format_to(ctx.out(), "F2n({})", static_cast<int>(ffit));
        case achilles::FormFactorInfo::Type::FA:
            return format_to(ctx.out(), "FA({})", static_cast<int>(ffit));
        case achilles::FormFactorInfo::Type::FCoh:
            return format_to(ctx.out(), "FCoh({})", static_cast<int>(ffit));
        default:
            return format_to(ctx.out(), "Unknown achilles::FormFactorInfo::Type({}) ",
                             static_cast<int>(ffit));
        }
    }
};
} // namespace fmt

#endif
