#include "Achilles/SpectralFunction.hh"
#include <memory>

extern "C" {
    std::shared_ptr<achilles::SpectralFunction> LoadSpectralFunction(const char *filename) {
        std::string filename_str(filename); 
        auto spectral =  std::make_shared<achilles::SpectralFunction>(filename_str);
        return spectral;
    }

    void DeleteSpectralFunction(std::shared_ptr<achilles::SpectralFunction> spectral) { }

    double SpectralNormalization(std::shared_ptr<achilles::SpectralFunction> spectral) {
        return spectral -> Normalization();
    }

    double SpectralFunction(std::shared_ptr<achilles::SpectralFunction> spectral, double p, double E) {
        return spectral -> operator()(p, E);
    }
}
