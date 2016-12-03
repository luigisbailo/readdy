#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <readdy/Simulation.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/common/nodelete.h>
#include "ExportSchemeApi.h"
#include "PyPotential.h"
#include "PyFunction.h"

namespace py = pybind11;

using rvp = py::return_value_policy;
using sim = readdy::Simulation;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::model::Vec3;
using pot2 = readdy::rpy::PotentialOrder2Wrapper;
using model = readdy::model::KernelStateModel;
using ctx = readdy::model::KernelContext;
using kern = readdy::model::Kernel;

// thin wrappers
void setBoxSize(sim &self, const vec &size) { /* explicitly choose void(vec) signature */ self.setBoxSize(size); }

std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.getSelectedKernelType(); }

void addParticle(sim &self, const std::string &type, const vec &pos) { self.addParticle(pos[0], pos[1], pos[2], type); }

void registerPotentialOrder2(sim &self, pot2 *potential, std::string type1, std::string type2) {
    self.registerPotentialOrder2(potential, type1, type2);
}

unsigned long
registerObservable_ParticlePositions(sim &self, unsigned int stride, pybind11::object callbackFun,
                                     std::vector<std::string> types) {
    auto pyFun = readdy::rpy::PyFunction<void(readdy::model::ParticlePositionObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::ParticlePositionObservable>(std::move(pyFun), stride, types);
}

unsigned long
registerObservable_RadialDistribution(sim &self, unsigned int stride, pybind11::object callbackFun,
                                      py::array_t<double> &binBorders, std::string typeCountFrom,
                                      std::string typeCountTo, double particleDensity) {
    auto pyFun = readdy::rpy::PyFunction<void(readdy::model::RadialDistributionObservable::result_t)>(callbackFun);
    const auto info = binBorders.request();
    std::vector<double> binBordersVec{};
    binBordersVec.reserve(info.shape[0]);
    const auto data = static_cast<double *>(info.ptr);
    for (auto i = 0; i < info.shape[0]; ++i) binBordersVec.push_back(data[i]);
    return self.registerObservable<readdy::model::RadialDistributionObservable>(std::move(pyFun), stride, binBordersVec,
                                                                                typeCountFrom, typeCountTo,
                                                                                particleDensity);
}

unsigned long
registerObservable_CenterOfMass(sim &self, unsigned int stride, const pybind11::object &callbackFun,
                                std::vector<std::string> types) {
    auto pyFun = readdy::rpy::PyFunction<void(readdy::model::CenterOfMassObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::CenterOfMassObservable>(
            std::move(pyFun), stride, types
    );
}

unsigned long
registerObservable_HistogramAlongAxisObservable(sim &self, unsigned int stride, const py::object &callbackFun,
                                                py::array_t<double> binBorders, unsigned int axis,
                                                std::vector<std::string> types) {
    const auto info = binBorders.request();
    const auto sizeBorders = info.shape[0];
    auto binBordersData = static_cast<double *>(info.ptr);
    std::vector<double> binBordersVec{};
    binBordersVec.reserve(sizeBorders);
    for (auto i = 0; i < sizeBorders; ++i) {
        binBordersVec.push_back(binBordersData[i]);
    }
    auto pyFun = readdy::rpy::PyFunction<void(readdy::model::HistogramAlongAxisObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::HistogramAlongAxisObservable>(std::move(pyFun), stride, binBordersVec,
                                                                                types, axis);
}

unsigned long
registerObservable_NParticles(sim &self, unsigned int stride, const py::object &callbackFun,
                              std::vector<std::string> types) {
    auto pyFun = readdy::rpy::PyFunction<void(readdy::model::NParticlesObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::NParticlesObservable>(std::move(pyFun), stride, types);
}

unsigned long registerObservable_ForcesObservable(sim &self, unsigned int stride, py::object callbackFun,
                                                  std::vector<std::string> types) {
    auto pyFun = readdy::rpy::PyFunction<void(readdy::model::ForcesObservable::result_t)>(callbackFun);
    return self.registerObservable<readdy::model::ForcesObservable>(std::move(pyFun), stride, types);
}

// module
PYBIND11_PLUGIN (api) {

    if (!readdy::log::console()) {
        spdlog::set_sync_mode();
        auto console = spdlog::stdout_color_mt("console");
        console->set_level(spdlog::level::debug);
        console->set_pattern("[          ] [%Y-%m-%d %H:%M:%S] [%t] [%l] %v");
    }

    py::module api("api", "ReaDDy c++-api python module");

    exportSchemeApi<readdy::api::ReaDDyScheme>(api, "ReaDDyScheme");

    py::class_<sim>(api, "Simulation")
            .def(py::init<>())
            .def_property("kbt", &sim::getKBT, &sim::setKBT)
            .def_property("periodic_boundary", &sim::getPeriodicBoundary, &sim::setPeriodicBoundary)
            .def_property("box_size", &sim::getBoxSize, &setBoxSize)
            .def("register_particle_type", &sim::registerParticleType)
            .def("add_particle", [](sim &self, const std::string &type, const vec &pos) {
                self.addParticle(pos[0], pos[1], pos[2], type);
            })
            .def("is_kernel_selected", &sim::isKernelSelected)
            .def("get_selected_kernel_type", &getSelectedKernelType)
            .def("register_potential_order_2", &registerPotentialOrder2)
            .def("register_potential_harmonic_repulsion", &sim::registerHarmonicRepulsionPotential)
            .def("register_potential_piecewise_weak_interaction",
                 &sim::registerWeakInteractionPiecewiseHarmonicPotential)
            .def("register_potential_box", &sim::registerBoxPotential)
            .def("register_potential_sphere", &sim::registerSpherePotential)
            .def("get_particle_positions", &sim::getParticlePositions)
            .def("register_observable_particle_positions", &registerObservable_ParticlePositions)
            .def("register_observable_radial_distribution", &registerObservable_RadialDistribution)
            .def("register_observable_histogram_along_axis", &registerObservable_HistogramAlongAxisObservable)
            .def("register_observable_center_of_mass", &registerObservable_CenterOfMass)
            .def("register_observable_n_particles", &registerObservable_NParticles)
            .def("register_observable_forces", &registerObservable_ForcesObservable)
            .def("deregister_observable", &sim::deregisterObservable)
            .def("register_reaction_conversion", &sim::registerConversionReaction, rvp::reference_internal)
            .def("register_reaction_enzymatic", &sim::registerEnzymaticReaction, rvp::reference_internal)
            .def("register_reaction_fission", &sim::registerFissionReaction, rvp::reference_internal)
            .def("register_reaction_fusion", &sim::registerFusionReaction, rvp::reference_internal)
            .def("register_reaction_decay", &sim::registerDecayReaction, rvp::reference_internal)
            .def("get_recommended_time_step", &sim::getRecommendedTimeStep)
            .def("set_kernel", &sim::setKernel)
            .def("set_time_step", &sim::setTimeStep)
            .def("run_scheme_readdy", [](sim &self, bool defaults) {
                     return std::make_unique<readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>>(
                             self.runScheme<readdy::api::ReaDDyScheme>(defaults)
                     );
                 }
            )
            .def("run", [](sim &self, const readdy::model::observables::time_step_type steps, const double timeStep) {
                py::gil_scoped_release release;
                self.run(steps, timeStep);
            });

    py::class_<kp, std::unique_ptr<kp, readdy::util::nodelete>>(api, "KernelProvider")
            .def_static("get", &kp::getInstance, rvp::reference)
            .def("load_from_dir", &kp::loadKernelsFromDirectory);

    py::class_<pot2>(api, "Pot2")
            .def(py::init<std::string, py::object, py::object>())
            .def("calc_energy", &pot2::calculateEnergy)
            .def("calc_force", &pot2::calculateForce);

    py::class_<kern>(api, "Kernel").def("get_name", &kern::getName, rvp::reference);

    return api.ptr();

}