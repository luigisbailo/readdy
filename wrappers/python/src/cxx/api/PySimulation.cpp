/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <readdy/api/Simulation.h>
#include <readdy/io/File.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/common/nodelete.h>
#include "ExportSchemeApi.h"
#include "PyPotential.h"
#include "PyFunction.h"

namespace py = pybind11;

using rvp = py::return_value_policy;
using sim = readdy::Simulation;
using obs_handle_t = readdy::ObservableHandle;
using kp = readdy::plugin::KernelProvider;
using vec = readdy::model::Vec3;
using pot2 = readdy::rpy::PotentialOrder2Wrapper;
using model = readdy::model::KernelStateModel;
using ctx = readdy::model::KernelContext;
using kern = readdy::model::Kernel;

void exportTopologies(py::module&);

// thin wrappers
void setBoxSize(sim &self, const vec &size) { /* explicitly choose void(vec) signature */ self.setBoxSize(size); }

std::string getSelectedKernelType(sim &self) { /* discard const reference */ return self.getSelectedKernelType(); }

void addParticle(sim &self, const std::string &type, const vec &pos) { self.addParticle(pos[0], pos[1], pos[2], type); }

void registerPotentialOrder2(sim &self, pot2 *potential) {
    self.registerPotentialOrder2(potential);
}

obs_handle_t
registerObservable_Positions(sim &self, unsigned int stride, pybind11::object callbackFun,
                             std::vector<std::string> types) {
    if(callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::Positions>(stride, types);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Positions::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::Positions>(std::move(pyFun), stride, types);
    }
}

obs_handle_t registerObservable_Particles(sim &self, unsigned int stride, pybind11::object callbackFun) {
    if(callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::Particles>(stride);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Particles::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::Particles>(std::move(pyFun), stride);
    }
}


obs_handle_t
registerObservable_RadialDistribution(sim &self, unsigned int stride, pybind11::object callbackFun,
                                      py::array_t<double> &binBorders, std::vector<std::string> typeCountFrom,
                                      std::vector<std::string> typeCountTo, double particleToDensity) {
    const auto info = binBorders.request();
    std::vector<double> binBordersVec{};
    binBordersVec.reserve(info.shape[0]);
    const auto data = static_cast<double *>(info.ptr);
    for (auto i = 0; i < info.shape[0]; ++i) binBordersVec.push_back(data[i]);
    if(callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::RadialDistribution>(
                stride, binBordersVec, typeCountFrom, typeCountTo, particleToDensity
        );
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::RadialDistribution::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::RadialDistribution>(
                std::move(pyFun), stride, binBordersVec, typeCountFrom, typeCountTo, particleToDensity
        );
    }
}

obs_handle_t
registerObservable_CenterOfMass(sim &self, unsigned int stride, const pybind11::object &callbackFun,
                                std::vector<std::string> types) {
    if(callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::CenterOfMass>(stride, types);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::CenterOfMass::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::CenterOfMass>(
                std::move(pyFun), stride, types
        );
    }
}

obs_handle_t
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
    if(callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::HistogramAlongAxis>(
                stride, binBordersVec, types, axis
        );
    } else {
        auto f = readdy::rpy::PyFunction<void(readdy::model::observables::HistogramAlongAxis::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::HistogramAlongAxis>(
                std::move(f), stride, binBordersVec, types, axis
        );
    }
}

obs_handle_t
registerObservable_NParticles(sim &self, unsigned int stride, const py::object &callbackFun,
                              std::vector<std::string> types) {
    if(callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::NParticles>(stride, types);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::NParticles::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::NParticles>(std::move(pyFun), stride, types);
    }
}

obs_handle_t registerObservable_ForcesObservable(sim &self, unsigned int stride, py::object callbackFun,
                                                 std::vector<std::string> types) {
    if(callbackFun.is_none()) {
        return self.registerObservable<readdy::model::observables::Forces>(stride, types);
    } else {
        auto pyFun = readdy::rpy::PyFunction<void(readdy::model::observables::Forces::result_t)>(callbackFun);
        return self.registerObservable<readdy::model::observables::Forces>(std::move(pyFun), stride, types);
    }
}

enum class ParticleTypeFlavor {
    NORMAL = 0, TOPOLOGY, MEMBRANE
};

// module
PYBIND11_PLUGIN (api) {

    py::module api("api", "ReaDDy c++-api python module");

    exportReaDDySchemeApi<readdy::api::ReaDDyScheme>(api, "ReaDDyScheme");
    exportAdvancedSchemeApi<readdy::api::AdvancedScheme>(api, "AdvancedScheme");

    auto topologyModule = api.def_submodule("top");
    exportTopologies(topologyModule);

    py::class_<obs_handle_t>(api, "ObservableHandle")
            .def(py::init<>())
            .def("enable_write_to_file", &obs_handle_t::enableWriteToFile)
            .def("flush", &obs_handle_t::flush)
            .def("__repr__", [](const obs_handle_t& self) {
                return "ObservableHandle(id=" + std::to_string(self.getId()) + ")";
            });

    py::enum_<ParticleTypeFlavor>(api, "ParticleTypeFlavor")
            .value("NORMAL", ParticleTypeFlavor::NORMAL)
            .value("TOPOLOGY", ParticleTypeFlavor::TOPOLOGY)
            .value("MEMBRANE", ParticleTypeFlavor::MEMBRANE);

    py::class_<sim>(api, "Simulation")
            .def(py::init<>())
            .def_property("kbt", &sim::getKBT, &sim::setKBT)
            .def_property("periodic_boundary", &sim::getPeriodicBoundary, &sim::setPeriodicBoundary)
            .def_property("box_size", &sim::getBoxSize, &setBoxSize)
            .def("register_particle_type", [](sim& self, const std::string &name, double diffusionCoefficient, double radius, ParticleTypeFlavor flavor) {
                readdy::model::Particle::flavor_t f = [=] {
                    switch(flavor) {
                        case ParticleTypeFlavor::NORMAL: return readdy::model::Particle::FLAVOR_NORMAL;
                        case ParticleTypeFlavor::MEMBRANE: return readdy::model::Particle::FLAVOR_MEMBRANE;
                        case ParticleTypeFlavor::TOPOLOGY: return readdy::model::Particle::FLAVOR_TOPOLOGY;
                    }
                }();
                return self.registerParticleType(name, diffusionCoefficient, radius);
            }, py::arg("name"), py::arg("diffusion_coefficient"), py::arg("radius"), py::arg("flavor") = ParticleTypeFlavor::NORMAL)
            .def("add_particle", [](sim &self, const std::string &type, const vec &pos) {
                self.addParticle(pos[0], pos[1], pos[2], type);
            })
            .def("is_kernel_selected", &sim::isKernelSelected)
            .def("get_selected_kernel_type", &getSelectedKernelType)
            .def("record_trajectory", &sim::recordTrajectory)
            .def("close_trajectory_file", &sim::closeTrajectoryFile)
            .def("register_potential_order_2", &registerPotentialOrder2)
            .def("register_potential_harmonic_repulsion", &sim::registerHarmonicRepulsionPotential)
            .def("register_potential_piecewise_weak_interaction",
                 &sim::registerWeakInteractionPiecewiseHarmonicPotential)
            .def("register_potential_box", &sim::registerBoxPotential)
            .def("register_potential_sphere", &sim::registerSpherePotential)
            .def("get_particle_positions", &sim::getParticlePositions)
            .def("register_observable_particle_positions", &registerObservable_Positions)
            .def("register_observable_particles", &registerObservable_Particles)
            .def("register_observable_radial_distribution", &registerObservable_RadialDistribution)
            .def("register_observable_histogram_along_axis", &registerObservable_HistogramAlongAxisObservable)
            .def("register_observable_center_of_mass", &registerObservable_CenterOfMass)
            .def("register_observable_n_particles", &registerObservable_NParticles)
            .def("register_observable_forces", &registerObservable_ForcesObservable)
            .def("deregister_observable", [](sim& self, const obs_handle_t& handle) {
                self.deregisterObservable(handle.getId());
            })
            .def("register_reaction_conversion", &sim::registerConversionReaction, rvp::reference_internal)
            .def("register_reaction_enzymatic", &sim::registerEnzymaticReaction, rvp::reference_internal)
            .def("register_reaction_fission", &sim::registerFissionReaction, rvp::reference_internal)
            .def("register_reaction_fusion", &sim::registerFusionReaction, rvp::reference_internal)
            .def("register_reaction_decay", &sim::registerDecayReaction, rvp::reference_internal)
            .def("register_compartment_sphere", &sim::registerCompartmentSphere)
            .def("register_compartment_plane", &sim::registerCompartmentPlane)
            .def("get_recommended_time_step", &sim::getRecommendedTimeStep)
            .def("kernel_supports_topologies", &sim::kernelSupportsTopologies)
            .def("create_topology_particle", &sim::createTopologyParticle)
            .def("add_topology", &sim::addTopology, rvp::reference)
            //.def("add_topology", &sim::addTopology, rvp::reference)
            .def("set_kernel", &sim::setKernel)
            .def("run_scheme_readdy", [](sim &self, bool defaults) {
                     return std::make_unique<readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>>(
                             self.runScheme<readdy::api::ReaDDyScheme>(defaults)
                     );
                 }
            )
            .def("run_scheme_advanced", [](sim &self, bool defaults) {
                return std::make_unique<readdy::api::AdvancedSchemeConfigurator<readdy::api::AdvancedScheme>>(
                        self.runAdvancedScheme<readdy::api::AdvancedScheme>(defaults)
                );
            })
            .def("run", [](sim &self, const readdy::model::observables::time_step_type steps, const double timeStep) {
                py::gil_scoped_release release;
                self.run(steps, timeStep);
            });

    py::class_<kp, std::unique_ptr<kp, readdy::util::nodelete>>(api, "KernelProvider")
            .def_static("get", &kp::getInstance, rvp::reference)
            .def("load_from_dir", &kp::loadKernelsFromDirectory);

    py::class_<pot2>(api, "Pot2")
            .def(py::init<std::string, std::string, py::object, py::object>())
            .def("calc_energy", &pot2::calculateEnergy)
            .def("calc_force", &pot2::calculateForce);

    py::class_<kern>(api, "Kernel").def("get_name", &kern::getName, rvp::reference);

    return api.ptr();

}
