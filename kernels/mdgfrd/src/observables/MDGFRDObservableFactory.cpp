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


/**
 * << detailed description >>
 *
 * @file ObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 15.02.18
 */

#include <readdy/kernel/mdgfrd/observables/MDGFRDObservableFactory.h>
#include <readdy/kernel/mdgfrd/MDGFRDKernel.h>
#include <readdy/kernel/mdgfrd/observables/MDGFRDObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUAggregators.h>

namespace readdy {
namespace kernel {
namespace mdgfrd {
namespace observables {
MDGFRDObservableFactory::MDGFRDObservableFactory(MDGFRDKernel *const kernel) : readdy::model::observables::ObservableFactory(kernel),
                                                                kernel(kernel) {
}

std::unique_ptr<model::observables::HistogramAlongAxis>
MDGFRDObservableFactory::histogramAlongAxis(stride_type stride, std::vector<scalar> binBorders,
                                         std::vector<std::string> typesToCount, unsigned int axis) const {
    return {std::make_unique<MDGFRDHistogramAlongAxis>(kernel, stride, binBorders, typesToCount, axis)};
}

std::unique_ptr<model::observables::NParticles>
MDGFRDObservableFactory::nParticles(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<MDGFRDNParticles>(kernel, stride, typesToCount)};
}

std::unique_ptr<model::observables::Forces>
MDGFRDObservableFactory::forces(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<MDGFRDForces>(kernel, stride, typesToCount)};
}

std::unique_ptr<model::observables::Positions>
MDGFRDObservableFactory::positions(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<MDGFRDPositions>(kernel, stride, typesToCount)};
}

std::unique_ptr<model::observables::RadialDistribution>
MDGFRDObservableFactory::radialDistribution(stride_type stride, std::vector<scalar> binBorders,
                                         std::vector<std::string> typeCountFrom, std::vector<std::string> typeCountTo,
                                         scalar particleDensity) const {
    return {std::make_unique<readdy::kernel::scpu::observables::SCPURadialDistribution<MDGFRDKernel>>(
            kernel, stride, binBorders, typeCountFrom, typeCountTo, particleDensity
    )};
}

std::unique_ptr<model::observables::Particles> MDGFRDObservableFactory::particles(stride_type stride) const {
    return {std::make_unique<MDGFRDParticles>(kernel, stride)};
}

std::unique_ptr<model::observables::MeanSquaredDisplacement>
MDGFRDObservableFactory::msd(stride_type stride, std::vector<std::string> typesToCount,
                          model::observables::Particles *particlesObservable) const {
    return {std::make_unique<readdy::kernel::scpu::observables::SCPUMeanSquaredDisplacement<MDGFRDKernel>>(
            kernel, stride, typesToCount, particlesObservable
    )};
}

std::unique_ptr<model::observables::Reactions> MDGFRDObservableFactory::reactions(stride_type stride) const {
    return {std::make_unique<MDGFRDReactions>(kernel, stride)};
}

std::unique_ptr<model::observables::ReactionCounts> MDGFRDObservableFactory::reactionCounts(stride_type stride) const {
    return {std::make_unique<MDGFRDReactionCounts>(kernel, stride)};
}

std::unique_ptr<model::observables::Virial>
MDGFRDObservableFactory::virial(stride_type stride) const {
    return {std::make_unique<MDGFRDVirial>(kernel, stride)};
}

}
}
}
}
