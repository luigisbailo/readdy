/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/

/**
 * << detailed description >>
 *
 * @file SingleCPUObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */


#include <readdy/kernel/singlecpu/SCPUKernel.h>

#include <readdy/kernel/singlecpu/observables/SCPUObservableFactory.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUAggregators.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace observables {
SCPUObservableFactory::SCPUObservableFactory(readdy::kernel::scpu::SCPUKernel *const kernel)
        : ObservableFactory(kernel), kernel(kernel) {
}

std::unique_ptr<readdy::model::observables::HistogramAlongAxis>
SCPUObservableFactory::histogramAlongAxis(stride_type stride, std::vector<scalar> binBorders,
                                          std::vector<std::string> typesToCount, unsigned int axis) const {
    return {std::make_unique<SCPUHistogramAlongAxis>(kernel, stride, binBorders, typesToCount, axis)};
}

std::unique_ptr<readdy::model::observables::NParticles>
SCPUObservableFactory::nParticles(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<SCPUNParticles>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::Forces>
SCPUObservableFactory::forces(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<SCPUForces>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::Positions>
SCPUObservableFactory::positions(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<SCPUPositions>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::RadialDistribution>
SCPUObservableFactory::radialDistribution(stride_type stride, std::vector<scalar> binBorders,
                                          std::vector<std::string> typeCountFrom, std::vector<std::string> typeCountTo,
                                          scalar particleDensity) const {
    return {std::make_unique<readdy::model::observables::RadialDistribution>(kernel, stride, binBorders, typeCountFrom,
                                                                             typeCountTo, particleDensity)};
}

std::unique_ptr<readdy::model::observables::Particles>
SCPUObservableFactory::particles(stride_type stride) const {
    return {std::make_unique<SCPUParticles>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::MeanSquaredDisplacement>
SCPUObservableFactory::msd(stride_type stride, std::vector<std::string> typesToCount,
                           readdy::model::observables::Particles *particlesObservable) const {
    return {std::make_unique<SCPUMeanSquaredDisplacement<>>(kernel, stride, typesToCount, particlesObservable)};
}

std::unique_ptr<readdy::model::observables::Reactions>
SCPUObservableFactory::reactions(readdy::model::observables::ObservableFactory::stride_type stride) const {
    return {std::make_unique<SCPUReactions>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::ReactionCounts>
SCPUObservableFactory::reactionCounts(readdy::model::observables::ObservableFactory::stride_type stride) const {
    return {std::make_unique<SCPUReactionCounts>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::Virial>
SCPUObservableFactory::virial(readdy::model::observables::ObservableFactory::stride_type stride) const {
    return {std::make_unique<SCPUVirial>(kernel, stride)};
}

}
}
}
}