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
 * @todo This is a draft
 *
 * @file SCPUMdgfrdIntegrator.h
 * @brief << brief description >>
 * @author chrisfroe
 * @author luigisbailo
 * @date 12.02.19
 */
#pragma once
#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/common/boundary_condition_operations.h>

namespace readdy::kernel::scpu::actions {

struct MdgfrdParticleData {
    scalar exitTime = 0;
    scalar constructionTime = 0;
    scalar domainSize = 0;
    Vec3 exitPosition  = {0,0,0};
    bool GFpropagation = false;
};

struct MdgfrdSpeciesData {
    scalar minDomainSize = 0;
    scalar maxDomainSize = 0;
};


class SCPUMdgfrdIntegrator : public readdy::model::actions::MdgfrdIntegrator {
public:

    SCPUMdgfrdIntegrator(SCPUKernel *kernel, scalar timeStep)
            : readdy::model::actions::MdgfrdIntegrator(timeStep), kernel(kernel) {
        const auto &context = kernel->context();
        const auto &particlesTypeIds = context.particleTypes().typesFlat();
        // set up species data
        for (const auto& speciesI : particlesTypeIds) {
            const auto &infoI = context.particleTypes().infoOf(speciesI);
            MdgfrdSpeciesData speciesDataI;
            // fixme: 4 is just a guess
            speciesDataI.minDomainSize = 4. * std::sqrt(infoI.diffusionConstant * timeStep);
            scalar maxInteractionRadius = 0.;
            for (const auto& speciesJ : particlesTypeIds) {
                const auto &infoJ = context.particleTypes().infoOf(speciesJ);
                // loop over reaction of i and j
                // loop over potentials of i and j
                scalar interactionDistance = 10.;
                if (interactionDistance > maxInteractionRadius) {
                    maxInteractionRadius = interactionDistance;
                }
            }
            const auto &nlCellSize = kernel->getSCPUKernelStateModel().getNeighborList()->cellSize();
            //const auto minCellSize = std::min(std::min(nlCellSize[0], nlCellSize[1]), nlCellSize[2]);
            const auto minCellSize = *std::min_element(nlCellSize.data.begin(), nlCellSize.data.end());
            speciesDataI.maxDomainSize = (minCellSize - maxInteractionRadius) / 2.;
            speciesData.emplace(std::make_pair(speciesI, speciesDataI));
        }
    };

    void perform() override {
        // @todo set up neighborlist again with skin=maxcutoff
        if (firstPerform) {
            initialize();
        }
        // todo bookkeep particleData (because particles might have disappeared or appeared)

        // fractional propagation for previously protected particles

        // -----
        // update neighbor list
        // update distances
        // ---- the two seps above can be skipped under the assumption that
        //      only rarely a particle bursts a domain after exiting its own domains,
        //      or the gap left between domains should be larger than just the reaction distance

        // burst
        // update neighbor list
        // update distances
        // propagate or construct domain, find out minDomainSize

        const auto &context = kernel->context();
        const auto &pbc = context.periodicBoundaryConditions().data();
        const auto &kbt = context.kBT();
        const auto &box = context.boxSize().data();
        auto& stateModel = kernel->getSCPUKernelStateModel();
        const auto pd = stateModel.getParticleData();
        const auto t = stateModel.time();
        const auto &nl = stateModel.getNeighborList();

        for (auto& entry : *pd) {

        }

        // fractional propagation after exiting the protective domain
        for(auto& entry : *pd) {
            if(!entry.is_deactivated()) {
                if (particleData.at(entry.id).GFpropagation && particleData.at(entry.id).exitTime < stateModel.time()){

                    const scalar D = context.particleTypes().diffusionConstantOf(entry.type);
                    const auto randomDisplacement = std::sqrt(2. * D * (stateModel.time()-particleData.at(entry.id).exitTime )) *
                                                    (readdy::model::rnd::normal3<readdy::scalar>());
                    entry.pos = particleData.at(entry.id).exitPosition + randomDisplacement;
                    bcs::fixPosition(entry.pos, box, pbc);
                    particleData.at(entry.id).GFpropagation = false;
                    particleData.at(entry.id).domainSize = 0;

                }
            }
        }

        if (nl->cutoff() > 0) {
            stateModel.updateNeighborList();
        }

        // bursts - iteration over cells, it might be more performant to iterate over particles
        for (auto cell = 0_z; cell < nl->nCells(); ++cell) {
            for (auto it = nl->particlesBegin(cell); it != nl->particlesEnd(cell); ++it) {
                auto pidx = *it;
                auto &&entry = pd->entry_at(pidx);
                if (particleData.at(entry.id).GFpropagation ) {
                    nl->forEachNeighbor(it, cell, [&](const std::size_t neighbor) {
                        auto &&neighborEntry = pd->entry_at(neighbor);
                        if (!particleData.at(neighborEntry.id).GFpropagation) {
                            auto dij = bcs::dist(entry.pos, neighborEntry.pos, box, pbc);
                            if (dij - maxInteractionRadius.at(std::make_pair(entry.type,neighborEntry.type)) <
                            particleData.at(entry.id).domainSize + speciesData.at(neighborEntry.type).minDomainSize){
                                //entry is burst
                                // @todo - WRONG SAMPLING - the radius should be sample from appropriate PDF
                                scalar radius = readdy::model::rnd::uniform_real() *  particleData.at(entry.id).domainSize;
                                Vec3 displacement = randomOrientation(radius);
                                entry.pos += displacement;
                                bcs::fixPosition(entry.pos, box, pbc);
                                particleData.at(entry.id).GFpropagation = false;
                                particleData.at(entry.id).domainSize = 0;
                                particleData.at(entry.id).exitTime = stateModel.time();

                            }
                        }
                    });
                }
            }
        }

        if (nl->cutoff() > 0) {
            stateModel.updateNeighborList();
            // @todo calc forces
        }

        // propagate or construct domain - it might be more performant to iterate over particles
        for (auto cell = 0_z; cell < nl->nCells(); ++cell) {
            for (auto it = nl->particlesBegin(cell); it != nl->particlesEnd(cell); ++it) {
                auto pidx = *it;
                auto &&entry = pd->entry_at(pidx);
                if (!particleData.at(entry.id).GFpropagation ) {
                    scalar domainSize =  speciesData.at(entry.type).maxDomainSize;
                    nl->forEachNeighbor(it, cell, [&](const std::size_t neighbor) {
                        auto &&neighborEntry = pd->entry_at(neighbor);
                        if (!particleData.at(neighborEntry.id).GFpropagation) {

                            auto dij = bcs::dist(entry.pos, neighborEntry.pos, box, pbc);
                            auto dijExit = bcs::dist(entry.pos, particleData.at(neighborEntry.id).exitPosition, box, pbc);
                            //@todo implement algorithm for domain size
                            scalar tempDomainSize = dij/2 - maxInteractionRadius.at(std::make_pair(entry.type,neighborEntry.type));

                            if (tempDomainSize<domainSize){
                                domainSize = tempDomainSize;
                            }
                        }
                    });

                    if (domainSize>speciesData.at(entry.type).minDomainSize){

//                      @todo the exit-time should be sampled from an appropriate PDF
                        scalar tau = 10 * timeStep();
                        particleData.at(entry.id).GFpropagation = true;
                        particleData.at(entry.id).domainSize = domainSize;
                        particleData.at(entry.id).exitTime = stateModel.time() + tau;
                        Vec3 displacement = randomOrientation(domainSize);
                        particleData.at(entry.id).exitPosition, box, pbc) = entry.pos + displacement;
                    }
                    else{

                        const scalar D = context.particleTypes().diffusionConstantOf(entry.type);
                        const auto randomDisplacement = std::sqrt(2. * D * _timeStep) *
                                                        (readdy::model::rnd::normal3<readdy::scalar>());
                        entry.pos += randomDisplacement;
                        //WARNING: forces not updated
                        const auto deterministicDisplacement = entry.force * _timeStep * D / kbt;
                        entry.pos += deterministicDisplacement;
                        bcs::fixPosition(entry.pos, box, pbc);

                    }
                }
            }
        }

    }

private:
    void initialize() {
        // Todo setup particle data
    }

    SCPUKernel *kernel;

    bool firstPerform = true;


    std::map<readdy::model::Particle::id_type, MdgfrdParticleData> particleData;
    std::unordered_map<readdy::model::Particle::type_type, MdgfrdSpeciesData> speciesData;
    util::particle_type_pair_unordered_map<scalar> maxInteractionRadius;

    Vec3 randomOrientation ( scalar radius ) {

        //Theta is defined in [0,2pi]
        //Phi is defined in [0,pi]
        auto u = readdy::model::rnd::uniform_real();
        auto v = readdy::model::rnd::uniform_real();
        // @todo check if there is internal pi definition
        double theta = 2 * M_PI * u;
        double phi = std::acos( 2*v - 1 );

        Vec3 vector;

        vector[0] = radius * cos(theta) * sin(phi);
        vector[1] = radius * sin(theta) * sin(phi);
        vector[2] = radius * cos(phi);

        return vector;
    }

};

}
