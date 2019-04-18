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
    scalar domainReduction = 0;
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
            speciesDataI.minDomainSize = alpha * std::sqrt(infoI.diffusionConstant * timeStep);
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
            speciesDataI.domainReduction = reduction_coefficient * std::sqrt(2 * infoI.diffusionConstant * timeStep);
        }
    };

    void perform() override {
        // @todo set up neighborlist again with skin=maxcutoff
        if (firstPerform) {
            initialize();
        }
        // todo bookkeep particleData (because particles might have disappeared or appeared)

        // fractional propagation for previously protected particles
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
                                scalar deltaT = stateModel.time() - particleData.at(entry.type).constructionTime;
                                scalar xi = readdy::model::rnd::uniform_real();
                                scalar b = particleData.at(entry.id).domainSize;
                                scalar iDiff = context.particleTypes().diffusionConstantOf(entry.type);
                                scalar radius = drawPosNewt(deltaT,b,iDiff,xi);
                                Vec3 displacement = randomOrientation(radius);
                                entry.pos += displacement;
                                bcs::fixPosition(entry.pos, box, pbc);
                                particleData.at(entry.id).GFpropagation = false;
                                particleData.at(entry.id).domainSize = 0;
                                particleData.at(entry.id).exitTime = stateModel.time();
                                // next line is optional. A not GF particle shold nevet access the exit position
                                particleData.at(entry.id).exitPosition = entry.pos;

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
                    auto iDiff = context.particleTypes().diffusionConstantOf(entry.type);

                    nl->forEachNeighbor(it, cell, [&](const std::size_t neighbor) {

                        auto &&neighborEntry = pd->entry_at(neighbor);
                        auto jDiff = context.particleTypes().diffusionConstantOf(neighborEntry.type);
                        auto dij = bcs::dist(entry.pos, neighborEntry.pos, box, pbc)
                                   - maxInteractionRadius.at(std::make_pair(entry.type,neighborEntry.type));

                        if ( dij < 0 ){
                            // particles are interacting - domain is not constructed
                            //@todo break this loop
                        }
                        scalar tempDomainSize;
                        if ( !particleData.at(neighborEntry.id).GFpropagation ){

                            tempDomainSize = dij / ( 1 + std::sqrt(iDiff/jDiff) );

                        }
                        else {

                            auto dijExit = bcs::dist(entry.pos, particleData.at(neighborEntry.id).exitPosition, box, pbc);
                            auto neighbDomainSize = particleData.at(neighborEntry.id).domainSize;
                            auto interactionR = maxInteractionRadius.at(std::make_pair(entry.type,neighborEntry.type));

                            // distance to the domain borders of the neighbor particle minus interaction
                            // this is the largest value the domain can assume
                            auto distCurrentDomain = dij - neighbDomainSize - interactionR;

                            // distance to the exit position of the neighbor particle minus interaction
                            auto distNextDomain = dijExit - neighbDomainSize - interactionR;

                            auto deltaT = stateModel.time() - particleData.at(entry.type).constructionTime;

                            if ( particleData.at(neighborEntry.id).exitTime - particleData.at(neighborEntry.id).exitTime >
                                distNextDomain*distNextDomain/6/context.particleTypes().diffusionConstantOf(entry.type) ){
                                // this conditions ensures that 'deltaRoot' (later defined) is larger than zero
                                // the largest value is taken
                                tempDomainSize = distCurrentDomain;

                            }
                            else if( entry.type == neighborEntry.type) {

                                tempDomainSize = distNextDomain/2 + (3*iDiff*(deltaT)) / deltaT;

                            }
                            else {

                                auto deltaRoot = (distNextDomain*distNextDomain)/(36*iDiff*jDiff) -
                                        ((1/(6*jDiff) - (6*iDiff)) * (distNextDomain*distNextDomain/(6*jDiff)) + deltaT);

                                tempDomainSize = (distNextDomain/(6*jDiff) - std::sqrt(deltaRoot) ) /
                                        (1/(6*jDiff) - 1/(6*iDiff) );

                            };

                            // distCurrentDomain is the largest value that can be given to the domain
                            tempDomainSize = std::min(tempDomainSize,distCurrentDomain);
                        }

                        //The domain size is furhter reduced to decrease bursts
                        tempDomainSize -= speciesData.at(entry.type).domainReduction;

                        if (tempDomainSize<domainSize){
                            domainSize = tempDomainSize;
                        }

                        if ( domainSize < speciesData.at(entry.type).minDomainSize){
                            //@todo break the loop
                        }

                    });

                    if (domainSize>speciesData.at(entry.type).minDomainSize){
                        // domain is constructed

                        auto xi = readdy::model::rnd::uniform_real();
                        auto tau = drawTimeNewt(domainSize,context.particleTypes().diffusionConstantOf(entry.type),xi);
                        particleData.at(entry.id).GFpropagation = true;
                        particleData.at(entry.id).domainSize = domainSize;
                        particleData.at(entry.id).exitTime = stateModel.time() + tau;
                        particleData.at(entry.id).constructionTime = stateModel.time();
                        Vec3 displacement = randomOrientation(domainSize);
                        particleData.at(entry.id).exitPosition = entry.pos + displacement;

                    }
                    else{

                        const scalar D = context.particleTypes().diffusionConstantOf(entry.type);
                        const auto randomDisplacement = std::sqrt(2. * D * _timeStep) *
                                                        (readdy::model::rnd::normal3<readdy::scalar>());
                        entry.pos += randomDisplacement;
                        //@todo forces not updated - optimally forces should be computed only here
                        const auto deterministicDisplacement = entry.force * _timeStep * D / kbt;
                        entry.pos += deterministicDisplacement;
                        bcs::fixPosition(entry.pos, box, pbc);

                    }
                }
            }
        }

    }

private:

    //alpha sets the minimum domain size - large alpha values makes more difficult to construct domains
    scalar alpha = 8;
    // reduction_coefficient sets how much a domain is reduced
    scalar reduction_coefficient = 4;
    // convergence for PDF samplers
    scalar draw_convergence = 0.001;
    // max number of iterations in the samplers
    int max_iterations = 100;

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

    double drawTimeNewt ( double b, double D, double xi ) {
        // samples the exit time for the following particle parameters
        // b - domain size
        // D - diffusion coefficient
        // xi - random number [0,1]

        double t,tmem;

        t = 0.0917517*b*b/D;
        tmem = -1;

        double S = 1-Sfunct (t,b,D);
        double dS = Sder (t,b,D);
        int count = 0;

        while ( fabs(t-tmem) > draw_convergence | fabs(S-xi) > draw_convergence ) {


            count++;
            if (count > max_iterations){

                return t;

            }

            tmem = t;
            t = t + (S-xi)/dS;

            S = 1-Sfunct (t,b,D);
            if ( S==1 ) return t;
            dS = Sder (t,b,D);

        }

        return t;

    }


    double drawPosNewt ( double t, double b, double D, double xi ) {
        // samples the position at bursts for the following particle parameters
        // t - sampling time
        // b - domain size
        // D - diffusion coefficient
        // xi - random number [0,1]

        double r;

        double t0=0.063;
        double t1=0.234;
        if (t<t0*b*b/D) r = sqrt(t*D)*2;
        else if (t<t1*b*b/D) {
            double R0=2*sqrt(t0*b*b);
            double R1=0.646*b;
            double beta = (R0*exp(pow(t0,0.5)+pow(t1,0.5))-R1*exp(pow(t0,0.5)+pow(t1,0.5)))/(exp(pow(t0,0.5))-exp(pow(t1,0.5)));
            double gamma = -(R0*(exp(pow(t1,0.5))-1)*exp(pow(t0,0.5))-R1*(exp(pow(t0,0.5))-1)*exp(pow(t1,0.5)))
                           / (exp(pow(t0,0.5))-exp(pow(t1,0.5)));
            r=beta*(1-exp(-pow(t*D/b/b,0.5)))+gamma;
        }
        else r = 0.646*b;

        double S = Sfunct (t,b,D);
        double P = Pfunct (r,t,b,D,S);
        double dP = Pder (r,t,b,D,S);
        int count = 0;
        double rMem=-1;

        while ( fabs(r-rMem) > draw_convergence | fabs(P-xi) > draw_convergence ) {

            count++;
            if (count > max_iterations){

                return r;

            }

            rMem = r;
            r = r - (P-xi)/dP;

            S = Sfunct (t,b,D);
            P = Pfunct (r,t,b,D,S);
            dP = Pder (r,t,b,D,S);

        }

        return r;

    }


    double Sder ( double t, double b, double D) {
        // probability distribution of the first exit-time of the particle
        // t - sampling time
        // b - domain size
        // D - diffusion coefficient

        double coeff1 = exp ( - M_PI*M_PI*D*t/(b*b) );
        double coeff2 = 2*D*M_PI*M_PI/b/b;
        double S;
        double term;
        double termA,termB,termC,termD;
        int m;

        double conv = 0.00001/(b*b/D);

        S = 0;
        m = 1;

        if ( t>= b*b/D/100 ) {

            do {

                termA = pow (coeff1,m*m);
                termB = m*m;
                termC = pow (coeff1,(m+1)*(m+1));
                termD = (m+1)*(m+1);
                term = coeff2 * ( termA*termB - termC*termD );

                S += term;
                m += 2;

            } while ( fabs (term) > conv | termC > termD/100 | S < 0 );
            //The second condition takes into account that the series, for short times, can be negative in the beginning,
            // because of the term m*m
            //The series, when the term m*m  is dominant at short times, grows to a peak, and then finally deceases to zero
            //The S<0 condition is to ensure a longer convergence for short times
        }
        else
            S = 0;

        return S;

    }

    double Sfunct ( double t, double b, double D) {
        // cumulative of above function

        double coeff1 = exp ( -M_PI*M_PI*D*t/(b*b) );
        double coeff2 = 2;
        double S = 0;
        double term,termA,termB;
        int m = 1;
        double conv = 0.0000001/(b*b/D);

        if ( t>= b*b/D/100 ) {

            do {

                termA = pow (coeff1,m*m);
                termB = pow (coeff1,(m+1)*(m+1));
                term = coeff2 * (termA-termB);


                S += term;
                m += 2;

            } while ( fabs (term) > conv );

            S = 1-S;
        }
        else
            S = 0;

        return S;

    }

    double Pder ( double radius, double t, double b, double D, double S ) {
        // probability distribution of particle position inside a domain, assuming the particle
        // has never escaped the domain before
        // radius - distance from the center of the domain
        // t - sampling time
        // b - domain size
        // D - diffusion coefficient


        double coeff1 = exp ( -M_PI*M_PI*D*t/(b*b)  );
        double coeff2 = 2*M_PI/b/b;
        double P;
        double term;
        double termA,termB;
        int m;
        double conv = 0.0001/b;

        P = 0;
        m = 1;

        do {

            termA = pow(coeff1,m*m);
            termB = sin (m*M_PI*radius/b) * m * radius;
            term = coeff2 * termA * termB / (1-S);

            P += term;
            m += 1;

        } while ( fabs(term) > conv  |  termA>termB/100 );

        return P;
    }



    double Pfunct ( double radius, double t, double b, double D, double S ) {
        // cumulative of above probability distribution

        double coeff1 = exp ( - M_PI*M_PI*D*t/(b*b));
        double coeff2 = 2/(b*M_PI);
        double P;
        double  term, termA, termB;
        int m;
        double conv = 0.00000001/b;

        P = 0;
        m = 1;

        do {

            termA =  pow(coeff1,m*m);
            termB =  b*sin(m*M_PI*radius/b)/m - M_PI*radius*cos(m*M_PI*radius/b);
            term =  coeff2 * termA * termB / (1-S);
            P += term;
            m += 1;

        } while ( fabs(term) > conv  |  termA>termB/100 );

        return P;

    }



};

}
