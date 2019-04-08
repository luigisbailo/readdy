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
 * @file TestSimulationLoop.cpp
 * @author clonker
 * @author chrisfro**
 * @date 23.08.16
 */


#include <catch2/catch.hpp>

#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

#include <readdy/plugin/KernelProvider.h>
#include <readdy/api/Simulation.h>

namespace api = readdy::api;

using namespace readdytesting::kernel;


TEMPLATE_TEST_CASE("Test simulation loop", "[loop]", SingleCPU, CPU) {
    readdy::Simulation simulation {create<TestType>()};

    SECTION("Correct number of timesteps") {
        unsigned int counter = 0;
        auto increment = [&counter](readdy::model::observables::NParticles::result_type result) {
            counter++;
        };
        auto obsHandle = simulation.registerObservable(simulation.observe().nParticles(1), increment);
        simulation.run(3, 0.1);
        REQUIRE(counter == 4);
    }
    SECTION("Simple stopping criterion") {
        unsigned int counter = 0;
        auto increment = [&counter](readdy::model::observables::NParticles::result_type result) {
            counter++;
        };
        auto obsHandle = simulation.registerObservable(simulation.observe().nParticles(1), increment);
        auto shallContinue = [](readdy::time_step_type currentStep) {
            return currentStep < 5;
        };
        auto loop = simulation.createLoop(.1);
        loop.run(shallContinue);
        REQUIRE(counter == 6);
    }
    SECTION("Complex stopping criterion") {
        simulation.context().particleTypes().add("A", 0.);
        // A -> A + A, with probability = 1 each timestep. After 3 timesteps there will be 8 particles.
        // The counter will be 4 by then.
        simulation.context().reactions().addFission("bla", "A", "A", "A", 1e8, 0.);
        simulation.addParticle("A", 0, 0, 0);
        unsigned int counter = 0;
        bool doStop = false;
        auto increment = [&counter, &doStop](readdy::model::observables::NParticles::result_type result) {
            counter++;
            if (result[0] >= 8) {
                doStop = true;
            }
        };
        auto obsHandle = simulation.registerObservable(simulation.observe().nParticles(1), increment);
        auto shallContinue = [&doStop](readdy::time_step_type currentStep) {
            return !doStop;
        };
        simulation.createLoop(1.).run(shallContinue);
        REQUIRE(counter == 4);
    }
    SECTION("Skin size sanity check") {
        simulation.context().particleTypes().add("A", 1.);
        simulation.context().boxSize() = {{10., 10., 10.}};
        simulation.context().periodicBoundaryConditions() = {{true, true, true}};
        simulation.context().potentials().addHarmonicRepulsion("A", "A", 1., 2.);
        simulation.addParticle("A", 0., 0., 0.);
        simulation.addParticle("A", 1.5, 0., 0.);
        auto loop = simulation.createLoop(.001);
        loop.neighborListCutoff() += 0.1; // adding a skin/padding
        loop.run(10);
    }
}
