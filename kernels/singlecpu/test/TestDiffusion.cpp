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
 * @file TestDiffusion.cpp
 * @brief << brief description >>
 * @author luigisbailo
 * @date 19.01.18
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/actions/Actions.h>
#include <readdy/api/Simulation.h>

//using namespace readdy;

TEST(SingleCPUTestDiffusion, SanityTest) {
    EXPECT_TRUE(1==1);
    std::cout << std::endl << std::endl << std::endl;

    readdy::Simulation simulation;
    simulation.setKernel("MDGFRD");
    unsigned int n_particles = 3;
    simulation.setBoxSize(10,10,10);
    readdy::scalar diffusionConstant = 0.01;
    simulation.registerParticleType("type", diffusionConstant);
    for (auto _ = 0; _ < n_particles; ++_) {
        simulation.addParticle("type", 0, 0, 0);
    }
    readdy::scalar timestep = 1;
    readdy::scalar n_steps = 5;

    std::cout << "hello world!" << std::endl;

//    int n_callbacks = 0;
//    simulation.registerObservable(simulation.observe().positions(1),
//                                  [&n_callbacks](const readdy::model::observables::Positions::result_type &result) {
//                                      ++n_callbacks;
//                                      EXPECT_EQ(10, result.size());
//                                  });
    simulation.run(n_steps, timestep);
//    EXPECT_EQ(101, n_callbacks);

    std::cout << simulation.getKBT() << std::endl;
    std::cout << simulation.getBoxSize() << std::endl;
    std::vector<readdy::Vec3 > Vpos = simulation.getAllParticlePositions();

    for (int count=0; count<n_particles; count++){
        std::cout << Vpos[count] << std::endl;
    }

//    auto conf = simulation.runScheme();
//    conf.withIntegrator("MFPIntegrator");
//    auto scheme = conf.configure(0.1);
//    scheme->run(1000);

    std::cout << std::endl << std::endl << std::endl;

}