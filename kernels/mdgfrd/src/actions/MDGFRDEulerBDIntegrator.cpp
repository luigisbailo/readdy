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
 * @file MDGFRDEulerBDIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 15.02.18
 */

#include <readdy/kernel/mdgfrd/actions/MDGFRDEulerBDIntegrator.h>

namespace readdy {
namespace kernel {
namespace mdgfrd {
namespace actions {

namespace rnd = readdy::model::rnd;

void MDGFRDEulerBDIntegrator::perform(const readdy::util::PerformanceNode &node) {
    auto t = node.timeit();
    auto data = kernel->getMDGFRDKernelStateModel().getParticleData();
    const auto size = data->size();

    const auto &context = kernel->context();
    using iter_t = data::EntryDataContainer::iterator;

    const auto dt = timeStep;

    auto worker = [&context, data, dt](std::size_t, std::size_t beginIdx, iter_t entry_begin, iter_t entry_end)  {
        const auto &fixPos = context.fixPositionFun();
        const auto kbt = context.kBT();
        std::size_t idx = beginIdx;
        for (auto it = entry_begin; it != entry_end; ++it, ++idx) {
            it->simulationTime += dt;
            if(!it->deactivated ) {
                const scalar D = context.particle_types().diffusionConstantOf(it->type);
                if ( it->GFpropagation && it->simulationTime > it->exitTime) {
                    it->GFpropagation = false;
                    const auto dt = it->simulationTime - it->exitTime;
                    auto randomDisplacement = std::sqrt(2. * D * dt) * rnd::normal3<readdy::scalar>(0, 1);
                    auto rnd1 = readdy::model::rnd::uniform_real();
                    auto rnd2 = readdy::model::rnd::uniform_real();
                    Vec3 domainExit;
                    domainExit.polarTransform(it->domainSize, rnd1, rnd2);
                    data->displace(idx, randomDisplacement + domainExit);

                }
                else{
                    auto randomDisplacement = std::sqrt(2. * D * dt) * rnd::normal3<readdy::scalar>(0, 1);
                    auto deterministicDisplacement = it->force * dt * D / kbt;
                    data->displace(idx, randomDisplacement + deterministicDisplacement);

                }

            }
        }
    };

    std::vector<util::thread::joining_future<void>> waitingFutures;
    waitingFutures.reserve(kernel->getNThreads());
    auto &pool  = kernel->pool();
    {
        auto it = data->begin();
        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(kernel->getNThreads());

        auto granularity = kernel->getNThreads();
        const std::size_t grainSize = size / granularity;

        std::size_t idx = 0;
        for (auto i = 0_z; i < granularity-1; ++i) {
            auto itNext = it + grainSize;
            if(it != itNext) {
                waitingFutures.emplace_back(pool.push(worker, idx, it, itNext));
            }
            it = itNext;
            idx += grainSize;
        }
        if(it != data->end()) {
            waitingFutures.emplace_back(pool.push(worker, idx, it, data->end()));
        }
    }

}

MDGFRDEulerBDIntegrator::MDGFRDEulerBDIntegrator(MDGFRDKernel *kernel, scalar timeStep)
        : readdy::model::actions::EulerBDIntegrator(timeStep), kernel(kernel), init(true) {}

}
}
}
}
