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
 * @file MDGFRDBurst.cpp
 * @brief << brief description >>
 * @author luigisbailo
 * @date 21.02.18
 */

#include "readdy/kernel/mdgfrd/actions/MDGFRDBurst.h"

namespace readdy {
namespace kernel {
namespace mdgfrd {
namespace actions {

void MDGFRDBurst::perform (const util::PerformanceNode &node) {

    auto t = node.timeit();
    auto &stateModel = kernel->getMDGFRDKernelStateModel();
    const auto &context = kernel->context();
    auto data = stateModel.getParticleData();
    const auto size = data->size();
    const MDGFRDStateModel::neighbor_list *nl = stateModel.getNeighborList();
    model::potentials::PotentialRegistry::potential_o2_registry pot2;

    using iter_t = data::EntryDataContainer::iterator;

    auto worker = [&context,data,nl](std::size_t, std::size_t beginIdx, iter_t entry_begin, iter_t entry_end){

        model::Context::shortest_dist_fun d = context.shortestDifferenceFun();
        std::size_t  idx = beginIdx;
        for (auto it = entry_begin; it != entry_end; ++it, ++idx) {
            if (it->GFpropagation){
                scalar Rint = 1;
                auto domainSize = it->domainSize;
                nl->forEachNeighbor(idx, [&](auto neighborIndex)){
                    auto &neighbor = data->entry_at(neighborIndex);
                    if (!neighbor.deactivated && !neighbor.GFpropagation){
                        auto x_ij = d(it->pos,neighbor.pos);
                        auto rInteraction = pot2.find(std::tie(it.type, neighbor.type))->getCutoffRadius();
                        if (x_ij-domainSize-rInteraction<neighbor.minDomainSize){
                            it->GFpropagation = false;
                            auto deltaT = it.simulationTime - it.constructionTime;
                            auto radius = domainSize / 2;
                            Vec3 displacement.polarTransform(radius);
                            data->displace(idx,displacement);
                            it->domainSize = 0;
                        }
                    }
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