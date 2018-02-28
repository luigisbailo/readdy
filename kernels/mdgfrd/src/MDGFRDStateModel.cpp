/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * @file MDGFRDStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 15/02/18
 */


#include <future>
#include <readdy/kernel/mdgfrd/MDGFRDStateModel.h>

namespace readdy {
namespace kernel {
namespace mdgfrd {

namespace thd = readdy::util::thread;

using entries_it = MDGFRDStateModel::data_type::Entries::iterator;
using topologies_it = std::vector<std::unique_ptr<readdy::model::top::GraphTopology>>::const_iterator;
using pot1Map = readdy::model::potentials::PotentialRegistry::potential_o1_registry;
using pot2Map = readdy::model::potentials::PotentialRegistry::potential_o2_registry;
using dist_fun = readdy::model::Context::shortest_dist_fun;

const std::vector<Vec3> MDGFRDStateModel::getParticlePositions() const {
    const auto data = getParticleData();
    std::vector<Vec3> target{};
    target.reserve(data->size());
    for (const auto &entry : *data) {
        if (!entry.deactivated) target.push_back(entry.pos);
    }
    return target;
}

const std::vector<readdy::model::Particle> MDGFRDStateModel::getParticles() const {
    const auto data = getParticleData();
    std::vector<readdy::model::Particle> result;
    result.reserve(data->size());
    for (const auto &entry : *data) {
        if (!entry.deactivated) {
            result.push_back(data->toParticle(entry));
        }
    }
    return result;
}

MDGFRDStateModel::MDGFRDStateModel(data_type &data, const readdy::model::Context &context,
                             thread_pool &pool)
        : _pool(pool), _context(context), _data(data) {
    _neighborList = std::make_unique<neighbor_list>(_data.get(), _context.get(), _pool.get());
    _reorderConnection = std::make_unique<readdy::signals::scoped_connection>(
            getParticleData()->registerReorderEventListener([this](const std::vector<std::size_t> &indices) -> void {
            }));
}


void MDGFRDStateModel::resetReactionCounts() {
    if(!reactionCounts().empty()) {
        for(auto &e : reactionCounts()) {
            e.second = 0;
        }
    } else {
        const auto &reactions = _context.get().reactions();
        for (const auto &entry : reactions.order1()) {
            for (auto reaction : entry.second) {
                reactionCounts()[reaction->id()] = 0;
            }
        }
        for (const auto &entry : reactions.order2()) {
            for (auto reaction : entry.second) {
                reactionCounts()[reaction->id()] = 0;
            }
        }
    }
}


}
}
}