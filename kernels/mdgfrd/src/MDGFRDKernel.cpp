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
 * @file MDGFRDKernel.cpp
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 15/02/18
 */

#include <readdy/kernel/mdgfrd/MDGFRDKernel.h>


namespace readdy {
namespace kernel {
namespace mdgfrd {

const std::string MDGFRDKernel::name = "MDGFRD";

readdy::model::Kernel *MDGFRDKernel::create() {
    return new MDGFRDKernel();
}

MDGFRDKernel::MDGFRDKernel() : readdy::model::Kernel(name), _pool(readdy_default_n_threads()),
                         _data(_context, _pool), _actions(this),
                         _observables(this),
                         _stateModel(_data, _context, _pool) {}

void MDGFRDKernel::initialize() {
    readdy::model::Kernel::initialize();

    const auto &fullConfiguration = context().kernelConfiguration();

    const auto &configuration = fullConfiguration.mdgfrd;
    // thread config
    setNThreads(static_cast<std::uint32_t>(configuration.threadConfig.getNThreads()));
    {
        // state model config
        _stateModel.configure(configuration);
    }
    _stateModel.reactionRecords().clear();
    _stateModel.resetReactionCounts();
    _stateModel.virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};


    auto pot2 = _context.potentials().potentialsOrder2();
    auto typ = _context.particle_types().typeMapping();
    auto x_boxSize = _context.boxSize()[0];

    unsigned long idx = 0;

    for (auto typit1=typ.begin(); typit1!=typ.end(); ++typit1  ){
        scalar maxDist = x_boxSize/2;
        auto D = _context.particle_types().diffusionConstantOf(typit1->second);
        minDomainRegistry.emplace(typit1->second,12*sqrt(D));
        for (auto typit2=typ.begin(); typit2!=typ.end(); ++typit2){

            auto pot2 = _context.potentials().potentialsOrder2().find(std::tie(typit1->second,typit2->second));
            scalar tempDist = (x_boxSize - pot2->getCuffRadius())/2;
            if ( tempDist<maxDist ){
                maxDist = tempDist;
            }

        }
        maxDomainRegistry.emplace(typit1->second, maxDist);
    }



}
}
}
}

const char *name() {
    return readdy::kernel::mdgfrd::MDGFRDKernel::name.c_str();
}

readdy::model::Kernel *createKernel() {
    return readdy::kernel::mdgfrd::MDGFRDKernel::create();
}
