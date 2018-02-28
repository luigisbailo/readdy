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
 * @file MDGFRDKernel.h
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 14/02/18
 */

#pragma once

#include <readdy/model/Kernel.h>

#include "pool.h"
#include "MDGFRDStateModel.h"
#include "observables/MDGFRDObservableFactory.h"
#include "actions/MDGFRDActionFactory.h"
//#include "actions/topologies/MDGFRDTopologyActionFactory.h"

namespace readdy {
namespace kernel {
namespace mdgfrd {

class MDGFRDKernel : public readdy::model::Kernel {
public:
    static const std::string name;

    MDGFRDKernel();

    ~MDGFRDKernel() override = default;

    MDGFRDKernel(const MDGFRDKernel &) = delete;

    MDGFRDKernel &operator=(const MDGFRDKernel &) = delete;

    MDGFRDKernel(MDGFRDKernel &&) = delete;

    MDGFRDKernel &operator=(MDGFRDKernel &&) = delete;

    // factory method
    static readdy::model::Kernel *create();

    const MDGFRDStateModel &getMDGFRDKernelStateModel() const {
        return _stateModel;
    };

    MDGFRDStateModel &getMDGFRDKernelStateModel() {
        return _stateModel;
    };

    const model::StateModel &stateModel() const override {
        return _stateModel;
    };

    model::StateModel &stateModel() override {
        return _stateModel;
    };

    void setNThreads(std::uint32_t n) {
        _pool.resize_wait(n);
    };

    std::size_t getNThreads() {
        return static_cast<std::size_t>(_pool.size());
    }

    const model::actions::ActionFactory &actions() const override {
        return _actions;
    };

    model::actions::ActionFactory &actions() override {
        return _actions;
    };

//    const model::top::TopologyActionFactory *const getTopologyActionFactory() const override {
//        return &_topologyActionFactory;
//    };
//
//    model::top::TopologyActionFactory *const getTopologyActionFactory() override {
//        return &_topologyActionFactory;
//    };

    const model::observables::ObservableFactory &observe() const override {
        return _observables;
    };

    model::observables::ObservableFactory &observe() override {
        return _observables;
    };

    void initialize() override;


    void finalize() override {
        readdy::model::Kernel::finalize();
    }

    thread_pool &pool() {
        return _pool;
    }

    const thread_pool &pool() const {
        return _pool;
    }

protected:

    std::unordered_map<particle_type_type, scalar> minDomainRegistry;
    std::unordered_map<particle_type_type, scalar> maxDomainRegistry;
    MDGFRDStateModel::data_type _data;
    actions::MDGFRDActionFactory _actions;
    observables::MDGFRDObservableFactory _observables;
//    actions::top::MDGFRDTopologyActionFactory _topologyActionFactory;
    MDGFRDStateModel _stateModel;
    thread_pool _pool;
};

}
}
}

extern "C" const char *name();

extern "C" readdy::model::Kernel *createKernel();
