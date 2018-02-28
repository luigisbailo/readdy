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
 * @file UpdateNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 14.02.18
 */


#pragma once
#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/mdgfrd/MDGFRDKernel.h>

namespace readdy {
namespace kernel {
namespace mdgfrd {
namespace actions {
class MDGFRDUpdateNeighborList : public readdy::model::actions::UpdateNeighborList {
    using super = readdy::model::actions::UpdateNeighborList;
public:

    MDGFRDUpdateNeighborList(MDGFRDKernel *kernel, super::Operation op, scalar skin) : super(op, skin), kernel(kernel) {}

    void perform(const util::PerformanceNode &node) override {
        auto t = node.timeit();
        switch (operation) {
            case init:
                kernel->getMDGFRDKernelStateModel().initializeNeighborList(skinSize, node);
                break;
            case clear:
                kernel->getMDGFRDKernelStateModel().clearNeighborList(node);
                break;
            case update:
                kernel->getMDGFRDKernelStateModel().updateNeighborList(node);
                break;
        }

    }

    bool supportsSkin() const override {
        return true;
    }

private:
    MDGFRDKernel *kernel;
};
}
}
}
}
