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
 * @file Burst.h
 * @brief << brief description >>
 * @author luigisbailo
 * @date 21.02.18
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/mdgfrd/MDGFRDKernel.h>
#include <readdy/common/thread/barrier.h>


namespace readdy {
namespace kernel {
namespace mdgfrd {
namespace actions {
class MDGFRDBurst : public readdy::model::actions::Burst {
    using data_bounds = std::tuple<data::EntryDataContainer::iterator, data::EntryDataContainer::iterator>;
    using nl_bounds = std::tuple<std::size_t, std::size_t>;

public:

    explicit MDGFRDBurst(MDGFRDKernel *kernel) : kernel(kernel) {}

    void perform(const util::PerformanceNode &node) override;

private:

    MDGFRDKernel *const kernel;
};
}
}
}
}
