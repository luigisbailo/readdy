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
 * @file MDGFRDProgramFactory.h
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 14.2.18
 */

#pragma once

#include <readdy/model/actions/ActionFactory.h>

namespace readdy {
namespace kernel {
namespace mdgfrd {
class MDGFRDKernel;
namespace actions {
class MDGFRDActionFactory : public readdy::model::actions::ActionFactory {
    MDGFRDKernel *const kernel;
public:
    explicit MDGFRDActionFactory(MDGFRDKernel *kernel);

    std::unique_ptr<model::actions::AddParticles> addParticles(const std::vector<model::Particle> &particles) const override;

    std::unique_ptr<model::actions::EulerBDIntegrator> eulerBDIntegrator(scalar timeStep) const override;

    std::unique_ptr<model::actions::CalculateForces> calculateForces() const override;

    std::unique_ptr<model::actions::Burst> burst() const override;

    std::unique_ptr<model::actions::UpdateNeighborList>
    updateNeighborList(model::actions::UpdateNeighborList::Operation operation, scalar skinSize) const override;

    std::unique_ptr<model::actions::reactions::UncontrolledApproximation>
    uncontrolledApproximation(scalar timeStep) const override;

    std::unique_ptr<model::actions::reactions::Gillespie> gillespie(scalar timeStep) const override;

};

}
}
}
}
