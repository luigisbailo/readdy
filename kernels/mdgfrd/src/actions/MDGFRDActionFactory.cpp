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
 * @file MDGFRDActionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 15.02.18
 */

#include <readdy/kernel/mdgfrd/actions/MDGFRDActionFactory.h>
#include <readdy/kernel/mdgfrd/actions/MDGFRDEulerBDIntegrator.h>
#include <readdy/kernel/mdgfrd/actions/MDGFRDUpdateNeighborList.h>
#include <readdy/kernel/mdgfrd/actions/MDGFRDCalculateForces.h>
#include <readdy/kernel/mdgfrd/actions/MDGFRDBurst.h>
#include <readdy/kernel/mdgfrd/actions/reactions/MDGFRDGillespie.h>
#include <readdy/kernel/mdgfrd/actions/reactions/MDGFRDUncontrolledApproximation.h>

namespace core_p = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace mdgfrd {
namespace actions {

MDGFRDActionFactory::MDGFRDActionFactory(MDGFRDKernel *const kernel) : kernel(kernel) { }

std::unique_ptr<model::actions::AddParticles>
MDGFRDActionFactory::addParticles(const std::vector<model::Particle> &particles) const {
    return {std::make_unique<readdy::model::actions::AddParticles>(kernel, particles)};
}

std::unique_ptr<model::actions::EulerBDIntegrator> MDGFRDActionFactory::eulerBDIntegrator(scalar timeStep) const {
    return {std::make_unique<MDGFRDEulerBDIntegrator>(kernel, timeStep)};
}

std::unique_ptr<model::actions::CalculateForces> MDGFRDActionFactory::calculateForces() const {
    return {std::make_unique<MDGFRDCalculateForces>(kernel)};
}

std::unique_ptr<model::actions::Burst> MDGFRDActionFactory::burst() const {
    return {std::make_unique<MDGFRDBurst>(kernel)};
}

std::unique_ptr<model::actions::UpdateNeighborList>
MDGFRDActionFactory::updateNeighborList(model::actions::UpdateNeighborList::Operation operation, scalar skinSize) const {
    return {std::make_unique<MDGFRDUpdateNeighborList>(kernel, operation, skinSize)};
}

std::unique_ptr<model::actions::reactions::UncontrolledApproximation>
MDGFRDActionFactory::uncontrolledApproximation(scalar timeStep) const {
    return {std::make_unique<reactions::MDGFRDUncontrolledApproximation>(kernel, timeStep)};
}

std::unique_ptr<model::actions::reactions::Gillespie> MDGFRDActionFactory::gillespie(scalar timeStep) const {
    return {std::make_unique<reactions::MDGFRDGillespie>(kernel, timeStep)};
}


}
}
}
}
