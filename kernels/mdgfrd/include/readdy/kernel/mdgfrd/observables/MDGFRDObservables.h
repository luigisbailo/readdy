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
 * @file Observables.h
 * @brief << brief description >>
 * @author clonker
 * @author luigisbailo
 * @date 14.1.18
 */

#pragma once
#include <readdy/model/observables/Observables.h>

namespace readdy {
namespace kernel {
namespace mdgfrd {
class MDGFRDKernel;

namespace observables {

class MDGFRDVirial : public readdy::model::observables::Virial {
public:
    MDGFRDVirial(MDGFRDKernel *kernel, stride_type stride);

    void evaluate() override;

protected:
    MDGFRDKernel *const kernel;
};

class MDGFRDPositions : public readdy::model::observables::Positions {
public:
    MDGFRDPositions(MDGFRDKernel* kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {});

    void evaluate() override;

protected:
    MDGFRDKernel *const kernel;
};

class MDGFRDParticles : public readdy::model::observables::Particles {
public:
    MDGFRDParticles(MDGFRDKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    MDGFRDKernel *const kernel;
};

class MDGFRDHistogramAlongAxis : public readdy::model::observables::HistogramAlongAxis {

public:
    MDGFRDHistogramAlongAxis(MDGFRDKernel* kernel, unsigned int stride,
                       const std::vector<scalar> &binBorders,
                       const std::vector<std::string> &typesToCount,
                       unsigned int axis);

    void evaluate() override;

protected:
    MDGFRDKernel *const kernel;
    size_t size;
};

class MDGFRDNParticles : public readdy::model::observables::NParticles {
public:

    MDGFRDNParticles(MDGFRDKernel* kernel, unsigned int stride, std::vector<std::string> typesToCount = {});


    void evaluate() override;

protected:
    MDGFRDKernel *const kernel;
};

class MDGFRDForces : public readdy::model::observables::Forces {
public:
    MDGFRDForces(MDGFRDKernel* kernel, unsigned int stride, std::vector<std::string> typesToCount = {});

    ~MDGFRDForces() override = default;

    void evaluate() override;


protected:
    MDGFRDKernel *const kernel;
};

class MDGFRDReactions : public readdy::model::observables::Reactions {
public:
    MDGFRDReactions(MDGFRDKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    MDGFRDKernel *const kernel;
};

class MDGFRDReactionCounts : public readdy::model::observables::ReactionCounts {
public:
    MDGFRDReactionCounts(MDGFRDKernel* kernel, unsigned int stride);

    void evaluate() override;

protected:
    MDGFRDKernel *const kernel;
};

}
}
}
}
