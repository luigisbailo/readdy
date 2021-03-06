/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file CompartmentRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.09.17
 * @copyright BSD-3
 */

#pragma once

#include <readdy/common/macros.h>
#include <readdy/model/ParticleTypeRegistry.h>
#include <readdy/model/_internal/Util.h>
#include "Compartment.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(compartments)

class CompartmentRegistry {
public:

    using CompartmentVector = std::vector<std::shared_ptr<readdy::model::compartments::Compartment>>;

    explicit CompartmentRegistry(const ParticleTypeRegistry &types);

    Compartment::id_type addSphere(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                                   const Vec3 &origin, scalar radius, bool largerOrLess);

    Compartment::id_type addSphere(const Compartment::label_conversion_map &conversions, const std::string &uniqueName,
                                   const Vec3 &origin, scalar radius, bool largerOrLess) {
        return addSphere(_internal::util::transformTypesMap(conversions, _types.get()), uniqueName, origin, radius,
                         largerOrLess);
    }

    Compartment::id_type addPlane(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                                  const Vec3 &normalCoefficients, scalar distance, bool largerOrLess);

    Compartment::id_type addPlane(const Compartment::label_conversion_map &conversions, const std::string &uniqueName,
                                  const Vec3 &normalCoefficients, scalar distance, bool largerOrLess) {
        return addPlane(_internal::util::transformTypesMap(conversions, _types.get()), uniqueName, normalCoefficients,
                        distance, largerOrLess);
    }


    const CompartmentVector &get() const {
        return _compartments;
    }

    CompartmentVector &get() {
        return _compartments;
    }

private:
     CompartmentVector _compartments;

    std::reference_wrapper<const ParticleTypeRegistry> _types;
};

NAMESPACE_END(compartments)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
