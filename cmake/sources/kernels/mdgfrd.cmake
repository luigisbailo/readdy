#####################################################################
# Copyright (c) 2016 Computational Molecular Biology Group,         #
#                    Freie Universitaet Berlin (GER)                #
#                                                                   #
# This file is part of ReaDDy.                                      #
#                                                                   #
# ReaDDy is free software: you can redistribute it and/or modify    #
# it under the terms of the GNU Lesser General Public License as    #
# published by the Free Software Foundation, either version 3 of    #
# the License, or (at your option) any later version.               #
#                                                                   #
# This program is distributed in the hope that it will be useful,   #
# but WITHOUT ANY WARRANTY; without even the implied warranty of    #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     #
# GNU Lesser General Public License for more details.               #
#                                                                   #
# You should have received a copy of the GNU Lesser General         #
# Public License along with this program. If not, see               #
# <http://www.gnu.org/licenses/>.                                   #
#####################################################################


SET(SOURCES_DIR "${READDY_GLOBAL_DIR}/kernels/mdgfrd/src")
SET(MDGFRD_INCLUDE_DIR "${READDY_GLOBAL_DIR}/kernels/mdgfrd/include")

LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/test.cpp")


# --- main sources ---
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/MDGFRDKernel.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/MDGFRDStateModel.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/observables/MDGFRDObservableFactory.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/observables/MDGFRDObservables.cpp")

# --- neighbor list ---
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/nl/CellLinkedList.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/nl/ContiguousCellLinkedList.cpp")

# --- actions ---
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/actions/MDGFRDActionFactory.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/actions/MDGFRDEulerBDIntegrator.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/actions/MDGFRDCalculateForces.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/actions/MDGFRDBurst.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/actions/reactions/ReactionUtils.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/actions/reactions/Event.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/actions/reactions/MDGFRDUncontrolledApproximation.cpp")
LIST(APPEND MDGFRD_SOURCES "${SOURCES_DIR}/actions/reactions/MDGFRDGillespie.cpp")
