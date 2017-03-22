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


#include <array>

#include <gtest/gtest.h>

#include <readdy/common/logging.h>
#include <readdy/common/ParticleTypeTuple.h>
#include <readdy/model/topologies/connectivity/Graph.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/model/topologies/TorsionPotential.h>
#include <readdy/common/numeric.h>
#include <readdy/testing/Utils.h>

struct TestTopologyGraphs : KernelTest {

};

/**
 * << detailed description >>
 *
 * @file TestTopologyGraphs.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

using particle_t = readdy::model::Particle;
using topology_particle_t = readdy::model::TopologyParticle;

using dihedral_bond = readdy::model::top::CosineDihedralPotential;

TEST(TestTopologyGraphs, TestQuadruple) {
    readdy::util::particle_type_quadruple_hasher hasher;
    std::array<readdy::particle_type_type, 4> range{1, 2, 3, 4};
    do {
        std::stringstream ss;
        ss << range[0] << ", " << range[1] << ", " << range[2] << ", " << range[3];
        EXPECT_EQ(readdy::util::sortTypeQuadruple(range[0], range[1], range[2], range[3]), std::make_tuple(1, 2, 3, 4))
                            << "failed for range " << ss.str() << ", should always yield a sorted tuple!";
    } while (std::next_permutation(range.begin(), range.end()));
    EXPECT_EQ(hasher(std::make_tuple(1, 2, 3, 4)), hasher(std::make_tuple(4, 3, 2, 1)));
    EXPECT_EQ(hasher(std::make_tuple(1, 3, 2, 4)), hasher(std::make_tuple(4, 2, 3, 1)));
}

TEST(TestTopologyGraphs, TestTriple) {
    readdy::util::particle_type_triple_hasher hasher;
    std::array<readdy::particle_type_type, 3> range{1, 2, 3};
    do {
        std::stringstream ss;
        ss << range[0] << ", " << range[1] << ", " << range[2];
        EXPECT_EQ(readdy::util::sortTypeTriple(range[0], range[1], range[2]), std::make_tuple(1, 2, 3))
                            << "failed for range " << ss.str() << ", should always yield a sorted tuple!";
    } while (std::next_permutation(range.begin(), range.end()));
    EXPECT_EQ(hasher(std::make_tuple(1, 2, 3)), hasher(std::make_tuple(3, 2, 1)));
    EXPECT_EQ(hasher(std::make_tuple(2, 1, 3)), hasher(std::make_tuple(3, 1, 2)));
}

TEST(TestTopologyGraphs, TestTuple) {
    readdy::util::particle_type_pair_hasher hasher;
    EXPECT_EQ(hasher(std::make_tuple(1, 2)), hasher(std::make_tuple(2, 1)));
    std::unordered_map<readdy::util::particle_type_pair, int, readdy::util::particle_type_pair_hasher, readdy::util::particle_type_pair_equal_to> map;
    int a = 1, b = 2;
    map[std::make_tuple(a, b)] = 5;
    EXPECT_EQ(map[std::make_tuple(1, 2)], 5);
    EXPECT_EQ(map[std::make_tuple(2, 1)], 5);
}

TEST(TestTopologyGraphs, TestGraphWithNamedVertices) {
    readdy::model::top::graph::Graph graph;
    graph.addVertex(0, 0, "myVertex");
    graph.addVertex(1, 0, "myVertex2");
    graph.addEdge("myVertex", "myVertex2");
    EXPECT_EQ(graph.vertices().size(), 2);
    EXPECT_EQ(graph.vertices().front().label, "myVertex");
    EXPECT_EQ(graph.vertices().back().label, "myVertex2");
    EXPECT_EQ(graph.namedVertex("myVertex").neighbors().size(), 1);
    EXPECT_EQ(graph.namedVertex("myVertex2").neighbors().size(), 1);
    EXPECT_EQ(graph.namedVertex("myVertex").neighbors().front()->label, "myVertex2");
    EXPECT_EQ(graph.namedVertex("myVertex2").neighbors().back()->label, "myVertex");
    graph.removeVertex("myVertex");
    EXPECT_EQ(graph.vertices().size(), 1);
    EXPECT_EQ(graph.vertices().front().label, "myVertex2");
    EXPECT_EQ(graph.namedVertex("myVertex2").neighbors().size(), 0);
}

TEST(TestTopologyGraphs, TestGraphWithIndices) {
    readdy::model::top::graph::Graph graph;
    graph.addVertex(0, 0);
    graph.addVertex(1, 0);
    graph.addEdge(graph.vertices().begin(), ++graph.vertices().begin());
    EXPECT_EQ(graph.vertices().size(), 2);
    EXPECT_EQ(graph.vertexForParticleIndex(0).particleIndex, 0);
    EXPECT_EQ(graph.vertexForParticleIndex(1).particleIndex, 1);
    EXPECT_EQ(graph.vertexForParticleIndex(0).neighbors().size(), 1);
    EXPECT_EQ(graph.vertexForParticleIndex(1).neighbors().size(), 1);
    EXPECT_EQ(graph.vertexForParticleIndex(0).neighbors()[0]->particleIndex, 1);
    EXPECT_EQ(graph.vertexForParticleIndex(1).neighbors()[0]->particleIndex, 0);
    graph.removeParticle(0);
    EXPECT_EQ(graph.vertices().size(), 1);
    EXPECT_EQ(graph.vertexForParticleIndex(1).particleIndex, 1);
    EXPECT_EQ(graph.vertexForParticleIndex(1).neighbors().size(), 0);
}

TEST(TestTopologyGraphs, TestTopologyWithGraph) {
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    auto &ctx = kernel->getKernelContext();

    ctx.registerParticleType("Topology A", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
    ctx.registerParticleType("Topology B", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);

    ctx.configureTopologyBondPotential("Topology A", "Topology B", {1.0, 1.0});
    ctx.configureTopologyAnglePotential("Topology B", "Topology A", "Topology A", {1.0, 1.0});
    ctx.configureTopologyTorsionPotential("Topology A", "Topology B", "Topology A", "Topology A", {1.0, 1.0, 3.0});

    ctx.setBoxSize(10, 10, 10);
    topology_particle_t x_i{-1, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_j{0, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_k{0, 0, 1, ctx.getParticleTypeID("Topology B")};
    topology_particle_t x_l{1, .1, 1, ctx.getParticleTypeID("Topology A")};

    auto top = kernel->getKernelStateModel().addTopology({x_i, x_j, x_k, x_l});
    EXPECT_EQ(top->graph().vertices().size(), 4);
    auto it = top->graph().vertices().begin();
    auto it2 = ++top->graph().vertices().begin();
    EXPECT_FALSE(top->graph().isConnected());
    top->graph().addEdge(it++, it2++);
    EXPECT_FALSE(top->graph().isConnected());
    top->graph().addEdge(it++, it2++);
    EXPECT_FALSE(top->graph().isConnected());
    top->graph().addEdge(it++, it2++);
    EXPECT_TRUE(top->graph().isConnected());

    top->graph().setVertexLabel(top->graph().firstVertex(), "begin");
    top->graph().setVertexLabel(top->graph().lastVertex(), "end");

    top->graph().addEdge("begin", "end");

    top->configure();
}


TEST_P(TestTopologyGraphs, DihedralPotentialSteeperAngle) {
    auto &ctx = kernel->getKernelContext();
    ctx.registerParticleType("Topology A", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
    ctx.setBoxSize(10, 10, 10);
    topology_particle_t x_i{-1, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_j{0, 0, 0, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_k{0, 0, 1, ctx.getParticleTypeID("Topology A")};
    topology_particle_t x_l{1, 3, 1, ctx.getParticleTypeID("Topology A")};
    auto top = kernel->getKernelStateModel().addTopology({x_i, x_j, x_k, x_l});
    auto it = top->graph().vertices().begin();
    auto it2 = ++top->graph().vertices().begin();
    while(it2 != top->graph().vertices().end()) {
        top->graph().addEdge(it, it2);
        std::advance(it, 1);
        std::advance(it2, 1);
    }
    kernel->getKernelContext().configureTopologyTorsionPotential("Topology A", "Topology A", "Topology A", "Topology A", {1.0, 3, readdy::util::numeric::pi()});
    top->configure();
    auto fObs = kernel->createObservable<readdy::model::observables::Forces>(1);
    std::vector<readdy::model::Vec3> collectedForces;
    fObs->setCallback([&collectedForces](const readdy::model::observables::Forces::result_t &result) {
        for (const auto &force : result) {
            collectedForces.push_back(force);
        }
    });

    auto conn = kernel->connectObservable(fObs.get());

    ctx.configure();
    kernel->getKernelStateModel().calculateForces();
    kernel->evaluateObservables(1);

    EXPECT_EQ(collectedForces.size(), 4);
    EXPECT_DOUBLE_EQ(kernel->getKernelStateModel().getEnergy(), 1.8221921916437787);
    readdy::model::Vec3 force_x_i{0., 1.70762994, 0.};
    readdy::model::Vec3 force_x_j{0., -1.70762994, 0.};
    readdy::model::Vec3 force_x_k{0.51228898, -0.17076299, 0.};
    readdy::model::Vec3 force_x_l{-0.51228898, 0.17076299, 0.};
    EXPECT_VEC3_NEAR(collectedForces[0], force_x_i, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[1], force_x_j, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[2], force_x_k, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[3], force_x_l, 1e-6);
}

INSTANTIATE_TEST_CASE_P(TestTopologyGraphsCore, TestTopologyGraphs, ::testing::Values("SingleCPU", "CPU"));

