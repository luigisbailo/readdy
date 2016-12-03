/**
 * Test different contexts w.r.t. boxsize and periodicity, perform setupBoxes() and see if that worked.
 * Then add a few particles and perform fillBoxes(). The result of this is the actual neighborlist, which is
 * checked as well.
 *
 * @file CPUTestNeighborList.cpp
 * @brief Test the neighborlist object of the CPU kernel.
 * @author chrisfroe
 * @author clonker
 * @date 23.08.16
 */

#include <gtest/gtest.h>
#include <readdy/kernel/cpu/model/NeighborList.h>
#include <readdy/kernel/cpu/Kernel.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/testing/NOOPPotential.h>

namespace cpu = readdy::kernel::cpu;
namespace scpu = readdy::kernel::singlecpu;
namespace cpum = cpu::model;
namespace scpum = scpu::model;
namespace m = readdy::model;

namespace {

using data_t = cpu::model::ParticleData;
using nl_t = cpu::model::NeighborList;

struct TestNeighborList : ::testing::Test {

    std::unique_ptr<cpu::Kernel> kernel;
    unsigned int typeIdA;

    TestNeighborList() : kernel(std::make_unique<cpu::Kernel>()) {
        auto &ctx = kernel->getKernelContext();
        ctx.setDiffusionConstant("A", 1.0);
        double eductDistance = 1.2;
        ctx.registerReaction(kernel->createFusionReaction("test", "A", "A", "A", 0., eductDistance));

        ctx.registerPotential(std::make_unique<readdy::testing::NOOPPotentialOrder2>(1.1, 0., 0.), "A", "A");
        typeIdA = ctx.getParticleTypeID("A");
        ctx.configure();
    }

};

auto isPairInList = [](nl_t *pairs, data_t &data,
                       unsigned long idx1, unsigned long idx2) {
    const auto &neighbors1 = pairs->find_neighbors(idx1);
    for (auto &neigh_idx : neighbors1) {
        if (neigh_idx == idx2) {
            return true;
        }
    }
    const auto &neighbors2 = pairs->find_neighbors(idx2);
    for (auto &neigh_idx : neighbors2) {
        if (neigh_idx == idx1) {
            return true;
        }
    }
    return false;
};

auto isIdPairInList = [](nl_t *pairs, data_t &data, std::size_t id1, std::size_t id2) {
    return isPairInList(pairs, data, data.getIndexForId(id1), data.getIndexForId(id2));
};

auto getNumberPairs = [](const readdy::kernel::cpu::model::NeighborList &pairs) {
    using val_t = decltype(*pairs.begin());
    return std::accumulate(pairs.begin(), pairs.end(), 0, [](int acc, val_t &x) {
        return acc + x.size();
    });
};

TEST_F(TestNeighborList, TestCellsDirty) {
    auto& ctx = kernel->getKernelContext();
    ctx.setBoxSize(10, 10, 10);
    ctx.configure();

    readdy::util::thread::Config conf;
    readdy::kernel::cpu::model::ParticleData data {&ctx};
    cpum::NeighborList list(&ctx, data, &conf);
    list.setSkinSize(1.0);

    // Add three particles, two are in one outer box, the third on the other end and thus no neighbor
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, 1.8, 0, typeIdA)
    };

    data.addParticles(particles);
    list.create();

    list.displace(data.entry_at(0), {1.2, .0, .0});
    auto dirtyCells = list.findDirtyCells();
    EXPECT_EQ(125, dirtyCells.size());

    list.create();
    dirtyCells = list.findDirtyCells();
    EXPECT_EQ(0, dirtyCells.size());

}

TEST_F(TestNeighborList, ThreeBoxesNonPeriodic) {
    // maxcutoff is 1.2, system is 1.5 x 4 x 1.5, non-periodic, three cells
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(1.5, 4, 1.5);
    ctx.setPeriodicBoundary(false, false, false);

    readdy::util::thread::Config conf;
    readdy::kernel::cpu::model::ParticleData data {&ctx};
    cpum::NeighborList list(&ctx, data, &conf);

    list.setupCells();
    // Add three particles, two are in one outer box, the third on the other end and thus no neighbor
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, 1.8, 0, typeIdA)
    };

    data.addParticles(particles);
    list.fillCells();
    int sum = getNumberPairs(list);
    EXPECT_EQ(sum, 2);
    EXPECT_TRUE(isPairInList(&list, data, 0, 1));
    EXPECT_TRUE(isPairInList(&list, data, 1, 0));
}

TEST_F(TestNeighborList, OneDirection) {
    // maxcutoff is 1.2, system is 4.8 x 5 x 5.1
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(1.2, 1.1, 2.8);
    ctx.setPeriodicBoundary(false, false, true);

    readdy::util::thread::Config conf;
    readdy::kernel::cpu::model::ParticleData data {&ctx};
    // Add three particles, one of which is in the neighborhood of the other two
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, -1.1, typeIdA), m::Particle(0, 0, .4, typeIdA), m::Particle(0, 0, 1.1, typeIdA)
    };
    std::vector<std::size_t> ids(particles.size());
    std::transform(particles.begin(), particles.end(), ids.begin(), [](const m::Particle& p) {return p.getId();});
    data.addParticles(particles);

    cpum::NeighborList list(&ctx, data, &conf);
    list.create();

    int sum = getNumberPairs(list);
    EXPECT_EQ(sum, 4);
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(0), ids.at(2)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(2), ids.at(0)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(1), ids.at(2)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(2), ids.at(1)));
    EXPECT_FALSE(isIdPairInList(&list, data, ids.at(0), ids.at(1)));
    EXPECT_FALSE(isIdPairInList(&list, data, ids.at(1), ids.at(0)));
}

TEST_F(TestNeighborList, AllNeighborsInCutoffSphere) {
    // maxcutoff is 1.2, system is 4 x 4 x 4, all directions periodic
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(4, 4, 4);
    ctx.setPeriodicBoundary(true, true, true);
    readdy::util::thread::Config conf;
    readdy::kernel::cpu::model::ParticleData data {&ctx};
    cpum::NeighborList list(&ctx, data, &conf);
    list.setupCells();
    // Create a few particles. In this box setup, all particles are neighbors.
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(.3, 0, 0, typeIdA),
            m::Particle(0, .3, -.3, typeIdA), m::Particle(-.3, 0, .3, typeIdA), m::Particle(.3, -.3, 0, typeIdA)
    };

    data.addParticles(particles);
    list.fillCells();
    int sum = getNumberPairs(list);
    EXPECT_EQ(sum, 30);
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = i + 1; j < 6; ++j) {
            EXPECT_TRUE(isPairInList(&list, data, i, j)) << "Particles " << i << " and " << j << " were not neighbors.";
            EXPECT_TRUE(isPairInList(&list, data, j, i)) << "Particles " << j << " and " << i << " were not neighbors.";
        }
    }
}
}