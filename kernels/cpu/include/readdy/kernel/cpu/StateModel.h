/**
 * << detailed description >>
 *
 * @file CPUStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_CPUKERNEL_CPUSTATEMODEL_H
#define READDY_CPUKERNEL_CPUSTATEMODEL_H


#include <readdy/model/KernelStateModel.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/kernel/cpu/model/NeighborList.h>
#include <readdy/common/thread/Config.h>
#include <readdy/kernel/cpu/model/ParticleData.h>

namespace readdy {
namespace kernel {
namespace cpu {
class StateModel : public readdy::model::KernelStateModel {

public:

    using data_t = readdy::kernel::cpu::model::ParticleData;

    StateModel(readdy::model::KernelContext *const context, readdy::util::thread::Config const *const config);

    ~StateModel();

    virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;

    virtual const std::vector<readdy::model::Particle> getParticles() const override;

    virtual void updateNeighborList() override;

    virtual void calculateForces() override;

    virtual void addParticle(const readdy::model::Particle &p) override;

    virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;

    virtual void removeParticle(const readdy::model::Particle &p) override;

    virtual void removeAllParticles() override;

    virtual double getEnergy() const override;

    data_t *const getParticleData() const;

    model::NeighborList *const getNeighborList() const;

    virtual void clearNeighborList() override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
    readdy::util::thread::Config const *const config;
};
}
}
}

#endif //READDY_CPUKERNEL_CPUSTATEMODEL_H