/**
 * << detailed description >>
 *
 * @file ObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.07.16
 */

#include <readdy/kernel/cpu/observables/ObservableFactory.h>
#include <readdy/kernel/singlecpu/observables/SingleCPUObservables.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace observables {
ObservableFactory::ObservableFactory(CPUKernel *const kernel) : readdy::model::_internal::ObservableFactory(kernel),
                                                                kernel(kernel) {
}

readdy::model::NParticlesObservable *
ObservableFactory::createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new readdy::kernel::singlecpu::observables::NParticlesObservable<CPUKernel>(kernel, stride, typesToCount);
}

readdy::model::HistogramAlongAxisObservable *
ObservableFactory::createAxisHistogramObservable(unsigned int stride, std::vector<double> binBorders,
                                                 std::vector<std::string> typesToCount,
                                                 unsigned int axis) const {
    return new readdy::kernel::singlecpu::observables::HistogramAlongAxisObservable<CPUKernel>(kernel, stride,
                                                                                               binBorders, typesToCount,
                                                                                               axis);
}

readdy::model::ForcesObservable *
ObservableFactory::createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new readdy::kernel::singlecpu::observables::ForcesObservable<CPUKernel>(kernel, stride, typesToCount);
}

readdy::model::ParticlePositionObservable *
ObservableFactory::createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new readdy::kernel::singlecpu::observables::ParticlePositionObservable<CPUKernel>(kernel, stride, typesToCount);
}

readdy::model::RadialDistributionObservable *
ObservableFactory::createRadialDistributionObservable(unsigned int stride, std::vector<double> binBorders, std::string typeCountFrom,
                                                      std::string typeCountTo, double particleToDensity) const {
    return new readdy::kernel::singlecpu::observables::RadialDistributionObservable<CPUKernel>(kernel, stride, binBorders, typeCountFrom, typeCountTo,
                                                                                               particleToDensity);
}

}
}
}
}