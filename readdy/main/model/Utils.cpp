/**
 * << detailed description >>
 *
 * @file Utils.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 17.11.16
 */


#include <algorithm>
#include <unordered_map>

#include <readdy/model/Utils.h>

namespace readdy {
namespace model {
namespace util {

double getMaximumDisplacement(KernelContext& context) {
    context.configure();

    double kbt = context.getKBT();
    double tau = context.getTimeStep();

    double maximum_displacement = 0;
    for (auto &&pI : context.getAllRegisteredParticleTypes()) {
        double D = context.getDiffusionConstant(pI);
        double fMax = 0;

        for (auto &&pJ : context.getAllRegisteredParticleTypes()) {

            for (auto &&pot : context.getOrder2Potentials(pI, pJ)) {
                if (pot->getCutoffRadius() > 0) {
                    fMax = std::max(pot->getMaximalForce(kbt), fMax);
                }
            }
        }
        // todo magic number?
        double local_maximum_displacement = std::sqrt(2 * D * tau) + D * kbt * fMax * tau;
        maximum_displacement = std::max(local_maximum_displacement, maximum_displacement);
    }

    return maximum_displacement;
}

double getRecommendedTimeStep(unsigned int N, KernelContext& context) {
    double tau_R = 0;

    context.configure();

    double kbt = context.getKBT();
    double kReactionMax = 0;

    for (auto &&reactionO1 : context.getAllOrder1Reactions()) {
        kReactionMax = std::max(kReactionMax, reactionO1->getRate());
    }
    for (auto &&reactionO2 : context.getAllOrder2Reactions()) {
        kReactionMax = std::max(kReactionMax, reactionO2->getRate());
    }

    double tDMin = 0;
    std::unordered_map<unsigned int, double> fMaxes;
    for (auto &&pI : context.getAllRegisteredParticleTypes()) {
        double D = context.getDiffusionConstant(pI);
        double tD = 0;
        double xi = 0; // 1/(beta*Fmax)
        double fMax = 0;
        double rMin = std::numeric_limits<double>::max();

        for (auto &&reaction : context.getOrder1Reactions(pI)) {
            if (reaction->getNProducts() == 2 && reaction->getProductDistance() > 0) {
                rMin = std::min(rMin, reaction->getProductDistance());
            }
        }

        for (auto &&pot : context.getOrder1Potentials(pI)) {
            fMax = std::max(pot->getMaximalForce(kbt), fMax);
            if (pot->getRelevantLengthScale() > 0) {
                rMin = std::min(rMin, pot->getRelevantLengthScale());
            }
        }

        for (auto &&pJ : context.getAllRegisteredParticleTypes()) {

            for (auto &&reaction : context.getOrder2Reactions(pI, pJ)) {
                if (reaction->getEductDistance() > 0) {
                    rMin = std::min(rMin, reaction->getEductDistance());
                }
                if (reaction->getNProducts() == 2 && reaction->getProductDistance() > 0) {
                    rMin = std::min(rMin, reaction->getProductDistance());
                }
            }

            for (auto &&pot : context.getOrder2Potentials(pI, pJ)) {
                if (pot->getCutoffRadius() > 0) {
                    rMin = std::min(rMin, pot->getCutoffRadius());
                    fMax = std::max(pot->getMaximalForce(kbt), fMax);
                } else {
                    log::console()->warn("Discovered potential with cutoff radius 0.");
                }
            }
        }
        double rho = rMin / 2;
        if (fMax > 0) {
            xi = 1. / (context.getKBT() * fMax);
            tD = (xi * xi / D) * (1 + rho / xi - std::sqrt(1 + 2 * rho / xi));
        } else if (D > 0) {
            tD = .5 * rho * rho / D;
        }
        fMaxes.emplace(pI, fMax);
        log::console()->trace(" tau for {}: {} ( xi = {}, rho = {})", context.getParticleName(pI), tD, xi, rho);
        if (tDMin == 0) {
            tDMin = tD;
        } else {
            tDMin = std::min(tDMin, tD);
        }
    }

    log::console()->debug("Maximal displacement for particle types per time step (stochastic + deterministic): ");
    for (auto &&pI : context.getAllRegisteredParticleTypes()) {
        double D = context.getDiffusionConstant(pI);
        double xmax = std::sqrt(2 * D * tDMin) + D * kbt * fMaxes[pI] * tDMin;
        log::console()->debug("\t - {}: {} + {} = {}" , context.getParticleName(pI), std::sqrt(2 * D * tDMin),
                              D * kbt * fMaxes[pI] * tDMin, xmax);
    }

    if (kReactionMax > 0) tau_R = 1. / kReactionMax;

    double tau = std::max(tau_R, tDMin);
    if (tau_R > 0) tau = std::min(tau_R, tau);
    if (tDMin > 0) tau = std::min(tDMin, tau);
    tau /= (double) N;
    log::console()->debug("Estimated time step: {}", tau);
    return tau;
}

}
}
}