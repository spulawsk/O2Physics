// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file heavyStrangeBaryonResonanceSearch.cxx
/// \brief Data-only search for excited Xi and Omega baryons in pp collisions.
///
/// Channels:
///   Xi(1530) and higher Xi* -> Xi pi
///   Xi*                     -> Lambda K
///   Omega*                  -> Xi K0S
///
/// The task consumes the standard PWGLF resonance-derived tables produced by
/// resonanceInitializer.  It deliberately makes no assumption about a signal
/// PDG code and stores continuous invariant-mass spectra.

#include "PWGLF/DataModel/LFResonanceTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>

#include <array>
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

namespace
{
using FourVector = ROOT::Math::PxPyPzEVector;

template <typename T>
FourVector fourVector(const T &candidate, float mass)
{
  const auto energy = std::sqrt(candidate.px() * candidate.px() +
                                candidate.py() * candidate.py() +
                                candidate.pz() * candidate.pz() + mass * mass);
  return FourVector{candidate.px(), candidate.py(), candidate.pz(), energy};
}

template <typename TrackT, typename V0T>
bool sharesTrack(const TrackT &track, const V0T &v0)
{
  const auto *indices = v0.indices();
  return track.trackId() == indices[0] || track.trackId() == indices[1];
}

template <typename TrackT, typename CascT>
bool sharesTrackWithCascade(const TrackT &track, const CascT &cascade)
{
  const auto *indices = cascade.cascadeIndices();
  return track.trackId() == indices[0] || track.trackId() == indices[1] ||
         track.trackId() == indices[2];
}

template <typename CascT, typename V0T>
bool sharesTrackCascadeV0(const CascT &cascade, const V0T &v0)
{
  const auto *cascadeIndices = cascade.cascadeIndices();
  const auto *v0Indices = v0.indices();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      if (cascadeIndices[i] == v0Indices[j]) {
        return true;
      }
    }
  }
  return false;
}
}  // namespace

struct HeavyStrangeBaryonResonanceSearch {
  SliceCache cache;
  HistogramRegistry registry{
    "registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> requireINELgt0{"requireINELgt0", false,
                                    "Require reconstructed INEL>0"};
  Configurable<float> maxAbsRapidity{"maxAbsRapidity", 0.5f,
                                     "Maximum absolute candidate rapidity"};

  Configurable<float> minTrackPt{"minTrackPt", 0.15f,
                                 "Minimum primary pion/kaon pT"};
  Configurable<float> maxTrackEta{"maxTrackEta", 0.8f,
                                  "Maximum primary pion/kaon |eta|"};
  Configurable<float> maxTrackDCAxy{"maxTrackDCAxy", 0.1f,
                                    "Maximum primary pion/kaon |DCAxy|"};
  Configurable<float> maxTrackDCAz{"maxTrackDCAz", 0.2f,
                                   "Maximum primary pion/kaon |DCAz|"};
  Configurable<int> minTPCClusters{"minTPCClusters", 70,
                                   "Minimum found TPC clusters"};
  Configurable<float> maxTPCNSigmaPion{"maxTPCNSigmaPion", 3.f,
                                       "Maximum pion TPC n-sigma"};
  Configurable<float> maxTOFNSigmaPion{"maxTOFNSigmaPion", 3.f,
                                       "Maximum pion TOF n-sigma"};
  Configurable<float> maxTPCNSigmaKaon{"maxTPCNSigmaKaon", 3.f,
                                       "Maximum kaon TPC n-sigma"};
  Configurable<float> maxTOFNSigmaKaon{"maxTOFNSigmaKaon", 3.f,
                                       "Maximum kaon TOF n-sigma"};
  Configurable<bool> requireTOFAbovePt{"requireTOFAbovePt", true,
                                       "Require TOF above PID threshold"};
  Configurable<float> tofPtThreshold{"tofPtThreshold", 0.8f,
                                     "pT above which TOF is required"};

  Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.006f,
                                       "Lambda mass window"};
  Configurable<float> lambdaMinCosPA{"lambdaMinCosPA", 0.995f,
                                     "Minimum Lambda cos(PA)"};
  Configurable<float> lambdaMaxDaughterDCA{"lambdaMaxDaughterDCA", 0.5f,
                                           "Maximum Lambda daughter DCA"};
  Configurable<float> lambdaMinRadius{"lambdaMinRadius", 0.5f,
                                      "Minimum Lambda decay radius"};
  Configurable<float> lambdaMaxRadius{"lambdaMaxRadius", 200.f,
                                      "Maximum Lambda decay radius"};
  Configurable<float> lambdaDaughterPID{"lambdaDaughterPID", 5.f,
                                        "Maximum Lambda daughter TPC n-sigma"};

  Configurable<float> xiMassWindow{"xiMassWindow", 0.008f, "Xi mass window"};
  Configurable<float> xiMinV0CosPA{"xiMinV0CosPA", 0.97f,
                                   "Minimum Xi V0 cos(PA)"};
  Configurable<float> xiMinCascCosPA{"xiMinCascCosPA", 0.97f,
                                     "Minimum Xi cascade cos(PA)"};
  Configurable<float> xiMaxV0DaughterDCA{"xiMaxV0DaughterDCA", 0.8f,
                                         "Maximum V0 daughter DCA in Xi"};
  Configurable<float> xiMaxCascDaughterDCA{"xiMaxCascDaughterDCA", 0.8f,
                                           "Maximum cascade daughter DCA"};
  Configurable<float> xiMinRadius{"xiMinRadius", 0.8f,
                                  "Minimum Xi decay radius"};
  Configurable<float> xiBachelorPID{"xiBachelorPID", 5.f,
                                    "Maximum bachelor-pion TPC n-sigma"};

  Configurable<float> k0sMassWindow{"k0sMassWindow", 0.015f, "K0S mass window"};
  Configurable<float> k0sMinCosPA{"k0sMinCosPA", 0.98f, "Minimum K0S cos(PA)"};
  Configurable<float> k0sMaxDaughterDCA{"k0sMaxDaughterDCA", 0.5f,
                                        "Maximum K0S daughter DCA"};
  Configurable<float> k0sMinRadius{"k0sMinRadius", 0.5f,
                                   "Minimum K0S decay radius"};
  Configurable<float> k0sMaxRadius{"k0sMaxRadius", 200.f,
                                   "Maximum K0S decay radius"};
  Configurable<float> k0sDaughterPID{"k0sDaughterPID", 5.f,
                                     "Maximum K0S daughter-pion TPC n-sigma"};
  Configurable<bool> rejectLambdaInK0S{
    "rejectLambdaInK0S", true, "Reject Lambda-compatible K0S candidates"};
  Configurable<float> k0sLambdaRejectionWindow{
    "k0sLambdaRejectionWindow", 0.006f,
    "Lambda rejection window in K0S selection"};

  Configurable<int> numberMixedEvents{"numberMixedEvents", 5,
                                      "Number of mixed events"};
  ConfigurableAxis vertexMixingBins{"vertexMixingBins",
                                    {VARIABLE_WIDTH, -10.f, -8.f, -6.f, -4.f,
                                     -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f},
                                    "Vertex-z mixing bins"};
  ConfigurableAxis multiplicityMixingBins{
    "multiplicityMixingBins",
    {VARIABLE_WIDTH, 0.f, 1.f, 5.f, 10.f, 20.f, 30.f, 40.f, 50.f, 70.f, 100.f,
     110.f},
    "Multiplicity-percentile mixing bins"};

  using MixingBinning =
    ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  MixingBinning mixingBinning{{vertexMixingBins, multiplicityMixingBins}, true};

  void init(InitContext const &)
  {
    const AxisSpec eventCounter{5, 0.5, 5.5, "selection step"};
    const AxisSpec vertexAxis{120, -12., 12., "z_{PV} (cm)"};
    const AxisSpec multAxis{110, 0., 110., "multiplicity percentile"};
    const AxisSpec ptAxis{200, 0., 20., "p_{T} (GeV/c)"};
    const AxisSpec rapidityAxis{40, -1., 1., "y"};
    const AxisSpec pairTypeAxis{4, -0.5, 3.5, "pair type"};
    const AxisSpec xiPiMassAxis{1250, 1.45, 2.70, "m_{#Xi#pi} (GeV/c^{2})"};
    const AxisSpec lambdaKMassAxis{1091, 1.609, 2.70,
                                   "m_{#Lambda K} (GeV/c^{2})"};
    const AxisSpec xiK0sMassAxis{890, 1.81, 2.70,
                                 "m_{#Xi K^{0}_{S}} (GeV/c^{2})"};

    registry.add("Event/counter", "Event counter", HistType::kTH1F,
                 {eventCounter});
    registry.add("Event/vertexZ", "Accepted-event vertex", HistType::kTH1F,
                 {vertexAxis});
    registry.add("Event/multiplicity", "Accepted-event multiplicity",
                 HistType::kTH1F, {multAxis});

    registry.add("QA/lambdaMass", "Selected Lambda and anti-Lambda candidates",
                 HistType::kTH1F, {{200, 1.08, 1.16}});
    registry.add("QA/xiMass", "Selected Xi candidates", HistType::kTH1F,
                 {{240, 1.27, 1.39}});
    registry.add("QA/k0sMass", "Selected K0S candidates", HistType::kTH1F,
                 {{200, 0.44, 0.56}});
    registry.add("QA/pionTPC", "Selected primary pion TPC PID", HistType::kTH2F,
                 {ptAxis, {100, -5., 5.}});
    registry.add("QA/kaonTPC", "Selected primary kaon TPC PID", HistType::kTH2F,
                 {ptAxis, {100, -5., 5.}});

    registry.add("XiPi/same", "Xi-pion same event", HistType::kTHnSparseF,
                 {xiPiMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
    registry.add("XiPi/mixed", "Xi-pion mixed event", HistType::kTHnSparseF,
                 {xiPiMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
    registry.add(
      "LambdaK/same", "Lambda-kaon same event", HistType::kTHnSparseF,
      {lambdaKMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
    registry.add(
      "LambdaK/mixed", "Lambda-kaon mixed event", HistType::kTHnSparseF,
      {lambdaKMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
    registry.add("XiK0S/same", "Xi-K0S same event", HistType::kTHnSparseF,
                 {xiK0sMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
    registry.add("XiK0S/mixed", "Xi-K0S mixed event", HistType::kTHnSparseF,
                 {xiK0sMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
  }

  template <typename CollisionT>
  bool selectEvent(const CollisionT &collision, bool fillQA)
  {
    if (fillQA) {
      registry.fill(HIST("Event/counter"), 1.);
    }
    if (std::abs(collision.posZ()) >= 10.f) {
      return false;
    }
    if (requireINELgt0 && !collision.isRecINELgt0()) {
      return false;
    }
    if (fillQA) {
      registry.fill(HIST("Event/counter"), 2.);
      registry.fill(HIST("Event/vertexZ"), collision.posZ());
      registry.fill(HIST("Event/multiplicity"), collision.cent());
    }
    return true;
  }

  template <typename TrackT>
  bool selectPrimaryTrack(const TrackT &track) const
  {
    return track.pt() >= minTrackPt && std::abs(track.eta()) <= maxTrackEta &&
           std::abs(track.dcaXY()) <= maxTrackDCAxy &&
           std::abs(track.dcaZ()) <= maxTrackDCAz &&
           track.tpcNClsFound() >= minTPCClusters && track.isGlobalTrackWoDCA();
  }

  template <typename TrackT>
  bool selectPion(const TrackT &track) const
  {
    if (!selectPrimaryTrack(track) ||
        std::abs(track.tpcNSigmaPi()) > maxTPCNSigmaPion) {
      return false;
    }
    if (track.hasTOF()) {
      return std::abs(track.tofNSigmaPi()) <= maxTOFNSigmaPion;
    }
    return !(requireTOFAbovePt && track.pt() >= tofPtThreshold);
  }

  template <typename TrackT>
  bool selectKaon(const TrackT &track) const
  {
    if (!selectPrimaryTrack(track) ||
        std::abs(track.tpcNSigmaKa()) > maxTPCNSigmaKaon) {
      return false;
    }
    if (track.hasTOF()) {
      return std::abs(track.tofNSigmaKa()) <= maxTOFNSigmaKaon;
    }
    return !(requireTOFAbovePt && track.pt() >= tofPtThreshold);
  }

  template <typename V0T>
  bool selectLambda(const V0T &v0, bool antiLambda) const
  {
    const float mass = antiLambda ? v0.mAntiLambda() : v0.mLambda();
    const bool daughterPID =
      antiLambda
        ? (std::abs(v0.daughterTPCNSigmaPosPi()) <= lambdaDaughterPID &&
           std::abs(v0.daughterTPCNSigmaNegPr()) <= lambdaDaughterPID)
        : (std::abs(v0.daughterTPCNSigmaPosPr()) <= lambdaDaughterPID &&
           std::abs(v0.daughterTPCNSigmaNegPi()) <= lambdaDaughterPID);
    return std::abs(mass - MassLambda0) <= lambdaMassWindow &&
           v0.v0CosPA() >= lambdaMinCosPA &&
           v0.daughDCA() <= lambdaMaxDaughterDCA &&
           v0.transRadius() >= lambdaMinRadius &&
           v0.transRadius() <= lambdaMaxRadius && daughterPID;
  }

  template <typename CascT>
  bool selectXi(const CascT &xi) const
  {
    return std::abs(xi.mXi() - MassXiMinus) <= xiMassWindow &&
           xi.v0CosPA() >= xiMinV0CosPA && xi.cascCosPA() >= xiMinCascCosPA &&
           xi.daughDCA() <= xiMaxV0DaughterDCA &&
           xi.cascDaughDCA() <= xiMaxCascDaughterDCA &&
           xi.cascTransRadius() >= xiMinRadius &&
           std::abs(xi.daughterTPCNSigmaBachPi()) <= xiBachelorPID;
  }

  template <typename V0T>
  bool selectK0S(const V0T &k0s) const
  {
    if (std::abs(k0s.mK0Short() - MassK0Short) > k0sMassWindow ||
        k0s.v0CosPA() < k0sMinCosPA || k0s.daughDCA() > k0sMaxDaughterDCA ||
        k0s.transRadius() < k0sMinRadius || k0s.transRadius() > k0sMaxRadius ||
        std::abs(k0s.daughterTPCNSigmaPosPi()) > k0sDaughterPID ||
        std::abs(k0s.daughterTPCNSigmaNegPi()) > k0sDaughterPID) {
      return false;
    }
    return !rejectLambdaInK0S ||
           (std::abs(k0s.mLambda() - MassLambda0) > k0sLambdaRejectionWindow &&
            std::abs(k0s.mAntiLambda() - MassLambda0) >
              k0sLambdaRejectionWindow);
  }

  template <bool Mixed, typename CollisionT, typename CascadesT,
            typename TracksT>
  void fillXiPi(const CollisionT &collision, const CascadesT &cascades,
                const TracksT &tracks)
  {
    for (const auto &xi : cascades) {
      if (!selectXi(xi)) {
        continue;
      }
      if constexpr (!Mixed) {
        registry.fill(HIST("QA/xiMass"), xi.mXi());
      }
      const auto pXi = fourVector(xi, xi.mXi());
      for (const auto &pion : tracks) {
        if (!selectPion(pion) || (!Mixed && sharesTrackWithCascade(pion, xi))) {
          continue;
        }
        const auto candidate = pXi + fourVector(pion, MassPiPlus);
        if (std::abs(candidate.Rapidity()) > maxAbsRapidity) {
          continue;
        }
        // 0: Xi- pi+ signal, 1: anti-Xi+ pi- signal, 2/3: like sign.
        const int pairType = xi.sign() < 0 ? (pion.sign() > 0 ? 0 : 2)
                                           : (pion.sign() < 0 ? 1 : 3);
        registry.fill(Mixed ? HIST("XiPi/mixed") : HIST("XiPi/same"),
                      candidate.M(), candidate.Pt(), candidate.Rapidity(),
                      collision.cent(), pairType);
        if constexpr (!Mixed) {
          registry.fill(HIST("QA/pionTPC"), pion.pt(), pion.tpcNSigmaPi());
        }
      }
    }
  }

  template <bool Mixed, typename CollisionT, typename V0sT, typename TracksT>
  void fillLambdaK(const CollisionT &collision, const V0sT &v0s,
                   const TracksT &tracks)
  {
    for (const auto &lambda : v0s) {
      const bool isLambda = selectLambda(lambda, false);
      const bool isAntiLambda = selectLambda(lambda, true);
      if (!isLambda && !isAntiLambda) {
        continue;
      }
      for (const auto &kaon : tracks) {
        if (!selectKaon(kaon) || (!Mixed && sharesTrack(kaon, lambda))) {
          continue;
        }
        for (int anti = 0; anti < 2; ++anti) {
          if ((anti == 0 && !isLambda) || (anti == 1 && !isAntiLambda)) {
            continue;
          }
          const float lambdaMass =
            anti ? lambda.mAntiLambda() : lambda.mLambda();
          const auto candidate =
            fourVector(lambda, lambdaMass) + fourVector(kaon, MassKPlus);
          if (std::abs(candidate.Rapidity()) > maxAbsRapidity) {
            continue;
          }
          // 0: Lambda K-, 1: anti-Lambda K+, 2/3: wrong sign.
          const int pairType =
            anti == 0 ? (kaon.sign() < 0 ? 0 : 2) : (kaon.sign() > 0 ? 1 : 3);
          registry.fill(Mixed ? HIST("LambdaK/mixed") : HIST("LambdaK/same"),
                        candidate.M(), candidate.Pt(), candidate.Rapidity(),
                        collision.cent(), pairType);
          if constexpr (!Mixed) {
            registry.fill(HIST("QA/lambdaMass"), lambdaMass);
            registry.fill(HIST("QA/kaonTPC"), kaon.pt(), kaon.tpcNSigmaKa());
          }
        }
      }
    }
  }

  template <bool Mixed, typename CollisionT, typename CascadesT, typename V0sT>
  void fillXiK0S(const CollisionT &collision, const CascadesT &cascades,
                 const V0sT &v0s)
  {
    for (const auto &xi : cascades) {
      if (!selectXi(xi)) {
        continue;
      }
      for (const auto &k0s : v0s) {
        if (!selectK0S(k0s) || (!Mixed && sharesTrackCascadeV0(xi, k0s))) {
          continue;
        }
        const auto candidate =
          fourVector(xi, xi.mXi()) + fourVector(k0s, MassK0Short);
        if (std::abs(candidate.Rapidity()) > maxAbsRapidity) {
          continue;
        }
        const int pairType = xi.sign() < 0 ? 0 : 1;
        registry.fill(Mixed ? HIST("XiK0S/mixed") : HIST("XiK0S/same"),
                      candidate.M(), candidate.Pt(), candidate.Rapidity(),
                      collision.cent(), pairType);
        if constexpr (!Mixed) {
          registry.fill(HIST("QA/k0sMass"), k0s.mK0Short());
        }
      }
    }
  }

  void processData(const aod::ResoCollision &collision,
                   const aod::ResoTracks &tracks, const aod::ResoV0s &v0s,
                   const aod::ResoCascades &cascades)
  {
    if (!selectEvent(collision, true)) {
      return;
    }
    fillXiPi<false>(collision, cascades, tracks);
    fillLambdaK<false>(collision, v0s, tracks);
    fillXiK0S<false>(collision, cascades, v0s);
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processData,
                 "Process same-event data", true);

  void processMixedXiPi(const aod::ResoCollisions &collisions,
                        const aod::ResoTracks &tracks,
                        const aod::ResoCascades &cascades)
  {
    auto daughters = std::make_tuple(cascades, tracks);
    Pair<aod::ResoCollisions, aod::ResoCascades, aod::ResoTracks, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto &[collision1, cascades1, collision2, tracks2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillXiPi<true>(collision1, cascades1, tracks2);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedXiPi,
                 "Mixed-event Xi-pion", false);

  void processMixedLambdaK(const aod::ResoCollisions &collisions,
                           const aod::ResoTracks &tracks,
                           const aod::ResoV0s &v0s)
  {
    auto daughters = std::make_tuple(v0s, tracks);
    Pair<aod::ResoCollisions, aod::ResoV0s, aod::ResoTracks, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto &[collision1, v0s1, collision2, tracks2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillLambdaK<true>(collision1, v0s1, tracks2);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedLambdaK,
                 "Mixed-event Lambda-kaon", false);

  void processMixedXiK0S(const aod::ResoCollisions &collisions,
                         const aod::ResoV0s &v0s,
                         const aod::ResoCascades &cascades)
  {
    auto daughters = std::make_tuple(cascades, v0s);
    Pair<aod::ResoCollisions, aod::ResoCascades, aod::ResoV0s, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto &[collision1, cascades1, collision2, v0s2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillXiK0S<true>(collision1, cascades1, v0s2);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedXiK0S,
                 "Mixed-event Xi-K0S", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const &cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HeavyStrangeBaryonResonanceSearch>(cfgc)};
}
