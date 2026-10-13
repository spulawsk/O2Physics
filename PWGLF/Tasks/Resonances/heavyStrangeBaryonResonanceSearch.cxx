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
#include <Framework/AnalysisDataModel.h>
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
#include <cstdint>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace heavyresocand
{
DECLARE_SOA_COLUMN(CollisionId, collisionId, int64_t);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Rapidity, rapidity, float);
DECLARE_SOA_COLUMN(PairType, pairType, uint8_t);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, float);
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, float);

DECLARE_SOA_COLUMN(XiMass, xiMass, float);
DECLARE_SOA_COLUMN(XiPt, xiPt, float);
DECLARE_SOA_COLUMN(XiV0CosPA, xiV0CosPA, float);
DECLARE_SOA_COLUMN(XiCascCosPA, xiCascCosPA, float);
DECLARE_SOA_COLUMN(XiV0DaughterDCA, xiV0DaughterDCA, float);
DECLARE_SOA_COLUMN(XiCascDaughterDCA, xiCascDaughterDCA, float);
DECLARE_SOA_COLUMN(XiRadius, xiRadius, float);
DECLARE_SOA_COLUMN(XiBachelorTPCNSigmaPi, xiBachelorTPCNSigmaPi, float);

DECLARE_SOA_COLUMN(TrackPt, trackPt, float);
DECLARE_SOA_COLUMN(TrackEta, trackEta, float);
DECLARE_SOA_COLUMN(TrackDCAxy, trackDCAxy, float);
DECLARE_SOA_COLUMN(TrackDCAz, trackDCAz, float);
DECLARE_SOA_COLUMN(TrackTPCNCls, trackTPCNCls, int16_t);
DECLARE_SOA_COLUMN(TrackTPCNSigma, trackTPCNSigma, float);
DECLARE_SOA_COLUMN(TrackTOFNSigma, trackTOFNSigma, float);
DECLARE_SOA_COLUMN(TrackHasTOF, trackHasTOF, uint8_t);

DECLARE_SOA_COLUMN(V0Mass, v0Mass, float);
DECLARE_SOA_COLUMN(V0AlternativeMass, v0AlternativeMass, float);
DECLARE_SOA_COLUMN(V0SecondAlternativeMass, v0SecondAlternativeMass, float);
DECLARE_SOA_COLUMN(V0Pt, v0Pt, float);
DECLARE_SOA_COLUMN(V0CosPA, v0CosPA, float);
DECLARE_SOA_COLUMN(V0DaughterDCA, v0DaughterDCA, float);
DECLARE_SOA_COLUMN(V0Radius, v0Radius, float);
DECLARE_SOA_COLUMN(V0PosTPCNSigmaPr, v0PosTPCNSigmaPr, float);
DECLARE_SOA_COLUMN(V0NegTPCNSigmaPr, v0NegTPCNSigmaPr, float);
DECLARE_SOA_COLUMN(V0PosTPCNSigmaPi, v0PosTPCNSigmaPi, float);
DECLARE_SOA_COLUMN(V0NegTPCNSigmaPi, v0NegTPCNSigmaPi, float);
DECLARE_SOA_COLUMN(K0SMass, k0sMass, float);
DECLARE_SOA_COLUMN(K0SPt, k0sPt, float);
DECLARE_SOA_COLUMN(K0SCosPA, k0sCosPA, float);
DECLARE_SOA_COLUMN(K0SDaughterDCA, k0sDaughterDCA, float);
DECLARE_SOA_COLUMN(K0SRadius, k0sRadius, float);
DECLARE_SOA_COLUMN(K0SPosTPCNSigmaPi, k0sPosTPCNSigmaPi, float);
DECLARE_SOA_COLUMN(K0SNegTPCNSigmaPi, k0sNegTPCNSigmaPi, float);
} // namespace heavyresocand

DECLARE_SOA_TABLE(HeavyXiPiCandidates, "AOD", "HEAVYXIPI",
                  heavyresocand::CollisionId, heavyresocand::Mass,
                  heavyresocand::Pt, heavyresocand::Rapidity,
                  heavyresocand::PairType, heavyresocand::VertexZ,
                  heavyresocand::Multiplicity, heavyresocand::XiMass,
                  heavyresocand::XiPt, heavyresocand::XiV0CosPA,
                  heavyresocand::XiCascCosPA,
                  heavyresocand::XiV0DaughterDCA,
                  heavyresocand::XiCascDaughterDCA,
                  heavyresocand::XiRadius,
                  heavyresocand::XiBachelorTPCNSigmaPi,
                  heavyresocand::TrackPt, heavyresocand::TrackEta,
                  heavyresocand::TrackDCAxy, heavyresocand::TrackDCAz,
                  heavyresocand::TrackTPCNCls,
                  heavyresocand::TrackTPCNSigma,
                  heavyresocand::TrackTOFNSigma,
                  heavyresocand::TrackHasTOF);

DECLARE_SOA_TABLE(HeavyXiKCandidates, "AOD", "HEAVYXIK",
                  heavyresocand::CollisionId, heavyresocand::Mass,
                  heavyresocand::Pt, heavyresocand::Rapidity,
                  heavyresocand::PairType, heavyresocand::VertexZ,
                  heavyresocand::Multiplicity, heavyresocand::XiMass,
                  heavyresocand::XiPt, heavyresocand::XiV0CosPA,
                  heavyresocand::XiCascCosPA,
                  heavyresocand::XiV0DaughterDCA,
                  heavyresocand::XiCascDaughterDCA,
                  heavyresocand::XiRadius,
                  heavyresocand::XiBachelorTPCNSigmaPi,
                  heavyresocand::TrackPt, heavyresocand::TrackEta,
                  heavyresocand::TrackDCAxy, heavyresocand::TrackDCAz,
                  heavyresocand::TrackTPCNCls,
                  heavyresocand::TrackTPCNSigma,
                  heavyresocand::TrackTOFNSigma,
                  heavyresocand::TrackHasTOF);

DECLARE_SOA_TABLE(HeavyLambdaKCandidates, "AOD", "HEAVYLAMK",
                  heavyresocand::CollisionId, heavyresocand::Mass,
                  heavyresocand::Pt, heavyresocand::Rapidity,
                  heavyresocand::PairType, heavyresocand::VertexZ,
                  heavyresocand::Multiplicity, heavyresocand::V0Mass,
                  heavyresocand::V0AlternativeMass, heavyresocand::V0Pt,
                  heavyresocand::V0CosPA, heavyresocand::V0DaughterDCA,
                  heavyresocand::V0Radius,
                  heavyresocand::V0PosTPCNSigmaPr,
                  heavyresocand::V0NegTPCNSigmaPr,
                  heavyresocand::V0PosTPCNSigmaPi,
                  heavyresocand::V0NegTPCNSigmaPi,
                  heavyresocand::TrackPt, heavyresocand::TrackEta,
                  heavyresocand::TrackDCAxy, heavyresocand::TrackDCAz,
                  heavyresocand::TrackTPCNCls,
                  heavyresocand::TrackTPCNSigma,
                  heavyresocand::TrackTOFNSigma,
                  heavyresocand::TrackHasTOF);

DECLARE_SOA_TABLE(HeavyXiK0SCandidates, "AOD", "HEAVYXIK0S",
                  heavyresocand::CollisionId, heavyresocand::Mass,
                  heavyresocand::Pt, heavyresocand::Rapidity,
                  heavyresocand::PairType, heavyresocand::VertexZ,
                  heavyresocand::Multiplicity, heavyresocand::XiMass,
                  heavyresocand::XiPt, heavyresocand::XiV0CosPA,
                  heavyresocand::XiCascCosPA,
                  heavyresocand::XiV0DaughterDCA,
                  heavyresocand::XiCascDaughterDCA,
                  heavyresocand::XiRadius,
                  heavyresocand::XiBachelorTPCNSigmaPi,
                  heavyresocand::V0Mass,
                  heavyresocand::V0AlternativeMass,
                  heavyresocand::V0SecondAlternativeMass,
                  heavyresocand::V0Pt,
                  heavyresocand::V0CosPA, heavyresocand::V0DaughterDCA,
                  heavyresocand::V0Radius,
                  heavyresocand::V0PosTPCNSigmaPi,
                  heavyresocand::V0NegTPCNSigmaPi);

DECLARE_SOA_TABLE(HeavyLambdaPiCandidates, "AOD", "HEAVYLAMPI",
                  heavyresocand::CollisionId, heavyresocand::Mass,
                  heavyresocand::Pt, heavyresocand::Rapidity,
                  heavyresocand::PairType, heavyresocand::VertexZ,
                  heavyresocand::Multiplicity, heavyresocand::V0Mass,
                  heavyresocand::V0AlternativeMass, heavyresocand::V0Pt,
                  heavyresocand::V0CosPA, heavyresocand::V0DaughterDCA,
                  heavyresocand::V0Radius,
                  heavyresocand::V0PosTPCNSigmaPr,
                  heavyresocand::V0NegTPCNSigmaPr,
                  heavyresocand::V0PosTPCNSigmaPi,
                  heavyresocand::V0NegTPCNSigmaPi,
                  heavyresocand::TrackPt, heavyresocand::TrackEta,
                  heavyresocand::TrackDCAxy, heavyresocand::TrackDCAz,
                  heavyresocand::TrackTPCNCls,
                  heavyresocand::TrackTPCNSigma,
                  heavyresocand::TrackTOFNSigma,
                  heavyresocand::TrackHasTOF);

DECLARE_SOA_TABLE(HeavyLambdaK0SCandidates, "AOD", "HEAVYLAMK0S",
                  heavyresocand::CollisionId, heavyresocand::Mass,
                  heavyresocand::Pt, heavyresocand::Rapidity,
                  heavyresocand::PairType, heavyresocand::VertexZ,
                  heavyresocand::Multiplicity, heavyresocand::V0Mass,
                  heavyresocand::V0AlternativeMass, heavyresocand::V0Pt,
                  heavyresocand::V0CosPA, heavyresocand::V0DaughterDCA,
                  heavyresocand::V0Radius,
                  heavyresocand::V0PosTPCNSigmaPr,
                  heavyresocand::V0NegTPCNSigmaPr,
                  heavyresocand::V0PosTPCNSigmaPi,
                  heavyresocand::V0NegTPCNSigmaPi,
                  heavyresocand::K0SMass, heavyresocand::K0SPt,
                  heavyresocand::K0SCosPA,
                  heavyresocand::K0SDaughterDCA,
                  heavyresocand::K0SRadius,
                  heavyresocand::K0SPosTPCNSigmaPi,
                  heavyresocand::K0SNegTPCNSigmaPi);
} // namespace o2::aod

namespace
{
using FourVector = ROOT::Math::PxPyPzEVector;

template <typename T>
FourVector fourVector(const T& candidate, float mass)
{
  const auto energy = std::sqrt(candidate.px() * candidate.px() +
                                candidate.py() * candidate.py() +
                                candidate.pz() * candidate.pz() + mass * mass);
  return FourVector{candidate.px(), candidate.py(), candidate.pz(), energy};
}

template <typename TrackT, typename TrackIdsT>
int trackSourceId(const TrackT& track, const TrackIdsT& trackIds)
{
  const auto rowIndex = track.globalIndex();
  if (rowIndex >= 0 && rowIndex < trackIds.size()) {
    return trackIds.rawIteratorAt(rowIndex).trackId();
  }
  return -1;
}

template <typename TrackT, typename V0T, typename TrackIdsT>
bool sharesTrack(const TrackT& track, const V0T& v0,
                 const TrackIdsT& trackIds)
{
  const int sourceTrackId = trackSourceId(track, trackIds);
  const auto* indices = v0.indices();
  return sourceTrackId >= 0 &&
         (sourceTrackId == indices[0] || sourceTrackId == indices[1]);
}

template <typename TrackT, typename CascT, typename TrackIdsT>
bool sharesTrackWithCascade(const TrackT& track, const CascT& cascade,
                            const TrackIdsT& trackIds)
{
  const int sourceTrackId = trackSourceId(track, trackIds);
  const auto* indices = cascade.cascadeIndices();
  return sourceTrackId >= 0 &&
         (sourceTrackId == indices[0] || sourceTrackId == indices[1] ||
          sourceTrackId == indices[2]);
}

template <typename CascT, typename V0T>
bool sharesTrackCascadeV0(const CascT& cascade, const V0T& v0)
{
  const auto* cascadeIndices = cascade.cascadeIndices();
  const auto* v0Indices = v0.indices();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 2; ++j) {
      if (cascadeIndices[i] == v0Indices[j]) {
        return true;
      }
    }
  }
  return false;
}

template <typename V0T1, typename V0T2>
bool sharesTrackV0V0(const V0T1& first, const V0T2& second)
{
  const auto* firstIndices = first.indices();
  const auto* secondIndices = second.indices();
  return firstIndices[0] == secondIndices[0] ||
         firstIndices[0] == secondIndices[1] ||
         firstIndices[1] == secondIndices[0] ||
         firstIndices[1] == secondIndices[1];
}
} // namespace

struct HeavyStrangeBaryonResonanceSearch {
  Produces<aod::HeavyXiPiCandidates> xiPiCandidates;
  Produces<aod::HeavyXiKCandidates> xiKCandidates;
  Produces<aod::HeavyLambdaKCandidates> lambdaKCandidates;
  Produces<aod::HeavyXiK0SCandidates> xiK0sCandidates;
  Produces<aod::HeavyLambdaPiCandidates> lambdaPiCandidates;
  Produces<aod::HeavyLambdaK0SCandidates> lambdaK0sCandidates;

  SliceCache cache;
  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> requireINELgt0{"requireINELgt0", false,
                                    "Require reconstructed INEL>0"};
  Configurable<float> maxAbsRapidity{"maxAbsRapidity", 0.5f,
                                     "Maximum absolute candidate rapidity"};
  Configurable<bool> fillCandidateTables{"fillCandidateTables", true,
                                         "Write same-event candidate tables"};
  Configurable<float> tableXiMassWindow{"tableXiMassWindow", 0.03f,
                                        "Loose Xi mass window for tables"};
  Configurable<float> tableLambdaMassWindow{
    "tableLambdaMassWindow", 0.03f, "Loose Lambda mass window for tables"};
  Configurable<float> tableK0SMassWindow{"tableK0SMassWindow", 0.05f,
                                         "Loose K0S mass window for tables"};

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

  void init(InitContext const&)
  {
    const AxisSpec eventCounter{5, 0.5, 5.5, "selection step"};
    const AxisSpec vertexAxis{120, -12., 12., "z_{PV} (cm)"};
    const AxisSpec multAxis{110, 0., 110., "multiplicity percentile"};
    const AxisSpec ptAxis{200, 0., 20., "p_{T} (GeV/c)"};
    const AxisSpec rapidityAxis{40, -1., 1., "y"};
    const AxisSpec pairTypeAxis{4, -0.5, 3.5, "pair type"};
    const AxisSpec xiPiMassAxis{1250, 1.45, 2.70, "m_{#Xi#pi} (GeV/c^{2})"};
    const AxisSpec xiKMassAxis{1200, 1.80, 3.00, "m_{#Xi K} (GeV/c^{2})"};
    const AxisSpec lambdaKMassAxis{1091, 1.609, 2.70,
                                   "m_{#Lambda K} (GeV/c^{2})"};
    const AxisSpec xiK0sMassAxis{890, 1.81, 2.70,
                                 "m_{#Xi K^{0}_{S}} (GeV/c^{2})"};
    const AxisSpec lambdaPiMassAxis{1450, 1.25, 2.70,
                                    "m_{#Lambda#pi} (GeV/c^{2})"};
    const AxisSpec lambdaK0sMassAxis{1100, 1.60, 2.70,
                                     "m_{#Lambda K^{0}_{S}} (GeV/c^{2})"};

    registry.add("QA/lambdaMassBefore", "Lambda hypothesis before cuts",
                 HistType::kTH1F, {{400, 1.05, 1.25}});

    registry.add("QA/antiLambdaMassBefore", "Anti-Lambda hypothesis before cuts",
                 HistType::kTH1F, {{400, 1.05, 1.25}});

    registry.add("QA/lambdaMassAfter", "Lambda hypothesis after cuts",
                 HistType::kTH1F, {{200, 1.08, 1.16}});

    registry.add("QA/pionTPCBefore", "TPC pion PID before cuts",
                 HistType::kTH2F, {ptAxis, {200, -10., 10.}});

    registry.add("QA/kaonTPCBefore", "TPC kaon PID before cuts",
                 HistType::kTH2F, {ptAxis, {200, -10., 10.}});

    registry.add("QA/pionTPCAfter", "TPC pion PID after track and PID cuts",
                 HistType::kTH2F, {ptAxis, {100, -5., 5.}});

    registry.add("QA/kaonTPCAfter", "TPC kaon PID after track and PID cuts",
                 HistType::kTH2F, {ptAxis, {100, -5., 5.}});

    registry.add("QA/pionTPCSelected", "TPC pion PID after complete PID selection",
                 HistType::kTH2F, {ptAxis, {100, -5., 5.}});
    registry.add("QA/kaonTPCSelected", "TPC kaon PID after complete PID selection",
                 HistType::kTH2F, {ptAxis, {100, -5., 5.}});
    registry.add("QA/trackCutFlow", "Primary-track selection cut flow",
                 HistType::kTH1F, {{7, 0.5, 7.5}});
    registry.add("QA/lambdaCutFlow", "Lambda selection cut flow",
                 HistType::kTH1F, {{7, 0.5, 7.5}});

    registry.add("Event/nTracks", "ResoTracks per collision",
                 HistType::kTH1F, {{501, -0.5, 500.5}});

    registry.add("Event/nV0s", "ResoV0s per collision",
                 HistType::kTH1F, {{101, -0.5, 100.5}});

    registry.add("Event/nCascades", "ResoCascades per collision",
                 HistType::kTH1F, {{51, -0.5, 50.5}});

    registry.add("Event/counter", "Event counter", HistType::kTH1F,
                 {eventCounter});
    registry.add("Event/vertexZ", "Accepted-event vertex", HistType::kTH1F,
                 {vertexAxis});
    registry.add("Event/multiplicity", "Accepted-event multiplicity",
                 HistType::kTH1F, {multAxis});

    registry.add("QA/xiMass", "Selected Xi candidates", HistType::kTH1F,
                 {{240, 1.27, 1.39}});
    registry.add("QA/k0sMass", "Selected K0S candidates", HistType::kTH1F,
                 {{200, 0.44, 0.56}});

    registry.add("XiPi/same", "Xi-pion same event", HistType::kTHnSparseF,
                 {xiPiMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
    registry.add("XiPi/mixed", "Xi-pion mixed event", HistType::kTHnSparseF,
                 {xiPiMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
    registry.add("XiK/same", "Xi-kaon same event", HistType::kTHnSparseF,
                 {xiKMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
    registry.add("XiK/mixed", "Xi-kaon mixed event", HistType::kTHnSparseF,
                 {xiKMassAxis, ptAxis, rapidityAxis, multAxis, pairTypeAxis});
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
    registry.add("LambdaPi/same", "Lambda-pion same event",
                 HistType::kTHnSparseF,
                 {lambdaPiMassAxis, ptAxis, rapidityAxis, multAxis,
                  pairTypeAxis});
    registry.add("LambdaPi/mixed", "Lambda-pion mixed event",
                 HistType::kTHnSparseF,
                 {lambdaPiMassAxis, ptAxis, rapidityAxis, multAxis,
                  pairTypeAxis});
    registry.add("LambdaK0S/same", "Lambda-K0S same event",
                 HistType::kTHnSparseF,
                 {lambdaK0sMassAxis, ptAxis, rapidityAxis, multAxis,
                  pairTypeAxis});
    registry.add("LambdaK0S/mixed", "Lambda-K0S mixed event",
                 HistType::kTHnSparseF,
                 {lambdaK0sMassAxis, ptAxis, rapidityAxis, multAxis,
                  pairTypeAxis});
  }

  template <typename CollisionT>
  bool selectEvent(const CollisionT& collision, bool fillQA)
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
  bool selectPrimaryTrack(const TrackT& track)
  {
    return track.pt() >= minTrackPt && std::abs(track.eta()) <= maxTrackEta &&
           std::abs(track.dcaXY()) <= maxTrackDCAxy &&
           std::abs(track.dcaZ()) <= maxTrackDCAz &&
           track.tpcNClsFound() >= minTPCClusters;
  }

  template <typename TrackT>
  bool selectPrimaryTrackWithQA(const TrackT& track)
  {
    registry.fill(HIST("QA/trackCutFlow"), 1.);
    if (track.pt() < minTrackPt) {
      return false;
    }
    registry.fill(HIST("QA/trackCutFlow"), 2.);
    if (std::abs(track.eta()) > maxTrackEta) {
      return false;
    }
    registry.fill(HIST("QA/trackCutFlow"), 3.);
    if (std::abs(track.dcaXY()) > maxTrackDCAxy) {
      return false;
    }
    registry.fill(HIST("QA/trackCutFlow"), 4.);
    if (std::abs(track.dcaZ()) > maxTrackDCAz) {
      return false;
    }
    registry.fill(HIST("QA/trackCutFlow"), 5.);
    if (track.tpcNClsFound() < minTPCClusters) {
      return false;
    }
    registry.fill(HIST("QA/trackCutFlow"), 6.);
    registry.fill(HIST("QA/trackCutFlow"), 7.);
    return true;
  }

  template <typename V0T>
  void fillLambdaCutFlow(const V0T& v0, bool antiLambda)
  {
    const float mass = antiLambda ? v0.mAntiLambda() : v0.mLambda();
    registry.fill(HIST("QA/lambdaCutFlow"), 1.);
    if (std::abs(mass - MassLambda0) > lambdaMassWindow) {
      return;
    }
    registry.fill(HIST("QA/lambdaCutFlow"), 2.);
    if (v0.v0CosPA() < lambdaMinCosPA) {
      return;
    }
    registry.fill(HIST("QA/lambdaCutFlow"), 3.);
    if (v0.daughDCA() > lambdaMaxDaughterDCA) {
      return;
    }
    registry.fill(HIST("QA/lambdaCutFlow"), 4.);
    if (v0.transRadius() < lambdaMinRadius ||
        v0.transRadius() > lambdaMaxRadius) {
      return;
    }
    registry.fill(HIST("QA/lambdaCutFlow"), 5.);
    const bool daughterPID =
      antiLambda
        ? (std::abs(v0.daughterTPCNSigmaPosPi()) <= lambdaDaughterPID &&
           std::abs(v0.daughterTPCNSigmaNegPr()) <= lambdaDaughterPID)
        : (std::abs(v0.daughterTPCNSigmaPosPr()) <= lambdaDaughterPID &&
           std::abs(v0.daughterTPCNSigmaNegPi()) <= lambdaDaughterPID);
    if (!daughterPID) {
      return;
    }
    registry.fill(HIST("QA/lambdaCutFlow"), 6.);
    registry.fill(HIST("QA/lambdaMassAfter"), mass);
  }

  template <typename TrackT>
  bool selectPion(const TrackT& track)
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
  bool selectKaon(const TrackT& track)
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
  bool selectLambda(const V0T& v0, bool antiLambda)
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
  bool selectXi(const CascT& xi)
  {
    return std::abs(xi.mXi() - MassXiMinus) <= xiMassWindow &&
           xi.v0CosPA() >= xiMinV0CosPA && xi.cascCosPA() >= xiMinCascCosPA &&
           xi.daughDCA() <= xiMaxV0DaughterDCA &&
           xi.cascDaughDCA() <= xiMaxCascDaughterDCA &&
           xi.cascTransRadius() >= xiMinRadius &&
           std::abs(xi.daughterTPCNSigmaBachPi()) <= xiBachelorPID;
  }

  template <typename V0T>
  bool selectK0S(const V0T& k0s)
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

  template <typename CollisionT, typename CascadesT, typename TracksT,
            typename TrackIdsT>
  void writeXiPiCandidates(const CollisionT& collision,
                           const CascadesT& cascades, const TracksT& tracks,
                           const TrackIdsT& trackIds)
  {
    for (const auto& xi : cascades) {
      if (std::abs(xi.mXi() - MassXiMinus) > tableXiMassWindow) {
        continue;
      }
      const auto pXi = fourVector(xi, xi.mXi());
      for (const auto& pion : tracks) {
        if (!selectPrimaryTrack(pion) ||
            sharesTrackWithCascade(pion, xi, trackIds)) {
          continue;
        }
        const auto candidate = pXi + fourVector(pion, MassPiPlus);
        const uint8_t pairType =
          xi.sign() < 0 ? (pion.sign() > 0 ? 0 : 2)
                        : (pion.sign() < 0 ? 1 : 3);
        xiPiCandidates(
          collision.globalIndex(), candidate.M(), candidate.Pt(),
          candidate.Rapidity(), pairType, collision.posZ(), collision.cent(),
          xi.mXi(), xi.pt(), xi.v0CosPA(), xi.cascCosPA(), xi.daughDCA(),
          xi.cascDaughDCA(), xi.cascTransRadius(),
          xi.daughterTPCNSigmaBachPi(), pion.pt(), pion.eta(), pion.dcaXY(),
          pion.dcaZ(), pion.tpcNClsFound(), pion.tpcNSigmaPi(),
          pion.hasTOF() ? pion.tofNSigmaPi() : 999.f,
          static_cast<uint8_t>(pion.hasTOF()));
      }
    }
  }

  template <typename CollisionT, typename CascadesT, typename TracksT,
            typename TrackIdsT>
  void writeXiKCandidates(const CollisionT& collision,
                          const CascadesT& cascades, const TracksT& tracks,
                          const TrackIdsT& trackIds)
  {
    for (const auto& xi : cascades) {
      if (std::abs(xi.mXi() - MassXiMinus) > tableXiMassWindow) {
        continue;
      }
      const auto pXi = fourVector(xi, xi.mXi());
      for (const auto& kaon : tracks) {
        if (!selectPrimaryTrack(kaon) ||
            sharesTrackWithCascade(kaon, xi, trackIds)) {
          continue;
        }
        const auto candidate = pXi + fourVector(kaon, MassKPlus);
        const uint8_t pairType =
          xi.sign() < 0 ? (kaon.sign() > 0 ? 0 : 2)
                        : (kaon.sign() < 0 ? 1 : 3);
        xiKCandidates(
          collision.globalIndex(), candidate.M(), candidate.Pt(),
          candidate.Rapidity(), pairType, collision.posZ(), collision.cent(),
          xi.mXi(), xi.pt(), xi.v0CosPA(), xi.cascCosPA(), xi.daughDCA(),
          xi.cascDaughDCA(), xi.cascTransRadius(),
          xi.daughterTPCNSigmaBachPi(), kaon.pt(), kaon.eta(), kaon.dcaXY(),
          kaon.dcaZ(), kaon.tpcNClsFound(), kaon.tpcNSigmaKa(),
          kaon.hasTOF() ? kaon.tofNSigmaKa() : 999.f,
          static_cast<uint8_t>(kaon.hasTOF()));
      }
    }
  }

  template <typename CollisionT, typename V0sT, typename TracksT,
            typename TrackIdsT>
  void writeLambdaKCandidates(const CollisionT& collision, const V0sT& v0s,
                              const TracksT& tracks,
                              const TrackIdsT& trackIds)
  {
    for (const auto& lambda : v0s) {
      for (uint8_t anti = 0; anti < 2; ++anti) {
        const float lambdaMass = anti ? lambda.mAntiLambda() : lambda.mLambda();
        if (std::abs(lambdaMass - MassLambda0) > tableLambdaMassWindow) {
          continue;
        }
        const float alternativeMass =
          anti ? lambda.mLambda() : lambda.mAntiLambda();
        const auto pLambda = fourVector(lambda, lambdaMass);
        for (const auto& kaon : tracks) {
          if (!selectPrimaryTrack(kaon) || sharesTrack(kaon, lambda, trackIds)) {
            continue;
          }
          const auto candidate = pLambda + fourVector(kaon, MassKPlus);
          const uint8_t pairType =
            anti == 0 ? (kaon.sign() < 0 ? 0 : 2)
                      : (kaon.sign() > 0 ? 1 : 3);
          lambdaKCandidates(
            collision.globalIndex(), candidate.M(), candidate.Pt(),
            candidate.Rapidity(), pairType, collision.posZ(), collision.cent(),
            lambdaMass, alternativeMass, lambda.pt(), lambda.v0CosPA(),
            lambda.daughDCA(), lambda.transRadius(),
            lambda.daughterTPCNSigmaPosPr(),
            lambda.daughterTPCNSigmaNegPr(),
            lambda.daughterTPCNSigmaPosPi(),
            lambda.daughterTPCNSigmaNegPi(), kaon.pt(), kaon.eta(),
            kaon.dcaXY(), kaon.dcaZ(), kaon.tpcNClsFound(),
            kaon.tpcNSigmaKa(),
            kaon.hasTOF() ? kaon.tofNSigmaKa() : 999.f,
            static_cast<uint8_t>(kaon.hasTOF()));
        }
      }
    }
  }

  template <typename CollisionT, typename CascadesT, typename V0sT>
  void writeXiK0SCandidates(const CollisionT& collision,
                            const CascadesT& cascades, const V0sT& v0s)
  {
    for (const auto& xi : cascades) {
      if (std::abs(xi.mXi() - MassXiMinus) > tableXiMassWindow) {
        continue;
      }
      const auto pXi = fourVector(xi, xi.mXi());
      for (const auto& k0s : v0s) {
        if (std::abs(k0s.mK0Short() - MassK0Short) > tableK0SMassWindow ||
            sharesTrackCascadeV0(xi, k0s)) {
          continue;
        }
        const auto candidate = pXi + fourVector(k0s, MassK0Short);
        const uint8_t pairType = xi.sign() < 0 ? 0 : 1;
        xiK0sCandidates(
          collision.globalIndex(), candidate.M(), candidate.Pt(),
          candidate.Rapidity(), pairType, collision.posZ(), collision.cent(),
          xi.mXi(), xi.pt(), xi.v0CosPA(), xi.cascCosPA(), xi.daughDCA(),
          xi.cascDaughDCA(), xi.cascTransRadius(),
          xi.daughterTPCNSigmaBachPi(), k0s.mK0Short(), k0s.mLambda(),
          k0s.mAntiLambda(), k0s.pt(), k0s.v0CosPA(), k0s.daughDCA(),
          k0s.transRadius(), k0s.daughterTPCNSigmaPosPi(),
          k0s.daughterTPCNSigmaNegPi());
      }
    }
  }

  template <typename CollisionT, typename V0sT, typename TracksT,
            typename TrackIdsT>
  void writeLambdaPiCandidates(const CollisionT& collision, const V0sT& v0s,
                               const TracksT& tracks,
                               const TrackIdsT& trackIds)
  {
    for (const auto& lambda : v0s) {
      for (uint8_t anti = 0; anti < 2; ++anti) {
        const float lambdaMass = anti ? lambda.mAntiLambda() : lambda.mLambda();
        if (std::abs(lambdaMass - MassLambda0) > tableLambdaMassWindow) {
          continue;
        }
        const float alternativeMass =
          anti ? lambda.mLambda() : lambda.mAntiLambda();
        const auto pLambda = fourVector(lambda, lambdaMass);
        for (const auto& pion : tracks) {
          if (!selectPrimaryTrack(pion) || sharesTrack(pion, lambda, trackIds)) {
            continue;
          }
          const auto candidate = pLambda + fourVector(pion, MassPiPlus);
          const uint8_t pairType =
            anti == 0 ? (pion.sign() > 0 ? 0 : 2)
                      : (pion.sign() < 0 ? 1 : 3);
          lambdaPiCandidates(
            collision.globalIndex(), candidate.M(), candidate.Pt(),
            candidate.Rapidity(), pairType, collision.posZ(), collision.cent(),
            lambdaMass, alternativeMass, lambda.pt(), lambda.v0CosPA(),
            lambda.daughDCA(), lambda.transRadius(),
            lambda.daughterTPCNSigmaPosPr(),
            lambda.daughterTPCNSigmaNegPr(),
            lambda.daughterTPCNSigmaPosPi(),
            lambda.daughterTPCNSigmaNegPi(), pion.pt(), pion.eta(),
            pion.dcaXY(), pion.dcaZ(), pion.tpcNClsFound(),
            pion.tpcNSigmaPi(),
            pion.hasTOF() ? pion.tofNSigmaPi() : 999.f,
            static_cast<uint8_t>(pion.hasTOF()));
        }
      }
    }
  }

  template <typename CollisionT, typename V0sT>
  void writeLambdaK0SCandidates(const CollisionT& collision,
                                const V0sT& v0s)
  {
    for (const auto& lambda : v0s) {
      for (uint8_t anti = 0; anti < 2; ++anti) {
        const float lambdaMass = anti ? lambda.mAntiLambda() : lambda.mLambda();
        if (std::abs(lambdaMass - MassLambda0) > tableLambdaMassWindow) {
          continue;
        }
        const float alternativeMass =
          anti ? lambda.mLambda() : lambda.mAntiLambda();
        const auto pLambda = fourVector(lambda, lambdaMass);
        for (const auto& k0s : v0s) {
          if (std::abs(k0s.mK0Short() - MassK0Short) > tableK0SMassWindow ||
              sharesTrackV0V0(lambda, k0s)) {
            continue;
          }
          const auto candidate = pLambda + fourVector(k0s, MassK0Short);
          lambdaK0sCandidates(
            collision.globalIndex(), candidate.M(), candidate.Pt(),
            candidate.Rapidity(), anti, collision.posZ(), collision.cent(),
            lambdaMass, alternativeMass, lambda.pt(), lambda.v0CosPA(),
            lambda.daughDCA(), lambda.transRadius(),
            lambda.daughterTPCNSigmaPosPr(),
            lambda.daughterTPCNSigmaNegPr(),
            lambda.daughterTPCNSigmaPosPi(),
            lambda.daughterTPCNSigmaNegPi(), k0s.mK0Short(), k0s.pt(),
            k0s.v0CosPA(), k0s.daughDCA(), k0s.transRadius(),
            k0s.daughterTPCNSigmaPosPi(),
            k0s.daughterTPCNSigmaNegPi());
        }
      }
    }
  }

  template <bool Mixed, typename CollisionT, typename CascadesT,
            typename TracksT, typename TrackIdsT>
  void fillXiPi(const CollisionT& collision, const CascadesT& cascades,
                const TracksT& tracks, const TrackIdsT& trackIds)
  {
    for (const auto& xi : cascades) {
      if (!selectXi(xi)) {
        continue;
      }
      if constexpr (!Mixed) {
        registry.fill(HIST("QA/xiMass"), xi.mXi());
      }
      const auto pXi = fourVector(xi, xi.mXi());
      for (const auto& pion : tracks) {
        if (!selectPion(pion)) {
          continue;
        }
        if constexpr (!Mixed) {
          if (sharesTrackWithCascade(pion, xi, trackIds)) {
            continue;
          }
        }
        const auto candidate = pXi + fourVector(pion, MassPiPlus);
        if (std::abs(candidate.Rapidity()) > maxAbsRapidity) {
          continue;
        }
        // 0: Xi- pi+ signal, 1: anti-Xi+ pi- signal, 2/3: like sign.
        const int pairType = xi.sign() < 0 ? (pion.sign() > 0 ? 0 : 2)
                                           : (pion.sign() < 0 ? 1 : 3);
        if constexpr (Mixed) {
          registry.fill(HIST("XiPi/mixed"), candidate.M(), candidate.Pt(),
                        candidate.Rapidity(), collision.cent(), pairType);
        } else {
          registry.fill(HIST("XiPi/same"), candidate.M(), candidate.Pt(),
                        candidate.Rapidity(), collision.cent(), pairType);
        }
      }
    }
  }

  template <bool Mixed, typename CollisionT, typename CascadesT,
            typename TracksT, typename TrackIdsT>
  void fillXiK(const CollisionT& collision, const CascadesT& cascades,
               const TracksT& tracks, const TrackIdsT& trackIds)
  {
    for (const auto& xi : cascades) {
      if (!selectXi(xi)) {
        continue;
      }
      const auto pXi = fourVector(xi, xi.mXi());
      for (const auto& kaon : tracks) {
        if (!selectKaon(kaon)) {
          continue;
        }
        if constexpr (!Mixed) {
          if (sharesTrackWithCascade(kaon, xi, trackIds)) {
            continue;
          }
        }
        const auto candidate = pXi + fourVector(kaon, MassKPlus);
        if (std::abs(candidate.Rapidity()) > maxAbsRapidity) {
          continue;
        }
        const int pairType =
          xi.sign() < 0 ? (kaon.sign() > 0 ? 0 : 2)
                        : (kaon.sign() < 0 ? 1 : 3);
        if constexpr (Mixed) {
          registry.fill(HIST("XiK/mixed"), candidate.M(), candidate.Pt(),
                        candidate.Rapidity(), collision.cent(), pairType);
        } else {
          registry.fill(HIST("XiK/same"), candidate.M(), candidate.Pt(),
                        candidate.Rapidity(), collision.cent(), pairType);
        }
      }
    }
  }

  template <bool Mixed, typename CollisionT, typename V0sT, typename TracksT,
            typename TrackIdsT>
  void fillLambdaK(const CollisionT& collision, const V0sT& v0s,
                   const TracksT& tracks, const TrackIdsT& trackIds)
  {
    for (const auto& lambda : v0s) {
      const bool isLambda = selectLambda(lambda, false);
      const bool isAntiLambda = selectLambda(lambda, true);
      if (!isLambda && !isAntiLambda) {
        continue;
      }
      for (const auto& kaon : tracks) {
        if (!selectKaon(kaon)) {
          continue;
        }
        if constexpr (!Mixed) {
          if (sharesTrack(kaon, lambda, trackIds)) {
            continue;
          }
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
          if constexpr (Mixed) {
            registry.fill(HIST("LambdaK/mixed"), candidate.M(), candidate.Pt(),
                          candidate.Rapidity(), collision.cent(), pairType);
          } else {
            registry.fill(HIST("LambdaK/same"), candidate.M(), candidate.Pt(),
                          candidate.Rapidity(), collision.cent(), pairType);
          }
        }
      }
    }
  }

  template <bool Mixed, typename CollisionT, typename V0sT, typename TracksT,
            typename TrackIdsT>
  void fillLambdaPi(const CollisionT& collision, const V0sT& v0s,
                    const TracksT& tracks, const TrackIdsT& trackIds)
  {
    for (const auto& lambda : v0s) {
      const bool isLambda = selectLambda(lambda, false);
      const bool isAntiLambda = selectLambda(lambda, true);
      if (!isLambda && !isAntiLambda) {
        continue;
      }
      for (const auto& pion : tracks) {
        if (!selectPion(pion)) {
          continue;
        }
        if constexpr (!Mixed) {
          if (sharesTrack(pion, lambda, trackIds)) {
            continue;
          }
        }
        for (int anti = 0; anti < 2; ++anti) {
          if ((anti == 0 && !isLambda) || (anti == 1 && !isAntiLambda)) {
            continue;
          }
          const float lambdaMass =
            anti ? lambda.mAntiLambda() : lambda.mLambda();
          const auto candidate =
            fourVector(lambda, lambdaMass) + fourVector(pion, MassPiPlus);
          if (std::abs(candidate.Rapidity()) > maxAbsRapidity) {
            continue;
          }
          const int pairType =
            anti == 0 ? (pion.sign() > 0 ? 0 : 2)
                      : (pion.sign() < 0 ? 1 : 3);
          if constexpr (Mixed) {
            registry.fill(HIST("LambdaPi/mixed"), candidate.M(),
                          candidate.Pt(), candidate.Rapidity(),
                          collision.cent(), pairType);
          } else {
            registry.fill(HIST("LambdaPi/same"), candidate.M(),
                          candidate.Pt(), candidate.Rapidity(),
                          collision.cent(), pairType);
          }
        }
      }
    }
  }

  template <bool Mixed, typename CollisionT, typename V0sT1, typename V0sT2>
  void fillLambdaK0S(const CollisionT& collision, const V0sT1& lambdas,
                     const V0sT2& k0sCandidates)
  {
    for (const auto& lambda : lambdas) {
      const bool isLambda = selectLambda(lambda, false);
      const bool isAntiLambda = selectLambda(lambda, true);
      if (!isLambda && !isAntiLambda) {
        continue;
      }
      for (const auto& k0s : k0sCandidates) {
        if (!selectK0S(k0s)) {
          continue;
        }
        if constexpr (!Mixed) {
          if (sharesTrackV0V0(lambda, k0s)) {
            continue;
          }
        }
        for (int anti = 0; anti < 2; ++anti) {
          if ((anti == 0 && !isLambda) || (anti == 1 && !isAntiLambda)) {
            continue;
          }
          const float lambdaMass =
            anti ? lambda.mAntiLambda() : lambda.mLambda();
          const auto candidate =
            fourVector(lambda, lambdaMass) + fourVector(k0s, MassK0Short);
          if (std::abs(candidate.Rapidity()) > maxAbsRapidity) {
            continue;
          }
          if constexpr (Mixed) {
            registry.fill(HIST("LambdaK0S/mixed"), candidate.M(),
                          candidate.Pt(), candidate.Rapidity(),
                          collision.cent(), anti);
          } else {
            registry.fill(HIST("LambdaK0S/same"), candidate.M(),
                          candidate.Pt(), candidate.Rapidity(),
                          collision.cent(), anti);
          }
        }
      }
    }
  }

  template <bool Mixed, typename CollisionT, typename CascadesT, typename V0sT>
  void fillXiK0S(const CollisionT& collision, const CascadesT& cascades,
                 const V0sT& v0s)
  {
    for (const auto& xi : cascades) {
      if (!selectXi(xi)) {
        continue;
      }
      for (const auto& k0s : v0s) {
        if (!selectK0S(k0s) || (!Mixed && sharesTrackCascadeV0(xi, k0s))) {
          continue;
        }
        const auto candidate =
          fourVector(xi, xi.mXi()) + fourVector(k0s, MassK0Short);
        if (std::abs(candidate.Rapidity()) > maxAbsRapidity) {
          continue;
        }
        const int pairType = xi.sign() < 0 ? 0 : 1;
        if constexpr (Mixed) {
          registry.fill(HIST("XiK0S/mixed"), candidate.M(), candidate.Pt(),
                        candidate.Rapidity(), collision.cent(), pairType);
        } else {
          registry.fill(HIST("XiK0S/same"), candidate.M(), candidate.Pt(),
                        candidate.Rapidity(), collision.cent(), pairType);
        }
        if constexpr (!Mixed) {
          registry.fill(HIST("QA/k0sMass"), k0s.mK0Short());
        }
      }
    }
  }

  void processData(const aod::ResoCollision& collision,
                   const aod::ResoTracks& tracks,
                   const aod::ResoTrackTracks& trackIds,
                   const aod::ResoV0s& v0s,
                   const aod::ResoCascades& cascades)
  {
    registry.fill(HIST("Event/nTracks"), tracks.size());
    registry.fill(HIST("Event/nV0s"), v0s.size());
    registry.fill(HIST("Event/nCascades"), cascades.size());

    // Input and selection QA are intentionally independent of pair building.
    for (const auto& track : tracks) {
      registry.fill(HIST("QA/pionTPCBefore"), track.pt(), track.tpcNSigmaPi());
      registry.fill(HIST("QA/kaonTPCBefore"), track.pt(), track.tpcNSigmaKa());

      if (!selectPrimaryTrackWithQA(track)) {
        continue;
      }

      // "After" monitors primary-track + TPC PID only.  The additional
      // "Selected" histograms include the TOF requirement as well.
      if (std::abs(track.tpcNSigmaPi()) <= maxTPCNSigmaPion) {
        registry.fill(HIST("QA/pionTPCAfter"), track.pt(), track.tpcNSigmaPi());
      }
      if (std::abs(track.tpcNSigmaKa()) <= maxTPCNSigmaKaon) {
        registry.fill(HIST("QA/kaonTPCAfter"), track.pt(), track.tpcNSigmaKa());
      }
      if (selectPion(track)) {
        registry.fill(HIST("QA/pionTPCSelected"), track.pt(), track.tpcNSigmaPi());
      }
      if (selectKaon(track)) {
        registry.fill(HIST("QA/kaonTPCSelected"), track.pt(), track.tpcNSigmaKa());
      }
    }

    for (const auto& v0 : v0s) {
      registry.fill(HIST("QA/lambdaMassBefore"), v0.mLambda());
      registry.fill(HIST("QA/antiLambdaMassBefore"), v0.mAntiLambda());
      fillLambdaCutFlow(v0, false);
      fillLambdaCutFlow(v0, true);
    }

    if (!selectEvent(collision, true)) {
      return;
    }

    if (fillCandidateTables) {
      writeXiPiCandidates(collision, cascades, tracks, trackIds);
      writeXiKCandidates(collision, cascades, tracks, trackIds);
      writeLambdaKCandidates(collision, v0s, tracks, trackIds);
      writeXiK0SCandidates(collision, cascades, v0s);
      writeLambdaPiCandidates(collision, v0s, tracks, trackIds);
      writeLambdaK0SCandidates(collision, v0s);
    }

    fillXiPi<false>(collision, cascades, tracks, trackIds);
    fillXiK<false>(collision, cascades, tracks, trackIds);
    fillLambdaK<false>(collision, v0s, tracks, trackIds);
    fillXiK0S<false>(collision, cascades, v0s);
    fillLambdaPi<false>(collision, v0s, tracks, trackIds);
    fillLambdaK0S<false>(collision, v0s, v0s);
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processData,
                 "Process same-event data", true);

  void processMixedXiPi(const aod::ResoCollisions& collisions,
                        const aod::ResoTracks& tracks,
                        const aod::ResoTrackTracks& trackIds,
                        const aod::ResoCascades& cascades)
  {
    auto daughters = std::make_tuple(cascades, tracks);
    Pair<aod::ResoCollisions, aod::ResoCascades, aod::ResoTracks, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto& [collision1, cascades1, collision2, tracks2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillXiPi<true>(collision1, cascades1, tracks2, trackIds);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedXiPi,
                 "Mixed-event Xi-pion", false);

  void processMixedXiK(const aod::ResoCollisions& collisions,
                       const aod::ResoTracks& tracks,
                       const aod::ResoTrackTracks& trackIds,
                       const aod::ResoCascades& cascades)
  {
    auto daughters = std::make_tuple(cascades, tracks);
    Pair<aod::ResoCollisions, aod::ResoCascades, aod::ResoTracks, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto& [collision1, cascades1, collision2, tracks2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillXiK<true>(collision1, cascades1, tracks2, trackIds);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedXiK,
                 "Mixed-event Xi-kaon", false);

  void processMixedLambdaK(const aod::ResoCollisions& collisions,
                           const aod::ResoTracks& tracks,
                           const aod::ResoTrackTracks& trackIds,
                           const aod::ResoV0s& v0s)
  {
    auto daughters = std::make_tuple(v0s, tracks);
    Pair<aod::ResoCollisions, aod::ResoV0s, aod::ResoTracks, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto& [collision1, v0s1, collision2, tracks2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillLambdaK<true>(collision1, v0s1, tracks2, trackIds);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedLambdaK,
                 "Mixed-event Lambda-kaon", false);

  void processMixedXiK0S(const aod::ResoCollisions& collisions,
                         const aod::ResoV0s& v0s,
                         const aod::ResoCascades& cascades)
  {
    auto daughters = std::make_tuple(cascades, v0s);
    Pair<aod::ResoCollisions, aod::ResoCascades, aod::ResoV0s, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto& [collision1, cascades1, collision2, v0s2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillXiK0S<true>(collision1, cascades1, v0s2);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedXiK0S,
                 "Mixed-event Xi-K0S", false);

  void processMixedLambdaPi(const aod::ResoCollisions& collisions,
                            const aod::ResoTracks& tracks,
                            const aod::ResoTrackTracks& trackIds,
                            const aod::ResoV0s& v0s)
  {
    auto daughters = std::make_tuple(v0s, tracks);
    Pair<aod::ResoCollisions, aod::ResoV0s, aod::ResoTracks, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto& [collision1, v0s1, collision2, tracks2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillLambdaPi<true>(collision1, v0s1, tracks2, trackIds);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedLambdaPi,
                 "Mixed-event Lambda-pion", false);

  void processMixedLambdaK0S(const aod::ResoCollisions& collisions,
                             const aod::ResoV0s& v0s)
  {
    auto daughters = std::make_tuple(v0s);
    SameKindPair<aod::ResoCollisions, aod::ResoV0s, MixingBinning>
      pairs{mixingBinning, numberMixedEvents, -1,
            collisions, daughters, &cache};
    for (const auto& [collision1, v0s1, collision2, v0s2] : pairs) {
      if (!selectEvent(collision1, false) || !selectEvent(collision2, false)) {
        continue;
      }
      fillLambdaK0S<true>(collision1, v0s1, v0s2);
    }
  }
  PROCESS_SWITCH(HeavyStrangeBaryonResonanceSearch, processMixedLambdaK0S,
                 "Mixed-event Lambda-K0S", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HeavyStrangeBaryonResonanceSearch>(cfgc)};
}
