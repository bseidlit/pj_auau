"""Raw release fixtures with full identities and deliberately interleaved rows."""
from copy import deepcopy
import ROOT
from fixtures import event_schema, photon_schema, truth_schema, ev, ph, tr, write_tree

EVENT = dict(event_schema, event_sequence="L", physical_event_sequence="L")
PHOTON = dict(EVENT, **{k: v for k, v in photon_schema.items() if k not in EVENT},
              candidate_id_hi="l", candidate_id_lo="l")
TRUTH = dict(truth_schema, event_id_hi="l", truth_photon_id_hi="l", truth_photon_id_lo="l")
LINK = dict(source_file_index="I", event_id_hi="l", event_id_lo="l",
            reco_type="I", truth_type="I", link_class="I", reco_index="I", truth_index="I",
            reco_id_hi="l", reco_id_lo="l", truth_id_hi="l", truth_id_lo="l")
JET = dict(source_file_index="I", event_id_hi="l", event_id_lo="l", truth_jet_pt="D")
KEY = ("source_file_index", "event_id_hi", "event_id_lo")


def fixture_rows():
    events = [dict(ev, source_entry=i, event_id_lo=i + 1, event_weight=1. if i == 3 else 2.,
                   total_calo_energy=float("nan"), event_sequence=i, physical_event_sequence=10 + i)
              for i in range(4)]
    # Event 3 has the same low ID/source as event 0, but a different high ID.
    events[3].update(event_id_hi=3, event_id_lo=1, vertex_z=45., terminal_status=2)
    photons = []
    for e, cid, et, ordinal in [(0, 101, 15., 9), (2, 103, 36., 0), (0, 102, 20., 14)]:
        row = dict(ph, **events[e])
        row.update(candidate_id_hi=11, candidate_id_lo=cid, photon_et=et,
                   photon_encounter_ordinal=ordinal, bdt_input_13=float("nan"))
        photons.append(row)
    photons[2].update(bdt_score=.9, iso_r04=5., bdt_is_tight=1, iso_r04_pass=0)
    truths = [dict(tr, **{key: events[e][key] for key in KEY}, truth_photon_id_hi=12,
                   truth_photon_id_lo=tid, generator_barcode=42 + e,
                   prompt_class=2 if e == 1 else 1, truth_photon_eta=.7 if e == 2 else 0.)
              for e, tid in [(0, 201), (1, 202), (2, 203)]]
    links = [dict({key: events[e][key] for key in KEY}, reco_type=1, truth_type=1, link_class=0,
                  reco_index=index, truth_index=0, reco_id_hi=11, reco_id_lo=cid,
                  truth_id_hi=12, truth_id_lo=tid)
             for e, index, cid, tid in [(0, 0, 101, 201), (0, 1, 102, 201), (2, 0, 103, 203)]]
    jets = [dict({key: events[e][key] for key in KEY}, truth_jet_pt=pt)
            for e, pt in [(0, 12.), (1, 22.), (0, 21.), (2, -1.)]]
    return dict(events=events, photons=photons, truthPhotons=truths, recoTruthLinks=links, truthJets=jets)


def write_part(path, rows=None, schemas=None, omit=(), auau=False):
    rows = deepcopy(fixture_rows() if rows is None else rows)
    schemas = {**dict(events=EVENT, photons=PHOTON, truthPhotons=TRUTH, recoTruthLinks=LINK, truthJets=JET),
               **(schemas or {})}
    if auau:
        for name in ("events", "photons"):
            schemas[name] = dict(schemas[name], scaled_bit22="I", minimum_bias_pass="I", nominal_event_selection_pass="I")
    file = ROOT.TFile(str(path), "RECREATE")
    for tree, schema in schemas.items():
        if tree not in omit:
            write_tree(file, tree, schema, rows[tree])
    file.Close()


# Unweighted counts; apply the configured sample weight in simulation tests. These expectations
# describe events by hand, independently of the histogram routing implementation.
GRID_EXPECTED = {
    (0, 0): dict(reco=3, truth=2, matched=2, iso=2, joint=2, response=2, unmatched=1),
    (0, 1): dict(reco=4, truth=5, matched=4, iso=4, joint=4, response=4, unmatched=0),
    (1, 0): dict(reco=4, truth=4, matched=3, iso=3, joint=2, response=2, unmatched=0),
    (1, 1): dict(reco=4, truth=4, matched=4, iso=4, joint=4, response=4, unmatched=0),
    (2, 0): dict(reco=2, truth=2, matched=2, iso=2, joint=2, response=2, unmatched=0),
    (2, 1): dict(reco=6, truth=5, matched=5, iso=5, joint=5, response=6, unmatched=0),
}
GRID_SELECTED_EVENTS = [3, 4, 3]


def grid_rows():
    """AuAu boundary, cross-eta, conditional-efficiency and zero-reco fixture.

    In each centrality class two ordinary events contain photons at the strict
    eta exterior edges, both interior sides, the shared edge zero, and NaN.
    Four extra events distinguish cross-cell matching, same-candidate ID/iso,
    multiple response candidates, and a truth denominator with no reco photon.
    """
    rows = {name: [] for name in ("events", "photons", "truthPhotons", "recoTruthLinks", "truthJets")}

    def event(centrality):
        index = len(rows["events"])
        result = dict(ev, source_entry=index, event_id_lo=index + 1,
                      event_sequence=index, physical_event_sequence=index,
                      centrality=centrality, scaled_bit22=1, minimum_bias_pass=1,
                      nominal_event_selection_pass=1)
        rows["events"].append(result)
        return result

    def truth(e, eta, index=0):
        row = dict(tr, **{key: e[key] for key in KEY}, truth_photon_eta=eta,
                   truth_photon_id_hi=12, truth_photon_id_lo=100 * e["event_id_lo"] + index,
                   generator_barcode=42 + index)
        rows["truthPhotons"].append(row)
        return row

    def photon(e, eta, index, target, truth_index=0, tight=True, isolated=True):
        row = dict(ph, **e, photon_eta=eta, candidate_id_hi=11,
                   candidate_id_lo=100 * e["event_id_lo"] + index,
                   photon_encounter_ordinal=100 + index, bdt_input_count=13,
                   bdt_score=.9 if tight else .6, bdt_is_tight=int(tight),
                   bdt_is_nontight=int(not tight), bdt_is_not_tight=int(not tight),
                   bdt_tight_threshold=.8, bdt_nontight_low_threshold=.2,
                   bdt_nontight_high_threshold=.7, iso_r03=0. if isolated else 5.,
                   iso_r04=0. if isolated else 5., iso_r03_threshold=1., iso_r04_threshold=1.,
                   iso_r03_nonisolated_threshold=2., iso_r04_nonisolated_threshold=2.,
                   iso_r03_pass=int(isolated), iso_r04_pass=int(isolated),
                   truth_barcode=target["generator_barcode"])
        rows["photons"].append(row)
        rows["recoTruthLinks"].append(dict({key: e[key] for key in KEY}, reco_type=1, truth_type=1,
            link_class=0, reco_index=index, truth_index=truth_index, reco_id_hi=11,
            reco_id_lo=row["candidate_id_lo"], truth_id_hi=12, truth_id_lo=target["truth_photon_id_lo"]))

    for centrality in (0., 20., 50., 80., -1., float("nan"), 19.999, 49.999, 79.999):
        e = event(centrality)
        for index, eta in enumerate((-.7, -.35, 0., .35, .7, float("nan"))):
            target = truth(e, eta, index)
            photon(e, eta, index, target, index)
    e = event(10.)  # Cross-eta: unmatched in the reco cell, failed truth efficiency in the other cell.
    photon(e, -.35, 0, truth(e, .35))
    e = event(30.)  # Passing ID and isolation on separate candidates does not pass joint efficiency.
    target = truth(e, -.35)
    photon(e, -.35, 0, target, isolated=False)
    photon(e, -.35, 1, target, tight=False)
    e = event(60.)  # One truth efficiency entry, two response fills.
    target = truth(e, .35)
    photon(e, .35, 0, target)
    photon(e, .35, 1, target)
    truth(event(40.), -.35)  # No reco photons; still a denominator and one selected event.
    return rows
