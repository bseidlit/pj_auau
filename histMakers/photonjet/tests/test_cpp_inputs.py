"""Raw PhotonJetTrees reader behavior; no prepared store or analysis output."""
from copy import deepcopy
import json
from pathlib import Path
import struct
import tempfile
import unittest

import ROOT
from cpp_fixtures import EVENT, PHOTON, TRUTH, fixture_rows, write_part

PJ = Path(__file__).resolve().parents[1]
ROOT.gROOT.SetBatch(True)
if not ROOT.gInterpreter.Declare('#include ' + json.dumps(str(PJ / 'histmakers/PhotonJetReader.h'))):
    raise RuntimeError('cannot load direct PhotonJetTrees reader')


class CppInputsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='ppg12-2-cpp-inputs-')
        self.addCleanup(self.temp.cleanup)
        self.source = Path(self.temp.name) / 'source.root'
        write_part(self.source)

    def load(self, simulation=True, require_mb=False, auau=False):
        file = ROOT.TFile.Open(str(self.source))
        self.assertTrue(file and not file.IsZombie())
        try:
            return ROOT.PJ.LoadPart(file, require_mb, simulation, auau)
        finally:
            file.Close()

    def test_retained_raw_values_grouping_and_full_event_identity(self):
        rows = fixture_rows()
        rows['events'][1]['vertex_z'] = -0.0
        write_part(self.source, rows)
        data = self.load()
        self.assertEqual((len(data.events), len(data.photons), len(data.truths)), (4, 3, 3))
        for old, new in zip(rows['events'], data.events):
            for name, kind in EVENT.items():
                expected = old.get(name, 0)
                actual = getattr(new, name)
                if kind == 'D':
                    self.assertEqual(struct.pack('=d', actual), struct.pack('=d', expected), name)
                else:
                    self.assertEqual(actual, expected, name)
        self.assertEqual([list(v) for v in data.reco_rows], [[0, 2], [], [1], []])
        self.assertEqual([list(v) for v in data.truth_rows], [[0], [1], [2], []])
        self.assertEqual([list(v) for v in data.links], [[0, 0], [], [0], []])
        self.assertEqual([p.candidate_id_lo for p in data.photons],
                         [p["candidate_id_lo"] for p in rows["photons"]])
        self.assertEqual([t.truth_photon_id_lo for t in data.truths],
                         [t["truth_photon_id_lo"] for t in rows["truthPhotons"]])
        self.assertEqual([p.original_index for p in data.photons], [0, 0, 1])
        # RNG indices follow within-event raw row order, independent of producer ordinals.
        self.assertEqual([p["photon_encounter_ordinal"] for p in rows["photons"]], [9, 0, 14])
        self.assertEqual([t.prompt_class for t in data.truths], [1, 2, 1])
        self.assertEqual([e.event_id_lo for e in data.events], [1, 2, 3, 1])
        self.assertEqual([e.event_id_hi for e in data.events], [2, 2, 2, 3])
        # A present negative jet stays raw; SelectEvent clamps its stitching pT.
        self.assertEqual(list(data.leading_jet), [21., 22., -1., 0.])
        self.assertEqual(list(data.jet_count), [2, 1, 1, 0])
        self.assertEqual((data.link_rows, data.photon_links, data.jet_rows), (3, 3, 4))

    def test_named_shapes_decode_pp_and_auau_positions(self):
        rows = fixture_rows()
        for photon in rows['photons']:
            for index in range(14):
                photon[f'bdt_input_{index:02d}'] = index + .125
        write_part(self.source, rows)
        names = ('weta', 'wphi', 'e11e33', 'et1', 'et2', 'et3', 'et4', 'e32e35')
        for auau, positions in ((False, (1, 2, 5, 6, 7, 8, 9, 10)), (True, (1, 2, 7, 8, 9, 10, 11, 12))):
            data = self.load(auau=auau)
            self.assertEqual([getattr(data.photons[0], name) for name in names],
                             [position + .125 for position in positions])

    def test_duplicate_or_orphan_identities_and_event_copies_fail(self):
        mutations = {
            'duplicate event': lambda r: r['events'].append(dict(r['events'][0])),
            'negative source index': lambda r: r['events'][0].update(source_file_index=-1),
            'negative source entry': lambda r: r['events'][0].update(source_entry=-1),
            'orphan photon high bits': lambda r: r['photons'][0].update(event_id_hi=99),
            'orphan truth low bits': lambda r: r['truthPhotons'][0].update(event_id_lo=99),
            'orphan jet source': lambda r: r['truthJets'][0].update(source_file_index=99),
            'orphan link': lambda r: r['recoTruthLinks'][0].update(event_id_hi=99),
            'event weight copy': lambda r: r['photons'][0].update(event_weight=123.),
            'event sequence copy': lambda r: r['photons'][0].update(physical_event_sequence=999),
            'duplicate candidate': lambda r: r['photons'].append(dict(r['photons'][0])),
            'duplicate truth': lambda r: r['truthPhotons'].append(dict(r['truthPhotons'][0])),
        }
        for name, mutate in mutations.items():
            with self.subTest(case=name):
                rows = fixture_rows()
                mutate(rows)
                write_part(self.source, rows)
                with self.assertRaises(Exception):
                    self.load()

    def test_invalid_duplicate_and_mismatched_links_fail(self):
        for name, value in (('reco_index', -1), ('truth_index', 6), ('reco_index', 5),
                            ('reco_id_hi', 99), ('reco_id_lo', 999), ('truth_id_hi', 99), ('truth_id_lo', 999)):
            with self.subTest(field=name):
                rows = fixture_rows()
                rows['recoTruthLinks'][0][name] = value
                write_part(self.source, rows)
                with self.assertRaises(Exception):
                    self.load()
        rows = fixture_rows()
        rows['recoTruthLinks'].append(dict(rows['recoTruthLinks'][0]))
        write_part(self.source, rows)
        with self.assertRaisesRegex(Exception, 'duplicate photon link'):
            self.load()

    def test_missing_associations_are_explicit(self):
        write_part(self.source, omit=('recoTruthLinks', 'truthJets'))
        data = self.load()
        self.assertFalse(data.links_available)
        self.assertFalse(data.jets_available)
        self.assertEqual([list(v) for v in data.links], [[-1, -1], [], [-1], []])
        self.assertEqual(list(data.jet_count), [0, 0, 0, 0])

    def test_missing_and_wrong_type_branches_fail(self):
        cases = [
            {'events': {key: value for key, value in EVENT.items() if key != 'event_id_hi'}},
            {'events': dict(EVENT, event_id_hi='D')},
            {'events': dict(EVENT, scaled_bit22='D')},
            {'photons': dict(PHOTON, bdt_input_13='I')},
            {'photons': dict(PHOTON, native_et1='I')},
            {'truthPhotons': dict(TRUTH, prompt_class='D')},
            # Diagnostics omitted from stored records still require exact raw types.
            {'photons': {key: value for key, value in PHOTON.items() if key != 'photon_phi'}},
            {'photons': dict(PHOTON, photon_phi='I')},
            {'photons': dict(PHOTON, photon_encounter_ordinal='D')},
            {'photons': dict(PHOTON, bdt_is_not_tight='D')},
            {'photons': dict(PHOTON, bdt_input_count='D')},
            {'photons': dict(PHOTON, truth_matched='D')},
            {'truthPhotons': dict(TRUTH, truth_photon_phi='I')},
            {'truthPhotons': dict(TRUTH, source_role='D')},
        ]
        for schemas in cases:
            with self.subTest(schemas=schemas):
                rows = fixture_rows()
                # Keep values convertible by the malformed scalar fixture writer.
                for row in rows['photons']:
                    row['bdt_input_13'] = 0
                for row in rows['truthPhotons']:
                    row['truth_photon_phi'] = 0
                write_part(self.source, rows, schemas=schemas)
                with self.assertRaises(Exception):
                    self.load()
        write_part(self.source, omit=('truthPhotons',))
        with self.assertRaisesRegex(Exception, 'truthPhotons'):
            self.load()
        write_part(self.source)
        with self.assertRaisesRegex(Exception, 'scaled_bit22'):
            self.load(require_mb=True, auau=True)

    def test_empty_part_and_optional_data_truth_are_valid(self):
        write_part(self.source, {name: [] for name in fixture_rows()})
        data = self.load()
        self.assertEqual((len(data.events), len(data.photons), len(data.truths)), (0, 0, 0))
        rows = fixture_rows()
        rows['recoTruthLinks'] = []
        write_part(self.source, rows, omit=('truthPhotons', 'truthJets', 'recoTruthLinks'))
        data = self.load(simulation=False)
        self.assertEqual(len(data.events), 4)
        self.assertEqual(len(data.truths), 0)


if __name__ == '__main__':
    unittest.main()
