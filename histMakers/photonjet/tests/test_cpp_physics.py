"""Direct-loop physics, exact keyed RNG, and centrality × eta cell behavior."""
from array import array
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import ROOT
from ruamel.yaml import YAML
from cpp_fixtures import fixture_rows, write_part, grid_rows, GRID_EXPECTED, GRID_SELECTED_EVENTS
from cpp_validation import compare, root_snapshot, histogram_snapshot
from fixtures import correction_file

PJ = Path(__file__).resolve().parents[1]
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('libyaml-cpp', '', True)
ROOT.gSystem.Load('libRooUnfold')
if not ROOT.gInterpreter.Declare(
    '#include ' + json.dumps(str(PJ / 'histmakers/PhotonJetHistMaker.C')) + '\n' +
    '#include ' + json.dumps(str(PJ / 'histmakers/InjectMbdEff.C')) + '\n' +
    '#include ' + json.dumps(str(PJ / 'plotting/PhotonJetPlotInput.h'))):
    raise RuntimeError('cannot load direct photonjet maker')
if not ROOT.gInterpreter.Declare(r'''
PJ::Cuts PJCppTestCuts(const std::string &path) { return PJ::LoadCuts(YAML::LoadFile(path)); }
std::vector<double> PJCppGaussianMoments(unsigned count) {
    double sum=0, sum2=0, sum3=0, sum4=0, cross=0, previous=0;
    unsigned tails=0;
    for (unsigned i=0; i<count; ++i) {
        const auto seed=PJ::CandidateSeed(42, "rng-fixture/simulation/pp_photonjet/part/1", 7, 1ULL<<43, i/4, i%4);
        const double x=PJ::Gaussian(seed);
        sum+=x; sum2+=x*x; sum3+=x*x*x; sum4+=x*x*x*x;
        if(i) cross+=x*previous;
        previous=x; tails+=std::abs(x)>3;
    }
    return {sum/count,sum2/count,sum3/count,sum4/count,cross/(count-1),double(tails)};
}
'''):
    raise RuntimeError('cannot load direct physics fixture helpers')


def seed_python(seed, part, source, hi, lo, candidate):
    """Independent byte encoding; do not rename the established RNG byte prefix."""
    mask = (1 << 64) - 1
    data = bytearray(b'PJRDF-RNG-v1')
    for value in (seed, len(part.encode())):
        data.extend(value.to_bytes(8, 'little'))
    data.extend(part.encode())
    for value in (source, hi, lo, candidate):
        data.extend(value.to_bytes(8, 'little'))
    result = 14695981039346656037
    for value in data:
        result = ((result ^ value) * 1099511628211) & mask
    result = (result + 0x9e3779b97f4a7c15) & mask
    result = ((result ^ (result >> 30)) * 0xbf58476d1ce4e5b9) & mask
    result = ((result ^ (result >> 27)) * 0x94d049bb133111eb) & mask
    result ^= result >> 31
    return result or 0x9e3779b97f4a7c15


def integral(hist):
    return hist.Integral(0, hist.GetNbinsX() + 1)


class CppPhysicsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='ppg12-2-cpp-physics-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.yaml = YAML()
        self.cfg = self.yaml.load((PJ / 'configs/pp/config_pj_pp_nom.yaml').read_text())
        self.cfg['analysis'].update(pT_bins=[16., 22., 36.], pT_bins_truth=[8., 22., 45.], apply_trigger_eff_correction=0)
        for kind in ('data', 'mc'):
            for index in range(3):
                self.cfg['analysis'][f'cluster_eres_{kind}_p{index}'] = 0.
        self.cfg['analysis']['unfold']['reweight'] = 0
        self.cfg['photonjet'].update(external_mbd_eff_file='', release='fixture-v1', weight_mode='stored')
        self.cfg['output']['var_type'] = 'PPG12_2_cpp_fixture'
        for key, stem in (('data_outfile', 'data'), ('eff_outfile', 'eff'), ('response_outfile', 'response'), ('final_outfile', 'final')):
            self.cfg['output'][key] = str(self.base / stem)
        self.rows = fixture_rows()
        self.config = self.base / 'config.yaml'
        self.listing = self.base / 'parts.list'
        self.source = self.base / 'source.root'

    def save(self):
        with self.config.open('w') as stream:
            self.yaml.dump(self.cfg, stream)
        return self.config

    def input(self, product='sim_signal', rows=None, omit=(), auau=False):
        raw = deepcopy(self.rows if rows is None else rows)
        if product == 'data':
            for name in ('truthPhotons', 'recoTruthLinks', 'truthJets'):
                raw[name] = []
        write_part(self.source, raw, omit=omit, auau=auau)
        master = self.base / 'files.txt'
        master.write_text(str(self.source) + '\n')
        self.listing.write_text('0 ' + str(self.source) + '\n')
        key = {'data': 'data_files_list', 'sim_signal': 'sim_signal_files_list', 'sim_inclusive': 'sim_inclusive_files_list'}[product]
        self.cfg['photonjet'][key] = str(master)

    def output(self, stem, tag=''):
        path = self.base / f'{stem}_PPG12_2_cpp_fixture.root'
        return path if not tag else path.parent / 'chunks' / f'{path.stem}.{tag}.root'

    def run_maker(self, product='sim_signal', tag='', audit=False):
        self.save()
        audit_path = str(self.base / ('audit-' + (tag or 'main') + '.txt')) if audit else ''
        ROOT.PJ.Run(str(self.config), product, str(self.listing), tag, audit_path)
        stem = 'data' if product == 'data' else 'eff_jet' if product == 'sim_inclusive' else 'eff'
        return self.output(stem, tag), self.output('response', tag)

    def test_all_shipped_configs_and_invalid_edges(self):
        for path in sorted((PJ / 'configs').glob('*/*.yaml')):
            with self.subTest(config=path.name):
                cuts = ROOT.PJCppTestCuts(str(path))
                self.assertEqual(cuts.random_seed, 42)
                bins = ROOT.PJ.BinLayout(cuts)
                self.assertEqual((bins.nCentrality(), bins.nEta()), (1, 1))
        for section, key, value in (
            ('analysis', 'eta_bins', []), ('analysis', 'eta_bins', [-.7]),
            ('analysis', 'eta_bins', [-.7, 0, 0, .7]), ('analysis', 'eta_bins', [-.7, float('nan'), .7]),
            ('analysis', 'eta_bins', [.7, -.7]), ('analysis', 'pT_bins', [16, 16, 36]),
            ('photonjet', 'centrality_bins', [0, 20]), ('photonjet', 'weight_mode', 'stord'),
            ('photonjet', 'random_seed', 0), ('photonjet', 'match_source', 'links_typo'),
            ('photonjet', 'system', 'AuAu'), ('photonjet', 'threshold_source', 'trees'),
            ('photonjet', 'flag_check', 'ignore'), ('analysis', 'vertex_cut', 'bad'),
            ('analysis', 'vertex_cut', float('nan')), ('analysis', 'use_topo_iso', 0),
        ):
            with self.subTest(key=key, value=value):
                original = deepcopy(self.cfg)
                self.cfg[section][key] = value
                with self.assertRaises(Exception):
                    ROOT.PJCppTestCuts(str(self.save()))
                self.cfg = original
        self.cfg['photonjet'].update(system='auau', threshold_source='tree')
        for edges in ([0], [-1, 20], [0, 101], [0, 20, 20], [0, float('inf')]):
            self.cfg['photonjet']['centrality_bins'] = edges
            with self.subTest(centrality=edges), self.assertRaises(Exception):
                ROOT.PJCppTestCuts(str(self.save()))

    def test_bin_boundaries_and_pp_inclusive_slot(self):
        pp = ROOT.PJ.BinLayout(ROOT.PJCppTestCuts(str(self.save())))
        for value in (-100., 0., 100., float('nan')):
            self.assertEqual(pp.CentralityBin(value), 0)
        self.cfg['photonjet'].update(system='auau', threshold_source='tree', centrality_bins=[0., 20., 50., 80.])
        self.cfg['analysis']['eta_bins'] = [-.7, 0., .7]
        bins = ROOT.PJ.BinLayout(ROOT.PJCppTestCuts(str(self.save())))
        for value, expected in ((-1., -1), (0., 0), (19.999, 0), (20., 1), (49.999, 1), (50., 2), (79.999, 2), (80., -1), (float('nan'), -1)):
            self.assertEqual(bins.CentralityBin(value), expected)
        for value, expected in ((-.7, -1), (-.6999, 0), (-.001, 0), (0., 1), (.6999, 1), (.7, -1), (float('nan'), -1)):
            self.assertEqual(bins.EtaBin(value), expected)
        self.assertEqual(str(bins.cell(2, 1).Name('h_all_cluster')), 'h_all_cluster_cent2_eta1')
        for cell in ((-1, 0), (3, 0), (0, 2)):
            with self.assertRaises(Exception):
                bins.cell(*cell)

    def test_sample_boundaries_sentinels_vertex_and_candidate_weights(self):
        self.input()
        cuts = ROOT.PJCppTestCuts(str(self.save()))
        bins = ROOT.PJ.BinLayout(cuts)
        sample = ROOT.PPG12.GetSampleConfig('photon10')
        vertex = ROOT.PJ.VertexTable()
        with ROOT.TFile.Open(str(self.source)) as file:
            data = ROOT.PJ.LoadPart(file, False, True, False)
        event = data.events[0]
        select = lambda: ROOT.PJ.SelectEvent(cuts, bins, True, sample, vertex, data, 0)
        self.assertEqual(select().weight, 2.)
        event.event_weight = 1.
        self.assertFalse(getattr(select(), 'pass'))
        cuts.weight_mode = 'sample_map'
        event.vertex_z = 0.
        self.assertFalse(getattr(select(), 'pass'))
        event.vertex_z = 1.
        self.assertTrue(getattr(select(), 'pass'))
        self.assertEqual(select().weight, sample.weight)
        truth = data.truths[0]
        for pt, accepted in ((13.999, False), (14., True), (21.999, True), (22., False)):
            truth.truth_photon_pt = pt
            self.assertEqual(getattr(select(), 'pass'), accepted, pt)
        truth.truth_photon_pt = 20.
        for z, accepted in ((30., True), (30.001, False), (-30., True), (-30.001, False)):
            event.vertex_z = z
            self.assertEqual(getattr(select(), 'pass'), accepted, z)
        event.vertex_z = 1.
        event.event_weight = 2.
        cuts.weight_mode = 'recompute'
        self.assertEqual(select().weight, sample.weight)
        cuts.weight_mode = 'sample_map'
        sample = ROOT.PPG12.GetSampleConfig('jet12')
        for leading, count, accepted in ((13.999, 1, False), (14., 1, True), (20.999, 1, True), (21., 1, False), (18., 0, False)):
            data.leading_jet[0] = leading
            data.jet_count[0] = count
            self.assertEqual(getattr(select(), 'pass'), accepted, leading)
        cuts.apply_shower_shape_windows = 0
        photon = data.photons[0]
        for et, accepted in ((23., True), (23.001, False)):
            photon.photon_et = et
            self.assertEqual(ROOT.PJ.SelectPhoton(photon, cuts, bins, True, sample, 3., 0).accepted, accepted)
        cuts.weight_mode = 'stored'
        self.assertEqual(ROOT.PJ.SelectPhoton(photon, cuts, bins, True, sample, 3., 0).weight, 3.)

    def test_zero_reco_split_id_iso_and_nonfiducial_match(self):
        self.input()
        self.cfg['photonjet']['match_source'] = 'links'
        main, response = self.run_maker()
        with ROOT.TFile.Open(str(main)) as file:
            reco, iso, ideff, all_eff = [file.Get('eff_' + name + '_eta_0') for name in ('reco', 'iso', 'id', 'all')]
            self.assertEqual(integral(reco.GetTotalHistogram()), 4.)
            self.assertEqual(integral(reco.GetPassedHistogram()), 2.)
            self.assertEqual(integral(iso.GetTotalHistogram()), 2.)
            self.assertEqual(integral(iso.GetPassedHistogram()), 2.)
            self.assertEqual(integral(ideff.GetPassedHistogram()), 0.)
            self.assertEqual(integral(all_eff.GetPassedHistogram()), 0.)
            self.assertEqual(file.Get('h_nontight_iso_cluster_notmatch_0').GetBinContent(3), 2.)
            self.assertEqual(file.Get('h_all_cluster_0').GetBinContent(0), 2.)
            self.assertEqual(file.Get('h_pj_cutflow').GetBinContent(1), 4.)
            for efficiency in (reco, iso, ideff, all_eff, file.Get('eff_converts_eta_0')):
                self.assertTrue(efficiency.UsesWeights())
                self.assertEqual(efficiency.GetStatisticOption(), ROOT.TEfficiency.kBUniform)
                self.assertEqual(efficiency.GetWeight(), 1.)
        with ROOT.TFile.Open(str(response)) as file:
            self.assertEqual(file.Get('h_response_full_0').GetEntries(), 0.)
            self.assertEqual(integral(file.Get('response_matrix_full_0').Hfakes()), 0.)
        self.assertEqual(len(root_snapshot(main)['objects']), 41)
        self.assertEqual(len(root_snapshot(response)['objects']), 5)
        # Plot input validation reads this same completed legacy-layout output.
        with ROOT.PJPlot.Open(str(main)) as plotted:
            self.assertTrue(plotted.Get('h_all_cluster_0'))
            self.assertTrue(plotted.Get('photonjet_provenance'))

    def test_multiple_response_candidates_and_response_only_smearing(self):
        for photon in self.rows['photons'][:1] + self.rows['photons'][2:]:
            photon.update(photon_et=20., bdt_score=.9, iso_r04=0.)
        self.input()
        self.cfg['photonjet']['match_source'] = 'links'
        main, response = self.run_maker(tag='plain', audit=True)
        baseline = root_snapshot(main)['objects']
        with ROOT.TFile.Open(str(response)) as file:
            self.assertEqual(file.Get('h_response_full_0').GetEntries(), 2.)
            self.assertEqual(integral(file.Get('h_pT_truth_response_0')), 4.)
        self.cfg['analysis']['cluster_eres_data_p2'] = .1
        smeared, _ = self.run_maker(tag='smear', audit=True)
        self.assertEqual(compare(baseline, root_snapshot(smeared)['objects']), [])
        cuts = ROOT.PJ.Cuts()
        cuts.eres_data_p2 = .1
        for gaussian, expected in ((-20., 0.), (0., 40.), (2., 44.), (20., 80.)):
            self.assertEqual(ROOT.PJ.ResponseET(cuts, 40., 20., gaussian), expected)
        # The audit uses original within-event positions 0 and 1, not ordinals 9 and 14.
        audit = (self.base / 'audit-plain.txt').read_text().splitlines()
        self.assertEqual(sorted(row.split()[0].rsplit('/', 1)[1] for row in audit), ['0', '1'])
        for row in audit:
            identity, seed, gaussian, reco, truth, weight = row.split()
            part, source, hi, lo, candidate = identity.rsplit('/', 4)
            self.assertEqual(int(seed), seed_python(42, part, int(source), int(hi), int(lo), int(candidate)))
            self.assertEqual(float.fromhex(gaussian).hex(), float(ROOT.PJ.Gaussian(int(seed))).hex())
            self.assertEqual((float.fromhex(reco), float.fromhex(truth), float.fromhex(weight)), (20., 20., 2.))

    def test_missing_links_bad_associations_and_nonfinite_weights_preserve_generation(self):
        self.input(omit=('recoTruthLinks',))
        main, response = self.run_maker()
        paths = (main, response, Path(str(main) + '.complete.yaml'))
        before = {p: p.read_bytes() for p in paths}
        self.cfg['photonjet']['match_source'] = 'links'
        with self.assertRaisesRegex(Exception, 'unavailable|missing'):
            self.run_maker()
        self.assertEqual({p: p.read_bytes() for p in paths}, before)
        # A malformed link is rejected even when its entire event fails vertex selection.
        rows = fixture_rows()
        for row in rows['events'] + rows['photons']:
            row['vertex_z'] = 45.
        rows['recoTruthLinks'][0]['truth_index'] = 99
        self.input(rows=rows)
        with self.assertRaises(Exception):
            self.run_maker()
        self.assertEqual({p: p.read_bytes() for p in paths}, before)
        rows = fixture_rows()
        rows['events'][1]['event_weight'] = float('nan')  # This is a zero-reco event.
        self.input(rows=rows)
        with self.assertRaisesRegex(Exception, 'nonfinite'):
            self.run_maker()
        self.assertEqual({p: p.read_bytes() for p in paths}, before)
        self.source.write_bytes(b'not a ROOT file')
        with self.assertRaises(Exception):
            self.run_maker()
        self.assertEqual({p: p.read_bytes() for p in paths}, before)

    def test_original_lists_and_invalid_product_preserve_outputs(self):
        self.input()
        main, response = self.run_maker()
        paths = (main, response, Path(str(main) + '.complete.yaml'))
        before = {path: path.read_bytes() for path in paths}
        master = self.base / 'files.txt'
        master_text = master.read_text()
        listing_text = self.listing.read_text()
        duplicate = self.base / 'copied-original.root'
        duplicate.write_bytes(self.source.read_bytes())  # Distinct path with the same ROOT UUID.
        missing = self.base / 'missing-original.root'
        cases = (
            ('duplicate index', master_text, f'0 {self.source}\n0 {self.source}\n'),
            ('negative index', master_text, f'-1 {self.source}\n'),
            ('out-of-range index', master_text, f'1 {self.source}\n'),
            ('noninteger index', master_text, f'oops {self.source}\n'),
            ('missing source', master_text + str(missing) + '\n', f'0 {self.source}\n1 {missing}\n'),
            ('duplicate ROOT UUID', master_text + str(duplicate) + '\n', f'0 {self.source}\n1 {duplicate}\n'),
        )
        try:
            for name, original_list, chosen_list in cases:
                with self.subTest(case=name):
                    master.write_text(original_list)
                    self.listing.write_text(chosen_list)
                    with self.assertRaises(Exception):
                        self.run_maker()
                    self.assertEqual({path: path.read_bytes() for path in paths}, before)
                    self.assertEqual(list(self.base.rglob('*.tmp.*')), [])
        finally:
            master.write_text(master_text)
            self.listing.write_text(listing_text)
        with self.assertRaisesRegex(Exception, 'unknown product'):
            ROOT.PJ.Run(str(self.config), 'sim_singal', str(self.listing), '')
        self.assertEqual({path: path.read_bytes() for path in paths}, before)

    def test_output_aliases_cannot_replace_original_config_or_input_list(self):
        self.input()
        main, response = self.run_maker()
        paths = (main, response, Path(str(main) + '.complete.yaml'))
        before = {path: path.read_bytes() for path in paths}
        original_cfg = deepcopy(self.cfg)
        for index, (protected, hardlink) in enumerate(((self.source, False), (self.config, False),
                                                      (self.listing, False), (self.source, True))):
            with self.subTest(protected=protected.name, hardlink=hardlink):
                self.cfg = deepcopy(original_cfg)
                self.cfg['output']['var_type'] = f'alias{index}'
                self.cfg['output']['eff_outfile'] = str(self.base / 'protected')
                self.save()
                alias = self.base / f'protected_alias{index}.root'
                if hardlink:
                    alias.hardlink_to(protected)
                else:
                    alias.symlink_to(protected)
                original = protected.read_bytes()
                try:
                    with self.assertRaisesRegex(Exception, 'collid'):
                        ROOT.PJ.Run(str(self.config), 'sim_signal', str(self.listing), '')
                    self.assertEqual(protected.read_bytes(), original)
                    self.assertEqual({path: path.read_bytes() for path in paths}, before)
                finally:
                    alias.unlink()
        self.cfg = original_cfg
        self.save()

    def test_mbd_import_does_not_create_missing_target(self):
        self.input()
        main, response = self.run_maker()
        paths = (main, response, Path(str(main) + '.complete.yaml'))
        before = {path: path.read_bytes() for path in paths}
        correction = self.base / 'correction.root'
        correction_file(correction, edges=self.cfg['analysis']['pT_bins_truth'], corrected=True)
        correction_before = correction.read_bytes()
        self.cfg['photonjet']['external_mbd_eff_file'] = str(correction)
        target = self.base / 'absent-target.root'
        with self.assertRaises(Exception):
            ROOT.PJ.InjectCorrection(str(self.save()), str(target))
        self.assertFalse(target.exists())
        self.assertEqual(correction.read_bytes(), correction_before)
        self.assertEqual({path: path.read_bytes() for path in paths}, before)
        self.assertEqual(list(self.base.rglob('*.tmp.*')), [])

    def test_response_finalizer_empty_weighted_and_flow_contract(self):
        cases = ([], [(18., 20., 2.), (24., 30., .4), (18.5, 20., 1.2)],
                 [(12., 7., .25), (40., 46., 1.25), (18., 20., 2.)])
        for index, points in enumerate(cases):
            with self.subTest(points=points):
                edges, truth_edges = array('d', [16., 22., 36.]), array('d', [8., 22., 45.])
                measured = ROOT.TH1D(f'fixed_reco_{index}', '', 2, edges)
                truth = ROOT.TH1D(f'fixed_truth_{index}', '', 2, truth_edges)
                matrix = ROOT.TH2D(f'fixed_matrix_{index}', '', 2, edges, 2, truth_edges)
                for histogram in (measured, truth, matrix):
                    histogram.Sumw2()
                    histogram.SetDirectory(0)
                reference = ROOT.RooUnfoldResponse(measured, truth, matrix, f'fixed_reference_{index}', '', False)
                for reco_pt, truth_pt, weight in points:
                    reference.Fill(reco_pt, truth_pt, weight)
                    measured.Fill(reco_pt, weight)
                    truth.Fill(truth_pt, weight)
                    matrix.Fill(reco_pt, truth_pt, weight)
                finalized = ROOT.PJ.FinalizeResponse(measured, truth, matrix)
                for method in ('Hresponse', 'Hmeasured', 'Htruth', 'Hfakes'):
                    expected = histogram_snapshot(getattr(reference, method)())
                    actual = histogram_snapshot(getattr(finalized, method)())
                    self.assertEqual(compare(expected, actual), [], method)
                self.assertFalse(finalized.HasFakes())

    def test_flags_floor_run_list_and_empty_data_efficiency(self):
        for photon in self.rows['photons']:
            if photon['bdt_score'] > .8:
                photon['bdt_is_nontight'] = 0
        self.input('data')
        main, _ = self.run_maker('data')
        before = main.read_bytes()
        self.cfg['analysis']['tight']['bdt_min_intercept'] = .1
        with self.assertRaisesRegex(Exception, 'stored flags'):
            self.run_maker('data')
        self.assertEqual(main.read_bytes(), before)
        self.cfg['photonjet'].update(flag_check='off', min_photon_et=16.)
        with self.assertRaisesRegex(Exception, 'minimum photon ET'):
            self.run_maker('data')
        self.assertEqual(main.read_bytes(), before)
        self.cfg['photonjet']['min_photon_et'] = 15.
        runs = self.base / 'runs.txt'
        runs.write_text('49999\n')
        self.cfg['analysis']['run_list_file'] = str(runs)
        rejected, _ = self.run_maker('data', tag='rejected')
        with ROOT.TFile.Open(str(rejected)) as file:
            self.assertEqual(file.Get('h_pj_cutflow').GetBinContent(1), 4.)
            self.assertEqual(file.Get('h_pj_cutflow').GetBinContent(2), 0.)
            for name in ('reco', 'iso', 'id', 'all', 'converts'):
                eff = file.Get('eff_' + name + '_eta_0')
                self.assertTrue(eff.UsesWeights())
                self.assertEqual(eff.GetStatisticOption(), ROOT.TEfficiency.kBUniform)
                self.assertEqual(integral(eff.GetTotalHistogram()), 0.)

    def test_common_window_failure_still_checks_stored_flags(self):
        rows = fixture_rows()
        rows['photons'] = rows['photons'][:1]
        rows['photons'][0].update(bdt_input_05=.99, bdt_is_tight=1, iso_r04_pass=0)
        self.input('data', rows=rows)
        self.cfg['photonjet'].update(apply_shower_shape_windows=1, flag_check='strict')
        with self.assertRaisesRegex(Exception, 'stored flags'):
            self.run_maker('data')
        self.assertFalse(self.output('data').exists())
        self.cfg['photonjet']['flag_check'] = 'report'
        main, _ = self.run_maker('data')
        with ROOT.TFile.Open(str(main)) as file:
            flags = file.Get('h_pj_flagcheck')
            self.assertEqual([flags.GetBinContent(i) for i in range(1, 5)], [1., 1., 0., 1.])
            self.assertEqual(integral(file.Get('h_common_cluster_0')), 0.)

    def test_native_rng_fixed_vectors_and_distribution(self):
        vectors = json.loads((PJ / 'tests/rng_vectors.json').read_text())
        self.assertEqual(ROOT.gROOT.GetVersion(), vectors['ROOT_version'])
        for vector in vectors['vectors']:
            self.assertEqual(float(ROOT.PJ.Gaussian(int(vector['seed']))).hex(), vector['gaussian'])
        mean, second, third, fourth, correlation, tails = ROOT.PJCppGaussianMoments(100000)
        self.assertLess(abs(mean), .016)
        self.assertLess(abs(second - 1), .025)
        self.assertLess(abs(third), .075)
        self.assertLess(abs(fourth - 3), .17)
        self.assertLess(abs(correlation), .016)
        self.assertGreater(tails, 185)
        self.assertLess(tails, 355)
        print('RNG moments:', mean, second, third, fourth, correlation, tails, flush=True)

    def test_original_high_bits_and_independent_keyed_seed_vectors(self):
        part = 'fixture-v1/simulation/pp_photonjet/part/17'
        seeds = []
        for values in ((42, part, 0, 2, 1, 0), (42, part, 0, 2 + (1 << 40), 1, 0),
                       (42, part, 0, 2, 1 + (1 << 40), 0), (42, part, 0, 2, 1, 1),
                       (42, part + '0', 0, 2, 1, 0), (43, part, 0, 2, 1, 0),
                       (42, part, 1, 2, 1, 0)):
            expected = seed_python(*values)
            self.assertEqual(ROOT.PJ.CandidateSeed(*values), expected)
            self.assertGreater(expected, 1 << 32)
            seeds.append(expected)
        self.assertEqual(len(set(seeds)), len(seeds))
        self.assertNotEqual(ROOT.PJ.Gaussian(1), ROOT.PJ.Gaussian((1 << 32) + 1))
        with self.assertRaises(Exception):
            ROOT.PJ.Gaussian(0)

    def grid_config(self):
        self.cfg['analysis'].update(eta_bins=[-.7, 0., .7], mc_iso_scale=1., mc_iso_shift=0.)
        self.cfg['photonjet'].update(system='auau', centrality_bins=[0., 20., 50., 80.],
            threshold_source='tree', match_source='links', apply_shower_shape_windows=0,
            weight_mode='sample_map', vertex_weight_file='')
        sample_map = self.base / 'samples.yaml'
        sample_map.write_text(json.dumps({'ranges': [{'first': 0, 'last': 0, 'sample': 'photon10'}]}))
        self.cfg['photonjet']['sample_map_signal'] = str(sample_map)
        self.grid_weight = float(ROOT.PPG12.GetSampleConfig('photon10').weight)
        self.input(rows=grid_rows(), auau=True)

    def test_three_centrality_two_eta_physics_and_global_accounting(self):
        self.grid_config()
        main, response = self.run_maker(audit=True)
        snapshots = [root_snapshot(path, titles=True) for path in (main, response)]
        self.assertEqual([len(s['objects']) for s in snapshots], [206, 25])
        provenance = self.yaml.load(snapshots[0]['metadata']['photonjet_provenance'])
        layout = provenance['layout']
        self.assertEqual(layout['centrality_edges'], [0., 20., 50., 80.])
        self.assertEqual(layout['eta_edges'], [-.7, 0., .7])
        self.assertEqual(layout['eta_policy'], 'outer-exclusive-interior-to-higher')
        self.assertEqual(provenance['selected_events_by_centrality'], GRID_SELECTED_EVENTS)
        with ROOT.TFile.Open(str(main)) as file:
            self.assertEqual(file.GetListOfKeys().GetEntries(), 209)
            self.assertEqual(file.Get('h_pj_cutflow').GetBinContent(1), 13.)
            self.assertEqual(file.Get('h_pj_cutflow').GetBinContent(2), 10.)
            self.assertEqual(file.Get('h_pj_vertex_z').GetEntries(), 10.)
            self.assertEqual(file.Get('h_pj_input_parts').GetBinContent(1), 1.)
            for (icent, ieta), counts in GRID_EXPECTED.items():
                suffix = f'_cent{icent}_eta{ieta}'
                self.assertEqual(integral(file.Get('h_all_cluster' + suffix)), counts['reco'] * self.grid_weight)
                self.assertEqual(integral(file.Get('h_truth_pT' + suffix)), counts['truth'] * self.grid_weight)
                self.assertEqual(integral(file.Get('h_tight_iso_cluster_notmatch' + suffix)), counts['unmatched'] * self.grid_weight)
                for name, denominator, numerator in (
                    ('reco', 'truth', 'matched'), ('iso', 'matched', 'iso'),
                    ('id', 'iso', 'joint'), ('all', 'truth', 'joint')):
                    efficiency = file.Get('eff_' + name + suffix)
                    self.assertEqual(integral(efficiency.GetTotalHistogram()), counts[denominator] * self.grid_weight, (icent, ieta, name))
                    self.assertEqual(integral(efficiency.GetPassedHistogram()), counts[numerator] * self.grid_weight, (icent, ieta, name))
                    self.assertTrue(efficiency.UsesWeights())
                    self.assertEqual(efficiency.GetStatisticOption(), ROOT.TEfficiency.kBUniform)
                names = [name for name in snapshots[0]['objects'] if name.endswith(suffix)]
                self.assertEqual(len(names), 33)
                for name in names:
                    self.assertIn('centrality', snapshots[0]['objects'][name]['title'])
                    self.assertIn('eta', snapshots[0]['objects'][name]['title'])
        with ROOT.TFile.Open(str(response)) as file:
            self.assertEqual(file.GetListOfKeys().GetEntries(), 28)
            for (icent, ieta), counts in GRID_EXPECTED.items():
                suffix = f'_cent{icent}_eta{ieta}'
                self.assertEqual(file.Get('h_response_full' + suffix).GetEntries(), counts['response'])
                self.assertEqual(integral(file.Get('h_pT_truth_response' + suffix)), counts['response'] * self.grid_weight)
                self.assertEqual(integral(file.Get('response_matrix_full' + suffix).Hfakes()), 0.)
        self.assertEqual(len((self.base / 'audit-main.txt').read_text().splitlines()), 20)
        # Current plotters must reject the recorded layout before choosing any cell.
        original = main.read_bytes()
        with self.assertRaisesRegex(Exception, 'plotting.*multiple'):
            ROOT.PJPlot.Open(str(main))
        self.assertEqual(main.read_bytes(), original)

    def test_each_centrality_row_matches_independent_run_and_preserves_rng(self):
        self.grid_config()
        full = self.run_maker(tag='all', audit=True)
        expected = [root_snapshot(path)['objects'] for path in full]
        audits = []
        edges = [0., 20., 50., 80.]
        for icent in range(3):
            self.cfg['photonjet']['centrality_bins'] = edges[icent:icent + 2]
            actual = [root_snapshot(path)['objects'] for path in self.run_maker(tag=f'cent{icent}', audit=True)]
            for original, single in zip(expected, actual):
                for ieta in range(2):
                    old_suffix, new_suffix = f'_cent{icent}_eta{ieta}', f'_cent0_eta{ieta}'
                    left = {name[:-len(old_suffix)]: value for name, value in original.items() if name.endswith(old_suffix)}
                    right = {name[:-len(new_suffix)]: value for name, value in single.items() if name.endswith(new_suffix)}
                    self.assertEqual(compare(left, right), [])
            audits.extend((self.base / f'audit-cent{icent}.txt').read_text().splitlines())
        self.assertEqual(sorted(audits), sorted((self.base / 'audit-all.txt').read_text().splitlines()))
        # Inserting another class does not alter original candidate seeds/draws.
        self.cfg['photonjet']['centrality_bins'] = [0., 10., 20., 50., 80., 90., 100.]
        main, response = self.run_maker(tag='extra', audit=True)
        enlarged = set((self.base / 'audit-extra.txt').read_text().splitlines())
        self.assertTrue(set(audits).issubset(enlarged))
        with ROOT.TFile.Open(str(main)) as file:
            for ieta in range(2):
                self.assertEqual(integral(file.Get(f'h_truth_pT_cent5_eta{ieta}')), 0.)
                self.assertEqual(integral(file.Get(f'eff_reco_cent5_eta{ieta}').GetTotalHistogram()), 0.)
        with ROOT.TFile.Open(str(response)) as file:
            for ieta in range(2):
                self.assertEqual(file.Get(f'h_response_full_cent5_eta{ieta}').GetEntries(), 0.)


if __name__ == '__main__':
    unittest.main()
