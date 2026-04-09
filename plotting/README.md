# Plotting

## PlotEvtCharacter.C

Compares EvtCharacter histograms for **data_aa**, **photon10_aa**, and **jet10_aa**.

**Input files (defined in the macro):**
- `/sphenix/u/bseidlitz/work/pj_auau/histMakers/results/evt_characterdata_aa.root`
- `/sphenix/u/bseidlitz/work/pj_auau/histMakers/results/evt_characterphoton10_aa.root`
- `/sphenix/u/bseidlitz/work/pj_auau/histMakers/results/evt_characterjet10_aa.root`

**Run from ROOT:**
```bash
root -l -b -q '/sphenix/u/bseidlitz/work/pj_auau/plotting/PlotEvtCharacter.C'
```
Or from inside ROOT:
```cpp
.x /sphenix/u/bseidlitz/work/pj_auau/plotting/PlotEvtCharacter.C
```

**Output:** PDF and PNG written to `/sphenix/u/bseidlitz/work/pj_auau/plotting/figs/`:
- evt_character_centrality
- evt_character_vertexz
- evt_character_centrality_cluster10
- evt_character_cluster_pt_above10 (log-y)

Data = points with errors; photon10_aa = blue line; jet10_aa = red line.
