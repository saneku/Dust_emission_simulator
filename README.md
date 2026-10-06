# Dust_emission_simulator

Mimicking WRF-Chem's dust emission schemes outside of WRF-Chem itself: `dust_opt=1` (GOCART) and `dust_opt=3` (AFWA), both coupled with the GOCART aerosol module. Useful for offline tuning/testing of dust-emission parameters against WRF output without rerunning the full model.

## Contents

- `gocart_source_dust.py` — GOCART dust emission scheme (`gocart_source_dust()`).
- `afwa_source_dust.py` — AFWA dust emission scheme (`afwa_source_dust()`), including a selectable soil-moisture correction (gravimetric, volumetric, or GOCART-simple).
- `gocart_python.py` / `gocart_python.ipynb` — driver script/notebook that reads a WRF output file and runs `gocart_source_dust`.
- `afwa_python.py` / `afwa_python.ipynb` — driver script/notebook that reads a WRF output file and runs `afwa_source_dust`.
- `gocart_plt_orgnl_wrfoutput.py` / `afwa_plt_orgnl_wrfoutput.py` — plot the dust flux already present in the WRF output file (`DUST_EMIS` et al.) for comparison against the recomputed flux.
- `utils.py` — shared grid/projection/plotting setup (domain size, cartopy projection, colormaps). Edit this file to match your own domain.
- `data/` — expected location of input WRF netCDF files (`grid.nc`, `gocart.nc`, `afwa.nc`).

## Requirements

Maps are rendered with [cartopy](https://scitools.org.uk/cartopy/). Install with conda (recommended — cartopy's GEOS/PROJ dependencies are easiest to get right this way):

```
conda env create -f environment.yml
conda activate dust-emission-simulator
```

or with pip:

```
pip install -r requirements.txt
```

## Usage

1. Edit `utils.py` to match your WRF domain: grid size (`nx`, `ny`), projection parameters (`cen_lat`, `cen_lon`, `true_lat1`, `true_lat2`, `dx`, `dy`), and `wrf_dir` (defaults to `./data/`).
2. Place your WRF output netCDF file(s) in `data/` (e.g. `afwa.nc`, `gocart.nc`), containing at minimum: `UST`, `SMOIS`, `ISLTYP`, `SNOWH`, `ZNT`, `ALT`, `EROD`, `CLAYFRAC`, `SANDFRAC`, `XLAND`, `Times`.
3. Run the notebook or script for the scheme you want:
   ```
   python afwa_python.py
   # or
   python gocart_python.py
   ```
   Both can also be run as Jupyter notebooks (`afwa_python.ipynb`, `gocart_python.ipynb`), which additionally expose widgets to interactively adjust the tuning parameters below.
4. Each run saves a PNG of the recomputed instantaneous dust flux (e.g. `afwa_inst_flux.png`).

### Tuning parameters

Both schemes accept keyword tuning parameters passed through `**tuning_params`:

- `afwa_source_dust`: `alpha` (global tuning constant), `gamma` (erodibility exponent), `smtune` (soil-moisture scaling), `ustune` (friction-velocity scaling).
- `gocart_source_dust`: `C_factor` (global emission scaling constant).

For AFWA, the soil-moisture correction scheme is selected via the `smois_opt` constant near the top of `afwa_source_dust()`: `0` = gravimetric SM, `1` = volumetric SM, `2` = GOCART-simple (current default).

## How to cite

If you use this code, please cite:

Ukhov, A., Ahmadov, R., Grell, G., and Stenchikov, G.: Improving dust simulations in WRF-Chem v4.1.3 coupled with the GOCART aerosol module, Geosci. Model Dev., 14, 473–493, https://doi.org/10.5194/gmd-14-473-2021, 2021.
