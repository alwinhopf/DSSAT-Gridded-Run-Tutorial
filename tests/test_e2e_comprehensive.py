#!/usr/bin/env python3
"""
Comprehensive provider integration test (Python) — tests weather & soil sources.
- Tests each weather source: Open-Meteo, Daymet, NASA-POWER, GridMET, AgERA5, NASA-POWER CHIRPS
- Tests each soil source: SoilGrids, SoilGrids Online, SSURGO, HWSD
- Uses isolated pytest fixtures to run each provider independently
- Calls live network providers; transient transport outages are retried and reported as skips
- Gracefully skips sources with missing optional deps or API keys

Run: pytest tests/test_e2e_comprehensive.py -v
"""

import io
import os
import sqlite3
import sys
import tempfile
import time
from contextlib import redirect_stdout
from pathlib import Path
import pandas as pd
import pytest

# Use the installed, CI-pinned package by default. Workspace developers can
# explicitly opt into a sibling checkout while changing dssatutils.
_HERE = Path(__file__).parent
_REPO = _HERE.parent
_WORKSPACE = _REPO.parent
if os.getenv("DSSAT_USE_WORKSPACE_SIBLINGS", "") == "1":
    sys.path.insert(0, str(_WORKSPACE / "dssatutils" / "python"))
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / "python_scripts"))

pytestmark = pytest.mark.integration

try:
    import dssatutils
except ImportError:
    pytest.skip("dssatutils not installed", allow_module_level=True)

# Check optional dependencies
HAS_SF = False
HAS_GEOPANDAS = False
try:
    import geopandas as gpd
    from shapely.geometry import Point
    HAS_GEOPANDAS = True
except ImportError:
    pass

try:
    import sf
    HAS_SF = True
except ImportError:
    pass

# Test fixtures
GLOBAL_POINTS = pd.DataFrame({
    "ID": ["GLOBAL_EU", "GLOBAL_AS", "GLOBAL_AF"],
    "LAT": [52.0, 30.9, -0.5],
    "LONG": [5.0, 75.8, 37.0],
})

US_POINTS = pd.DataFrame({
    "ID": ["US_CORN", "US_WHEAT"],
    "LAT": [40.0, 38.5],
    "LONG": [-97.0, -99.5],
})

def make_gdf(df):
    """Convert DataFrame to GeoDataFrame if geopandas is available."""
    if HAS_GEOPANDAS:
        return gpd.GeoDataFrame(
            df,
            geometry=[Point(xy) for xy in zip(df["LONG"], df["LAT"])],
            crs="EPSG:4326",
        )
    return df

@pytest.fixture
def work_dir():
    """Create a temporary work directory."""
    with tempfile.TemporaryDirectory(prefix="dssat_e2e_") as d:
        yield d


def _read_log_tail(path, lines=30):
    path = Path(path)
    if not path.exists():
        return ""
    try:
        return "\n".join(path.read_text(encoding="utf-8", errors="replace").splitlines()[-lines:])
    except OSError:
        return ""


def _is_transient_provider_failure(diagnostic):
    """Return True only for provider/transport failures worth retrying."""
    text = str(diagnostic).strip().lower()
    markers = (
        "could not resolve host",
        "temporary failure in name resolution",
        "connection reset",
        "connection refused",
        "connection timed out",
        "operation timed out",
        "read timed out",
        "remote disconnected",
        "service unavailable",
        "bad gateway",
        "gateway timeout",
        "too many requests",
        "http 429",
        "http 500",
        "http 502",
        "http 503",
        "http 504",
        "429 client error",
        "500 server error",
        "502 server error",
        "503 server error",
        "504 server error",
        "status code 429",
        "status code 500",
        "status code 502",
        "status code 503",
        "status code 504",
        "ssl connect error",
        "downloaded file failed netcdf validation",
        "invalid/corrupt netcdf cache file",
    )
    return any(marker in text for marker in markers)


@pytest.mark.parametrize(
    "diagnostic",
    (
        "500 Server Error: INTERNAL SERVER ERROR for url: https://daymet.example/api",
        "HTTP 503 while requesting weather data",
        "status code 429 returned by provider",
    ),
)
def test_transient_provider_error_wording(diagnostic):
    """Recognize common HTTP-client messages so live outages do not fail CI."""
    assert _is_transient_provider_failure(diagnostic)


class _Tee:
    """Mirror provider output to pytest while retaining it for diagnostics."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def _run_weather_provider(name, call, expected_paths, log_file):
    """Run a live provider without allowing missing outputs to pass silently."""
    attempts = max(1, int(os.getenv("DSSAT_PROVIDER_RETRIES", "3")))
    errors = []
    provider_output = []
    missing = list(expected_paths)

    for attempt in range(1, attempts + 1):
        captured = io.StringIO()
        try:
            with redirect_stdout(_Tee(sys.stdout, captured)):
                call()
        except ImportError:
            raise
        except Exception as exc:  # provider diagnostics are evaluated below
            errors.append(f"{type(exc).__name__}: {exc}")
        finally:
            output = captured.getvalue().strip()
            if output:
                provider_output.append(f"Attempt {attempt} output:\n{output}")

        missing = [Path(path) for path in expected_paths if not Path(path).exists()]
        if not missing:
            for path in expected_paths:
                _check_wth_file(Path(path))
            return

        diagnostic = "\n".join(errors + provider_output + [_read_log_tail(log_file)])
        if not _is_transient_provider_failure(diagnostic):
            break
        if attempt < attempts:
            time.sleep(min(5 * attempt, 15))

    diagnostic = "\n".join(errors + provider_output + [_read_log_tail(log_file)]).strip()
    missing_names = ", ".join(path.name for path in missing)
    if diagnostic and _is_transient_provider_failure(diagnostic):
        pytest.skip(
            f"{name} remained unavailable after {attempts} attempts; "
            f"missing: {missing_names}. Diagnostics: {diagnostic}"
        )
    pytest.fail(
        f"{name} did not produce every expected weather file; missing: "
        f"{missing_names}. Diagnostics: {diagnostic or 'none'}"
    )

# --- WEATHER SOURCE TESTS ---

class TestWeatherSources:
    """Test all weather data sources."""

    def test_weather_openmeteo_global(self, work_dir):
        """Test Open-Meteo (global, free, no key)."""
        wth_dir = Path(work_dir) / "weather_openmeteo"
        wth_dir.mkdir(parents=True)
        
        gdf = make_gdf(GLOBAL_POINTS)
        log_file = Path(work_dir) / "openmeteo.log"
        
        expected = [wth_dir / f"{pid}.WTH" for pid in GLOBAL_POINTS["ID"]]
        _run_weather_provider(
            "Open-Meteo",
            lambda: dssatutils.process_weather_openmeteo(
                shapefile=gdf,
                start_year=2010,
                end_year=2010,
                output_dir=str(wth_dir),
                id_col="ID",
                lat_col="LAT",
                lon_col="LONG",
                n_cores=1,
                log_file=str(log_file),
            ),
            expected,
            log_file,
        )

    @pytest.mark.skipif(not hasattr(dssatutils, "process_weather_daymet"), 
                        reason="process_weather_daymet not available")
    def test_weather_daymet_us(self, work_dir):
        """Test Daymet (US only, free, no key)."""
        wth_dir = Path(work_dir) / "weather_daymet"
        wth_dir.mkdir(parents=True)
        
        gdf = make_gdf(US_POINTS)
        log_file = Path(work_dir) / "daymet.log"
        
        try:
            expected = [wth_dir / f"{pid}.WTH" for pid in US_POINTS["ID"]]
            _run_weather_provider(
                "Daymet",
                lambda: dssatutils.process_weather_daymet(
                    shapefile=gdf,
                    start_year=2010,
                    end_year=2010,
                    output_dir=str(wth_dir),
                    id_col="ID",
                    lat_col="LAT",
                    lon_col="LONG",
                    n_cores=1,
                    log_file=str(log_file),
                ),
                expected,
                log_file,
            )
        except ImportError as e:
            pytest.skip(f"Daymet optional dep missing: {e}")

    @pytest.mark.skipif(not hasattr(dssatutils, "process_weather_nasapower"),
                        reason="process_weather_nasapower not available")
    def test_weather_nasapower_global(self, work_dir):
        """Test NASA-POWER (global, free, no key)."""
        wth_dir = Path(work_dir) / "weather_nasapower"
        wth_dir.mkdir(parents=True)
        
        gdf = make_gdf(GLOBAL_POINTS)
        log_file = Path(work_dir) / "nasapower.log"
        
        expected = [wth_dir / f"{pid}.WTH" for pid in GLOBAL_POINTS["ID"]]
        _run_weather_provider(
            "NASA POWER",
            lambda: dssatutils.process_weather_nasapower(
                shapefile=gdf,
                start_year=2010,
                end_year=2010,
                output_dir=str(wth_dir),
                id_col="ID",
                lat_col="LAT",
                lon_col="LONG",
                n_cores=1,
                log_file=str(log_file),
            ),
            expected,
            log_file,
        )

    @pytest.mark.skipif(not hasattr(dssatutils, "process_weather_gridmet"),
                        reason="process_weather_gridmet not available")
    def test_weather_gridmet_us(self, work_dir):
        """Test GridMET (US only, free, no key)."""
        wth_dir = Path(work_dir) / "weather_gridmet"
        wth_dir.mkdir(parents=True)
        
        gdf = make_gdf(US_POINTS)
        log_file = Path(work_dir) / "gridmet.log"
        cache_dir = Path(work_dir) / "gridmet_cache"
        cache_dir.mkdir(parents=True)
        
        try:
            expected = [wth_dir / f"{pid}.WTH" for pid in US_POINTS["ID"]]
            _run_weather_provider(
                "GridMET",
                lambda: dssatutils.process_weather_gridmet(
                    shapefile=gdf,
                    start_year=2010,
                    end_year=2010,
                    output_dir=str(wth_dir),
                    id_col="ID",
                    lat_col="LAT",
                    lon_col="LONG",
                    n_cores=1,
                    log_file=str(log_file),
                    gridmet_cache_dir=str(cache_dir),
                ),
                expected,
                log_file,
            )
        except ImportError as e:
            pytest.skip(f"GridMET optional dep missing: {e}")

    @pytest.mark.skipif(not hasattr(dssatutils, "process_weather_agera5"),
                        reason="process_weather_agera5 not available")
    def test_weather_agera5_global(self, work_dir):
        """Test AgERA5 (global, requires CDS key)."""
        if not os.getenv("CDSAPI_RC"):
            pytest.skip("CDS API key not configured (CDSAPI_RC env var)")
        
        wth_dir = Path(work_dir) / "weather_agera5"
        wth_dir.mkdir(parents=True)
        
        gdf = make_gdf(GLOBAL_POINTS)
        log_file = Path(work_dir) / "agera5.log"
        cache_dir = Path(work_dir) / "agera5_cache"
        cache_dir.mkdir(parents=True)
        
        try:
            expected = [wth_dir / f"{pid}.WTH" for pid in GLOBAL_POINTS["ID"]]
            _run_weather_provider(
                "AgERA5",
                lambda: dssatutils.process_weather_agera5(
                    shapefile=gdf,
                    start_year=2010,
                    end_year=2010,
                    output_dir=str(wth_dir),
                    id_col="ID",
                    lat_col="LAT",
                    lon_col="LONG",
                    n_cores=1,
                    log_file=str(log_file),
                    agera5_cache_dir=str(cache_dir),
                ),
                expected,
                log_file,
            )
        except ImportError as e:
            pytest.skip(f"AgERA5 optional dep missing: {e}")

    @pytest.mark.skipif(not hasattr(dssatutils, "process_weather_nasapower_chirps"),
                        reason="process_weather_nasapower_chirps not available")
    def test_weather_nasapower_chirps_global(self, work_dir):
        """Test NASA-POWER CHIRPS (global, combines NASA-POWER with CHIRPS rainfall)."""
        wth_dir = Path(work_dir) / "weather_nasapower_chirps"
        wth_dir.mkdir(parents=True)
        
        gdf = make_gdf(GLOBAL_POINTS)
        log_file = Path(work_dir) / "chirps.log"
        cache_dir = Path(work_dir) / "chirps_cache"
        cache_dir.mkdir(parents=True)
        
        try:
            expected = [wth_dir / f"{pid}.WTH" for pid in GLOBAL_POINTS["ID"]]
            _run_weather_provider(
                "NASA POWER + CHIRPS",
                lambda: dssatutils.process_weather_nasapower_chirps(
                    shapefile=gdf,
                    start_year=2010,
                    end_year=2010,
                    output_dir=str(wth_dir),
                    id_col="ID",
                    lat_col="LAT",
                    lon_col="LONG",
                    n_cores=1,
                    log_file=str(log_file),
                    chirps_cache_dir=str(cache_dir),
                ),
                expected,
                log_file,
            )
        except ImportError as e:
            pytest.skip(f"NASA-POWER CHIRPS optional dep missing: {e}")

# --- SOIL SOURCE TESTS ---

class TestSoilSources:
    """Test all soil data sources."""

    @pytest.mark.skipif(not hasattr(dssatutils, "process_soils_soilgrids_online"),
                        reason="process_soils_soilgrids_online not available")
    def test_soil_soilgrids_online_global(self, work_dir):
        """Test SoilGrids Online (global, free, REST).
        
        WARNING: REST API is slow (~10-30s per point). This test may timeout.
        """
        sol_dir = Path(work_dir) / "soil_soilgrids_online"
        sol_dir.mkdir(parents=True)
        
        gdf = make_gdf(GLOBAL_POINTS)
        csv_path = sol_dir / "soil_map.csv"
        
        try:
            dssatutils.process_soils_soilgrids_online(
                gridfile=gdf,
                soilfile_csv_path=str(csv_path),
                output_sol_dir=str(sol_dir),
                id_col="ID",
            )
            
            sols = list(sol_dir.glob("*.SOL"))
            if len(sols) == 0:
                pytest.skip("No .SOL files written (SoilGrids REST coverage gap or timeout)")
            
            for sol_file in sols:
                _check_sol_file(sol_file)
        except RuntimeError as e:
            if "No soil data extracted" in str(e):
                pytest.skip(f"SoilGrids provider returned no data: {e}")
            raise

    @pytest.mark.skipif(not hasattr(dssatutils, "process_soils_soilgrids"),
                        reason="process_soils_soilgrids not available")
    def test_soil_soilgrids_global(self, work_dir):
        """Test SoilGrids (global, offline VRT or online)."""
        sol_dir = Path(work_dir) / "soil_soilgrids"
        sol_dir.mkdir(parents=True)
        
        gdf = make_gdf(GLOBAL_POINTS)
        csv_path = sol_dir / "soil_map.csv"
        
        master = sol_dir / "dummy_master.SOL"
        master.write_text(
            "*SOILS: Test Soils Master\n"
            "*TEST0001     TEST_SOIL       52.000    5.000\n"
            "@SITE        COUNTRY          LAT     LONG SCS FAMILY\n"
            " TEST_SOIL   World         52.000    5.000\n"
            "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE\n"
            "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001\n"
            "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC\n"
            "     5   -99 0.100 0.200 0.300  1.00  10.0  1.40  1.00  20.0  40.0   0.0   -99   -99   -99   -99   -99\n"
            "    15   -99 0.100 0.200 0.300  0.80  10.0  1.40  0.80  20.0  40.0   0.0   -99   -99   -99   -99   -99\n",
            encoding="utf-8",
        )
        dssatutils.process_soils_soilgrids(
            grid_points=gdf,
            source_sol_file=str(master),
            output_csv_path=str(csv_path),
            output_sol_dir=str(sol_dir),
            id_col="ID",
            numeric_only_ids=False,
        )
        sols = list(sol_dir.glob("*.SOL"))
        assert sols, "offline SoilGrids fixture produced no .SOL files"
        for sol_file in sols:
            _check_sol_file(sol_file)

    @pytest.mark.skipif(not hasattr(dssatutils, "process_soils_ssurgo"),
                        reason="process_soils_ssurgo not available")
    def test_soil_ssurgo_us(self, work_dir):
        """Test SSURGO (US only)."""
        sol_dir = Path(work_dir) / "soil_ssurgo"
        sol_dir.mkdir(parents=True)
        
        gdf = make_gdf(US_POINTS)
        csv_path = sol_dir / "soil_map.csv"
        
        failure_csv = sol_dir / "soil_map_download_failures.csv"
        diagnostics = []
        last_error = None
        sols = []

        for attempt in range(3):
            try:
                dssatutils.process_soils_ssurgo(
                    grid_points=gdf,
                    output_dir_csv=str(csv_path),
                    output_dir_individual=str(sol_dir),
                    n_cores=1,
                    id_col="ID",
                    lat_col="LAT",
                    long_col="LONG",
                )
                last_error = None
            except Exception as exc:
                last_error = exc

            sols = list(sol_dir.glob("*.SOL"))
            if sols:
                break

            diagnostics = _read_ssurgo_failure_reasons(failure_csv)
            if last_error is not None:
                diagnostics.append(str(last_error))
            provider_unavailable = bool(diagnostics) and all(
                _is_transient_ssurgo_failure(reason) for reason in diagnostics
            )
            if not provider_unavailable or attempt == 2:
                break
            time.sleep(5 * (attempt + 1))

        if not sols and diagnostics and all(
            _is_transient_ssurgo_failure(reason) for reason in diagnostics
        ):
            pytest.skip(
                "SSURGO SDA remained unavailable after 3 attempts: "
                + " | ".join(diagnostics)
            )
        assert sols, (
            "SSURGO produced no .SOL files. Diagnostics: "
            + (" | ".join(diagnostics) if diagnostics else "none")
        )
        for sol_file in sols:
            _check_sol_file(sol_file)

    @pytest.mark.skipif(not hasattr(dssatutils, "process_soils_hwsd"),
                        reason="process_soils_hwsd not available")
    def test_soil_hwsd_global(self, work_dir):
        """Test HWSD (global, SQLite-based)."""
        sol_dir = Path(work_dir) / "soil_hwsd"
        sol_dir.mkdir(parents=True)
        
        gdf = make_gdf(GLOBAL_POINTS)
        csv_path = sol_dir / "soil_map.csv"
        
        rasterio = pytest.importorskip("rasterio")
        import numpy as np
        from rasterio.transform import from_origin

        raster = sol_dir / "dummy_hwsd.tif"
        with rasterio.open(
            raster, "w", driver="GTiff", width=10, height=10, count=1,
            dtype="int16", crs="EPSG:4326", transform=from_origin(-180, 90, 36, 18),
        ) as dst:
            dst.write(np.ones((10, 10), dtype="int16"), 1)
        database = sol_dir / "dummy_hwsd.sqlite"
        with sqlite3.connect(database) as con:
            con.execute(
                "CREATE TABLE HWSD2_LAYERS (HWSD2_SMU_ID INTEGER, SEQUENCE INTEGER, "
                "SHARE REAL, TOPDEP REAL, BOTDEP REAL, SAND REAL, CLAY REAL, SILT REAL, "
                "BULK REAL, ORG_CARBON REAL, COARSE REAL)"
            )
            con.execute(
                "INSERT INTO HWSD2_LAYERS VALUES (1,1,100,0,200,40,20,40,1.4,1.0,0)"
            )
        dssatutils.process_soils_hwsd(
            grid_points=gdf,
            hwsd_raster_file=str(raster),
            hwsd_db_file=str(database),
            output_csv_path=str(csv_path),
            output_sol_dir=str(sol_dir),
            id_col="ID", lat_col="LAT", long_col="LONG",
        )
        sols = list(sol_dir.glob("*.SOL"))
        assert sols, "synthetic HWSD fixture produced no .SOL files"
        for sol_file in sols:
            _check_sol_file(sol_file)

# --- VALIDATION HELPERS ---

def _read_ssurgo_failure_reasons(path):
    """Return per-point SSURGO failure reasons written by dssatutils."""
    if not path.exists():
        return []
    try:
        failures = pd.read_csv(path)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return []
    if "reason" not in failures.columns:
        return []
    return failures["reason"].dropna().astype(str).tolist()


def _is_transient_ssurgo_failure(reason):
    """Identify SDA transport failures without hiding data/schema defects."""
    text = str(reason).strip().lower()
    if text.startswith((
        "network: sda spatial query failed",
        "network: sda horizon query failed",
    )):
        return True
    transport_markers = (
        "could not resolve host",
        "connection reset",
        "connection refused",
        "connection timed out",
        "operation timed out",
        "read timed out",
        "service unavailable",
        "too many requests",
        "http 429",
        "http 500",
        "http 502",
        "http 503",
        "http 504",
    )
    return any(marker in text for marker in transport_markers)

def _check_wth_file(path):
    """Validate a .WTH file."""
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    
    assert len(lines) >= 5, f"{path}: too short"
    assert lines[0].startswith("$WEATHER"), f"{path}: missing $WEATHER header"
    
    has_date_header = any("@  DATE" in ln for ln in lines[:5])
    assert has_date_header, f"{path}: missing DATE column header"
    
    data_lines = [ln for ln in lines if ln.strip()[0].isdigit()]
    assert len(data_lines) >= 100, f"{path}: too few daily rows ({len(data_lines)})"
    
    all_data = "\n".join(data_lines)
    assert "nan" not in all_data.lower(), f"{path}: NaN found in data"

def _check_sol_file(path):
    """Validate a .SOL file."""
    with open(path) as f:
        txt = f.read().lower()
    
    assert "*soils" in txt or txt.lstrip().startswith("*"), f"{path}: missing *SOILS header"
    assert "slb" in txt or "layer" in txt, f"{path}: missing layer table"
    assert "nan" not in txt, f"{path}: NaN found in soil file"

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
