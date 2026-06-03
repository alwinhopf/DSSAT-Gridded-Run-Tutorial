#!/usr/bin/env python3
"""
run_bioenergy_analysis.py
=========================
Master script to run DSSAT gridded simulations across the Ames, Iowa example field
using the adapted RYAS8201.WHX template (Cereal Rye proxy) and perform
spatial, environmental, and economic analysis.
"""

import os
import sys
import shutil
import subprocess
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path

# Set up python path to import dssatutils if needed, though it's installed in the env
sys.path.insert(0, "/Users/alwinhopf/Documents/GitHub/dssatutils/python")
from dssatutils.weather_daymet import process_weather_daymet

def main():
    # 1. Setup paths
    base_dir = Path(__file__).resolve().parent
    conda_env_python = "/opt/homebrew/Caskroom/miniconda/base/envs/dssat_grid/bin/python3"
    dssat_exe = "/Users/alwinhopf/Documents/GitHub/DSSAT48/dscsm048"
    
    # 2. Generate the grid points first
    print("Generating field grid points...")
    # Centroid of the field in Ames, Iowa
    center_lat = 42.0200
    center_lon = -93.6200
    
    UTM_CRS = "EPSG:26915"
    WGS84_CRS = "EPSG:4326"
    
    # Define a 500m x 1000m field at 50m spacing
    half_width = 250.0
    half_height = 500.0
    spacing = 50.0
    
    from pyproj import Transformer
    from shapely.geometry import Point
    
    transformer_to_utm = Transformer.from_crs(WGS84_CRS, UTM_CRS, always_xy=True)
    center_x, center_y = transformer_to_utm.transform(center_lon, center_lat)
    
    xs = np.arange(center_x - half_width, center_x + half_width + 0.1, spacing)
    ys = np.arange(center_y - half_height, center_y + half_height + 0.1, spacing)
    
    points = []
    for y in ys:
        for x in xs:
            points.append(Point(x, y))
            
    gdf = gpd.GeoDataFrame(geometry=points, crs=UTM_CRS).to_crs(WGS84_CRS)
    gdf["LAT"] = gdf.geometry.y.round(6)
    gdf["LONG"] = gdf.geometry.x.round(6)
    gdf["ID"] = [f"{i+1:08d}" for i in range(len(gdf))]
    
    gridpoints_dir = base_dir / "gridpoints"
    gridpoints_dir.mkdir(exist_ok=True)
    shp_path = gridpoints_dir / "example_field_50m.shp"
    gdf.to_file(shp_path)
    print(f"Generated {len(gdf)} grid points at: {shp_path}")
    
    # 3. Download centroid weather and replicate to avoid Daymet rate limits
    weather_dir = base_dir / "weather" / "bioenergy_field_variability_50m_DAYMET"
    weather_dir.mkdir(parents=True, exist_ok=True)
    
    # Clear old weather files to ensure fresh run
    for f in weather_dir.glob("*.WTH"):
        f.unlink()
        
    print("Downloading weather for the field centroid...")
    centroid_gdf = gdf.iloc[len(gdf)//2 : len(gdf)//2 + 1].copy() # Middle point as centroid
    centroid_gdf["ID"] = "00000001"
    
    log_file = weather_dir / "download_errors.log"
    process_weather_daymet(
        shapefile=centroid_gdf,
        start_year=1982,
        end_year=1983,
        output_dir=str(weather_dir),
        id_col="ID",
        lat_col="LAT",
        lon_col="LONG",
        n_cores=1,
        log_file=str(log_file)
    )
    
    centroid_wth = weather_dir / "00000001.WTH"
    if not centroid_wth.exists():
        print("CRITICAL: Centroid weather download failed.")
        sys.exit(1)
        
    print("Replicating centroid weather for all field points...")
    for _, row in gdf.iterrows():
        pid = str(row["ID"])
        if pid == "00000001":
            continue
        shutil.copy2(centroid_wth, weather_dir / f"{pid}.WTH")
        
    # 4. Write custom config.yml for our run
    config_dir = base_dir / ".experiment_configs"
    config_dir.mkdir(exist_ok=True)
    config_file = config_dir / "bioenergy_config.yml"
    
    config_content = f"""# Unified configuration for bioenergy cover crop simulation
project_name:         "bioenergy_field_variability"
grid_spacing_meters:  50
crop_extension:       "WH"
run_tag:              ""
run_name_style:       "grid"
run_name_override:    "bioenergy_field_variability"

use_existing_point_shapefile: true
existing_point_shapefile_path: "{shp_path}"

weather_source:       "DAYMET"
weather_start_year:   1982
weather_end_year:     1983

soil_source:          "SSURGO"

template_file_name:   "RYAS8201.WHX"
run_mode:             "experiment"
treatment_start:      1
treatment_end:        2

bundle_genotype_files: true
"""
    
    with open(config_file, "w") as fh:
        fh.write(config_content)
    print(f"Wrote configuration to: {config_file}")
    
    # 5. Run dssat_main_pipeline.py
    env = os.environ.copy()
    env["DSSAT_CONFIG_FILE"] = str(config_file)
    env["DSSAT_EXE"] = dssat_exe
    
    pipeline_script = base_dir / "dssat_main_pipeline.py"
    
    print("=" * 60)
    print("RUNNING GRIDDED DSSAT PIPELINE FOR CEREAL RYE COVER CROP...")
    print("=" * 60)
    
    cmd = [conda_env_python, str(pipeline_script)]
    result = subprocess.run(cmd, env=env, check=False)
    
    if result.returncode != 0:
        print(f"CRITICAL: Pipeline run failed with exit code: {result.returncode}")
        sys.exit(result.returncode)
        
    print("=" * 60)
    print("PIPELINE RUN COMPLETED. STARTING POST-PROCESSING & ANALYSIS...")
    print("=" * 60)
    
    # 6. Load Results
    results_path = base_dir / "results" / "bioenergy_field_variability_results.csv"
    if not results_path.exists():
        print(f"CRITICAL: Merged results file not found at: {results_path}")
        sys.exit(1)
        
    df = pd.read_csv(results_path)
    print(f"Loaded {len(df)} simulation result rows from: {results_path}")
    
    # Clean up duplicate points if any (caused by multiprocessing artifacts)
    df = df.drop_duplicates(subset=["point_id", "treatment", "year_planting"])
    
    df_early = df[df["treatment"] == 1].copy()
    df_late = df[df["treatment"] == 2].copy()
    
    if df_early.empty or df_late.empty:
        print("CRITICAL: One or both treatments have no data. Verify DSSAT simulation outputs.")
        sys.exit(1)
        
    comparison = df_early.merge(df_late, on=["point_id", "latitude", "longitude"], suffixes=("_early", "_late"))
    
    comparison["biomass_diff"] = comparison["top_weight_kg_ha_early"] - comparison["top_weight_kg_ha_late"]
    comparison["n_leaching_diff"] = comparison["nitrate_leaching_kg_ha_late"] - comparison["nitrate_leaching_kg_ha_early"]
    
    # Economic assumptions:
    # Biofuel feedstock price: $55 per metric ton of dry matter -> $0.055 per kg
    price_per_kg = 0.055
    cost_per_ha = 111.0
    
    comparison["net_profit_early"] = comparison["top_weight_kg_ha_early"] * price_per_kg - cost_per_ha
    comparison["net_profit_late"] = comparison["top_weight_kg_ha_late"] * price_per_kg - cost_per_ha
    comparison["profit_diff"] = comparison["net_profit_early"] - comparison["net_profit_late"]
    
    # Delineate Suitability Zones based on early planting biomass and N leaching:
    # Zone A (High Suitability): Biomass > 2000 kg/ha and N Leaching < 15 kg N/ha
    # Zone B (Moderate Suitability): Biomass 1000 - 2000 kg/ha
    # Zone C (Low Suitability/Unprofitable): Biomass < 1000 kg/ha (cannot cover costs)
    def assign_zone(row):
        biomass = row["top_weight_kg_ha_early"]
        leaching = row["nitrate_leaching_kg_ha_early"]
        if pd.isna(biomass):
            return 3
        if biomass > 2000 and (pd.isna(leaching) or leaching < 15):
            return 1  # High Suitability
        elif biomass >= 1000:
            return 2  # Moderate Suitability
        else:
            return 3  # Low Suitability / Unprofitable
            
    comparison["zone"] = comparison.apply(assign_zone, axis=1)
    
    # Print summary table
    print("\n" + "=" * 50)
    print("SUMMARY STATISTICS ACROSS THE SUBFIELD")
    print("=" * 50)
    print(f"{'Metric':<30} | {'Early Planting':<15} | {'Late Planting':<15}")
    print("-" * 66)
    print(f"{'Mean Biomass (kg/ha)':<30} | {comparison['top_weight_kg_ha_early'].mean():.1f} | {comparison['top_weight_kg_ha_late'].mean():.1f}")
    print(f"{'Max Biomass (kg/ha)':<30} | {comparison['top_weight_kg_ha_early'].max():.1f} | {comparison['top_weight_kg_ha_late'].max():.1f}")
    print(f"{'Mean N Leaching (kg N/ha)':<30} | {comparison['nitrate_leaching_kg_ha_early'].mean():.2f} | {comparison['nitrate_leaching_kg_ha_late'].mean():.2f}")
    print(f"{'Mean SOC Delta (kg C/ha)':<30} | {comparison['soil_organic_carbon_delta_kg_C_ha_early'].mean():.2f} | {comparison['soil_organic_carbon_delta_kg_C_ha_late'].mean():.2f}")
    print(f"{'Mean Net Profit ($/ha)':<30} | ${comparison['net_profit_early'].mean():.2f} | ${comparison['net_profit_late'].mean():.2f}")
    print(f"{'Profitable Area Fraction':<30} | {(comparison['net_profit_early'] > 0).mean() * 100:.1f}% | {(comparison['net_profit_late'] > 0).mean() * 100:.1f}%")
    print("=" * 50 + "\n")
    
    # 7. Plot Spatial Maps
    print("Generating spatial maps...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    vmin_bio = min(comparison["top_weight_kg_ha_early"].min(), comparison["top_weight_kg_ha_late"].min())
    vmax_bio = max(comparison["top_weight_kg_ha_early"].max(), comparison["top_weight_kg_ha_late"].max())
    
    vmin_n = min(comparison["nitrate_leaching_kg_ha_early"].min(), comparison["nitrate_leaching_kg_ha_late"].min())
    vmax_n = max(comparison["nitrate_leaching_kg_ha_early"].max(), comparison["nitrate_leaching_kg_ha_late"].max())
    
    # Plot 1: Biomass Early
    sc1 = axes[0, 0].scatter(
        comparison["longitude"], comparison["latitude"],
        c=comparison["top_weight_kg_ha_early"], cmap="YlGn",
        s=40, vmin=vmin_bio, vmax=vmax_bio, edgecolors="grey", linewidths=0.2
    )
    plt.colorbar(sc1, ax=axes[0, 0], label="Aboveground Biomass (kg/ha)")
    axes[0, 0].set_title("Biomass at Termination: Early Planting (Sept 15)")
    axes[0, 0].set_aspect("equal")
    
    # Plot 2: Biomass Late
    sc2 = axes[0, 1].scatter(
        comparison["longitude"], comparison["latitude"],
        c=comparison["top_weight_kg_ha_late"], cmap="YlGn",
        s=40, vmin=vmin_bio, vmax=vmax_bio, edgecolors="grey", linewidths=0.2
    )
    plt.colorbar(sc2, ax=axes[0, 1], label="Aboveground Biomass (kg/ha)")
    axes[0, 1].set_title("Biomass at Termination: Late Planting (Oct 15)")
    axes[0, 1].set_aspect("equal")
    
    # Plot 3: N Leaching Early
    sc3 = axes[1, 0].scatter(
        comparison["longitude"], comparison["latitude"],
        c=comparison["nitrate_leaching_kg_ha_early"], cmap="OrRd",
        s=40, vmin=vmin_n, vmax=vmax_n, edgecolors="grey", linewidths=0.2
    )
    plt.colorbar(sc3, ax=axes[1, 0], label="Cumulative Nitrate Leached (kg N/ha)")
    axes[1, 0].set_title("N Leaching: Early Planting (Sept 15)")
    axes[1, 0].set_aspect("equal")
    
    # Plot 4: Delineated Suitability Zones
    from matplotlib.colors import ListedColormap
    colors = ["#2ca02c", "#bcbd22", "#d62728"]  # Green, Yellow-Green, Red
    cmap_zone = ListedColormap(colors)
    
    sc4 = axes[1, 1].scatter(
        comparison["longitude"], comparison["latitude"],
        c=comparison["zone"], cmap=cmap_zone,
        s=40, edgecolors="grey", linewidths=0.2
    )
    cbar4 = plt.colorbar(sc4, ax=axes[1, 1], ticks=[1.33, 2.0, 2.66])
    cbar4.ax.set_yticklabels(["High", "Moderate", "Unprofitable"])
    cbar4.set_label("Suitability / Feasibility Zone")
    axes[1, 1].set_title("Subfield Cover Crop Management Zones")
    axes[1, 1].set_aspect("equal")
    
    plt.tight_layout()
    plot_path = base_dir / "results" / "bioenergy_spatial_analysis.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"Saved spatial maps to: {plot_path}")
    
    # Generate standalone net profit map
    fig, ax = plt.subplots(figsize=(8, 7))
    sc_prof = ax.scatter(
        comparison["longitude"], comparison["latitude"],
        c=comparison["net_profit_early"], cmap="RdYlGn",
        s=50, edgecolors="grey", linewidths=0.2
    )
    plt.colorbar(sc_prof, ax=ax, label="Net Profit ($/ha)")
    ax.set_title("Net Profit Map ($/ha) - Early Planting Scenario")
    ax.set_aspect("equal")
    plt.tight_layout()
    profit_plot_path = base_dir / "results" / "bioenergy_profit_map.png"
    fig.savefig(profit_plot_path, dpi=150)
    plt.close(fig)
    print(f"Saved net profit map to: {profit_plot_path}")
    
    print("\nAll bioenergy cover crop analysis complete.")

if __name__ == "__main__":
    main()
