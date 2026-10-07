import importlib.util
from pathlib import Path
import pandas as pd


def test_hpc_matches_engine_missing_stocks_and_fluxes(tmp_path):
    spec=importlib.util.spec_from_file_location('mpi_review',Path(__file__).parents[1]/'hpc/dssat_mpi_runner.py')
    h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    from dssatengine.engine import _read_csv_safe,_merge_supplemental,_build_result_rows
    pd.DataFrame([{'RUNNO':1,'TRNO':3,'PDAT':2020001,'HYEAR':2020,'HWAM':-99,'CO2EM':12}]).to_csv(tmp_path/'summary.csv',index=False)
    pd.DataFrame({'RUN':[1,1],'SOMCT':[-99,30]}).to_csv(tmp_path/'soilorg.csv',index=False)
    mpi=h.process_dssat_outputs_no_pandas(str(tmp_path),'p',3,{})[0]
    summary=_read_csv_safe(tmp_path/'summary.csv')
    local=_build_result_rows('p',summary,_merge_supplemental(str(tmp_path),summary[['RUNNO']])).iloc[0]
    for name in ['final_grain_kg_ha','soil_organic_carbon_start_kg_C_ha','soil_organic_carbon_delta_kg_C_ha']:
        assert mpi[name] is None and pd.isna(local[name])
    assert mpi['cumulative_net_co2_emissions_kg_CO2_ha']==local['cumulative_net_co2_emissions_kg_CO2_ha']==44
    assert mpi['treatment']==3
