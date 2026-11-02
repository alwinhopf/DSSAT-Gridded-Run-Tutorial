"""Exercise the full provider selection chain, not isolated branch bodies."""
import ast
import copy
from pathlib import Path
import subprocess
import pytest
from tests.helpers.discovery import find_rscript

ROOT = Path(__file__).resolve().parents[1]

@pytest.mark.parametrize("root", [ROOT])
@pytest.mark.parametrize("language", ["python", "R"])
def test_soil_routing(root, language, tmp_path):
    expected = {"SSURGO":1,"GNATSGO":1,"LUCAS":1,"SOILGRIDS_10K":2,
                "AGMIP":2,"SOILGRIDS_ONLINE":3,"POLARIS":4}
    if language == "python":
        tree = ast.parse((root / "dssat_main_pipeline.py").read_text(encoding="utf-8"))
        assignments = [n for n in tree.body if isinstance(n, ast.Assign) and
                       any(isinstance(t, ast.Name) and t.id == "_PER_POINT_SOIL" for t in n.targets)]
        stage = next(
            n for n in ast.walk(tree)
            if isinstance(n, ast.If)
            and ast.unparse(n.test) == "RUN_STEP_1_SOILS"
            and any(isinstance(c, ast.If) and "SOIL_SOURCE" in ast.unparse(c.test) for c in n.body)
        )
        chain = copy.deepcopy(next(n for n in stage.body if isinstance(n, ast.If) and "SOIL_SOURCE" in ast.unparse(n.test)))
        current = chain
        index = 1
        while True:
            current.body = ast.parse(f"selected = {index}").body
            if not current.orelse or not isinstance(current.orelse[0],ast.If):break
            current = current.orelse[0];index += 1
        code = compile(ast.fix_missing_locations(ast.Module(body=assignments+[chain],type_ignores=[])),"routing","exec")
        for source, branch in expected.items():
            env = dict(SOIL_SOURCE=source)
            exec(code,env)
            assert env["selected"] == branch,source
    else:
        rscript = find_rscript()
        if not rscript:pytest.skip("R unavailable")
        script = tmp_path / "routing.R"
        script.write_text(r'''args<-commandArgs(TRUE)
e<-new.env()
has_soil_source_if<-function(x) {
 is.call(x) && identical(x[[1]],as.name("if")) && grepl("SOIL_SOURCE",paste(deparse(x[[2]]),collapse=" "))
}
find_stage<-function(x) {
 if(missing(x)||!is.call(x))return(NULL)
 if(identical(x[[1]],as.name("if")) && identical(x[[2]],as.name("RUN_STEP_1_SOILS")) &&
    length(x)>=3 && is.call(x[[3]]) && any(vapply(as.list(x[[3]])[-1],has_soil_source_if,logical(1))))return(x)
 for(z in as.list(x)[-1]) {r<-find_stage(z);if(!is.null(r))return(r)}
 NULL
}
for(x in parse(args[1])) {
 if(is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name("PER_POINT_SOIL")))eval(x,e)
 stage<-find_stage(x);if(!is.null(stage))break
}
chain<-Filter(has_soil_source_if,as.list(stage[[3]])[-1])[[1]]
replace<-function(x,i=1L) {x[[3]]<-i;if(length(x)>3 && is.call(x[[4]]) && identical(x[[4]][[1]],as.name("if")))x[[4]]<-replace(x[[4]],i+1L);x}
chain<-replace(chain)
for(source in c("SSURGO","GNATSGO","LUCAS","SOILGRIDS_10K","AGMIP","SOILGRIDS_ONLINE","POLARIS")) {
 e$SOIL_SOURCE<-source
 cat(source,eval(chain,e),"\n")
}
''')
        run = subprocess.run([str(rscript),"--vanilla",str(script),str(root/"dssat_main_pipeline.R")],capture_output=True,text=True,timeout=30)
        assert run.returncode==0,run.stdout+run.stderr
        actual={line.split()[0]:int(line.split()[1]) for line in run.stdout.splitlines() if len(line.split())==2}
        assert actual==expected


@pytest.mark.parametrize("root", [ROOT])
@pytest.mark.parametrize("language", ["python", "R"])
def test_per_point_soil_id_na_fallback(root, language, tmp_path):
    """Per-point soil sources must backfill missing SOIL_ID rows from POINT_ID_COLUMN."""
    import pandas as pd
    if language == "python":
        tree = ast.parse((root / "dssat_main_pipeline.py").read_text(encoding="utf-8"))
        assignments = [
            n for n in tree.body
            if isinstance(n, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == "_PER_POINT_SOIL" for t in n.targets)
        ]
        fallback_if = next(
            n for n in ast.walk(tree)
            if isinstance(n, ast.If)
            and ast.unparse(n.test) == "'SOIL_ID' not in points.columns"
        )
        code = compile(
            ast.fix_missing_locations(ast.Module(body=assignments + [fallback_if], type_ignores=[])),
            "soil_id_fallback",
            "exec",
        )
        points = pd.DataFrame({"ID": ["00000001", "00000002"], "SOIL_ID": ["00000001", None]})
        env = dict(points=points, POINT_ID_COLUMN="ID", SOIL_SOURCE="SSURGO")
        exec(code, env)
        assert env["points"]["SOIL_ID"].tolist() == ["00000001", "00000002"]
    else:
        rscript = find_rscript()
        if not rscript:
            pytest.skip("R unavailable")
        script = tmp_path / "soil_id_fallback.R"
        script.write_text(r'''args<-commandArgs(TRUE)
e<-new.env()
for(x in parse(args[1])) {
 if(is.call(x) && identical(x[[1]],as.name("<-")) && identical(x[[2]],as.name("PER_POINT_SOIL")))eval(x,e)
 if(is.call(x) && identical(x[[1]],as.name("if")) && grepl('SOIL_ID.*names\\(points\\)',paste(deparse(x[[2]]),collapse=" "))) {
  fallback_if<-x; break
 }
}
e$POINT_ID_COLUMN<-"ID"
e$SOIL_SOURCE<-"SSURGO"
e$points<-data.frame(ID=c("00000001","00000002"),SOIL_ID=c("00000001",NA_character_),stringsAsFactors=FALSE)
eval(fallback_if,e)
cat(paste(e$points$SOIL_ID,collapse=","),"\n")
''')
        run = subprocess.run([str(rscript), "--vanilla", str(script), str(root / "dssat_main_pipeline.R")],
                             capture_output=True, text=True, timeout=30)
        assert run.returncode == 0, run.stdout + run.stderr
        assert run.stdout.strip() == "00000001,00000002"
