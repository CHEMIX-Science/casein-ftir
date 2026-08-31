"""Install a built wheel in a clean temporary environment and test outside the source tree."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile
import venv

root = Path(__file__).resolve().parents[1]
wheels = [Path(sys.argv[1]).resolve()] if len(sys.argv) > 1 else list((root / "dist").glob("casein_ftir-*.whl"))
if len(wheels) != 1:
    raise SystemExit("Build exactly one wheel in dist/ first")
with tempfile.TemporaryDirectory(prefix="casein-wheel-") as temporary:
    folder = Path(temporary)
    env = folder / "venv"
    venv.EnvBuilder(with_pip=True).create(env)
    python = env / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    subprocess.run([str(python), "-m", "pip", "install", str(wheels[0])], cwd=folder, check=True)
    subprocess.run([str(python), "-c",
                    "from casein_ftir.defaults import get_galalithe_dft_reference; "
                    "r = get_galalithe_dft_reference(); "
                    "assert len(r.wavenumber) == 1801; "
                    "assert r.metadata['type'] == 'dft_packaged_reference'; "
                    "assert r.metadata['license'] == 'MIT'"], cwd=folder, check=True)
    commands = [
        ["--version"],
        ["info"],
        ["simulate", "--output", "demo.csv"],
        ["analyze", "demo.csv", "--deconv", "--reference", "demo.csv", "--report", "report.html"],
        ["analyze", "demo.csv", "--use-default-galalithe-ref", "--no-preprocess-ref", "--report", "dft-report.html"],
    ]
    for command in commands:
        subprocess.run([str(python), "-m", "casein_ftir", *command], cwd=folder, check=True)
    assert (folder / "report.html").is_file()
    assert "fragment_reference" in (folder / "dft-report.html").read_text(encoding="utf-8").lower()
print("Wheel installed and exercised outside the source tree.")
