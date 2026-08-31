"""Install a built wheel in a clean temporary environment and test outside the source tree."""
from pathlib import Path
import os
import subprocess
import tempfile
import venv

root = Path(__file__).resolve().parents[1]
wheels = list((root / "dist").glob("casein_ftir-*.whl"))
if len(wheels) != 1:
    raise SystemExit("Build exactly one wheel in dist/ first")
with tempfile.TemporaryDirectory(prefix="casein-wheel-") as temporary:
    folder = Path(temporary)
    env = folder / "venv"
    venv.EnvBuilder(with_pip=True).create(env)
    python = env / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    subprocess.run([str(python), "-m", "pip", "install", str(wheels[0])], cwd=folder, check=True)
    commands = [
        ["--version"],
        ["info"],
        ["simulate", "--output", "demo.csv"],
        ["analyze", "demo.csv", "--deconv", "--reference", "demo.csv", "--report", "report.html"],
    ]
    for command in commands:
        subprocess.run([str(python), "-m", "casein_ftir", *command], cwd=folder, check=True)
    assert (folder / "report.html").is_file()
print("Wheel installed and exercised outside the source tree.")
