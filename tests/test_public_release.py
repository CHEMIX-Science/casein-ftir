"""Regression tests for public inputs, reproducibility and report safety."""
import numpy as np
import pytest
from casein_ftir.io_module import Spectrum, load_spectrum, save_spectrum_csv
from casein_ftir.compare import resample_to_common_axis
from casein_ftir.preprocess import baseline_als, smooth_savgol
from casein_ftir.cli import main
from casein_ftir.report import html_report, text_report
from casein_ftir.reference_module import load_reference, compare_to_reference
from casein_ftir.io_gaussian import gaussian_log_to_spectrum
from casein_ftir.simulate import simulate_casein, simulate_galalithe


def test_unsorted_axis_preserves_pairs():
    spec = Spectrum([1000, 3000, 2000], [1, 3, 2])
    np.testing.assert_array_equal(spec.wavenumber, [3000, 2000, 1000])
    np.testing.assert_array_equal(spec.absorbance, [3, 2, 1])


@pytest.mark.parametrize("x,y", [([1, 1], [2, 3]), ([1, 2], [np.nan, 1]), ([1, np.inf], [1, 2]), ([[1,2]], [[1,2]])])
def test_invalid_spectrum_rejected(x, y):
    with pytest.raises(ValueError):
        Spectrum(x, y)


@pytest.mark.parametrize("body", ["wn,ab\n3000,0.1\n2000,bad\n1000,0.3", "3000,0.1\n2000,inf", "3000,0.1,extra\n2000,0.2,extra", "3000,0.1\n3000,0.2"])
def test_csv_never_silently_discards_bad_rows(tmp_path, body):
    path=tmp_path / "bad.csv"
    path.write_text(body, encoding="utf-8")
    with pytest.raises(ValueError):
        load_spectrum(path)


@pytest.mark.parametrize("sep", [",", ";", "\t", " "])
def test_csv_accepts_supported_delimiters(tmp_path, sep):
    p=tmp_path / "spectrum.csv"
    p.write_text(sep.join(["wn", "ab"])+"\n"+sep.join(["2000", "0.5"])+"\n"+sep.join(["1000", "0.2"]), encoding="utf-8")
    assert load_spectrum(p).n_points == 2


def test_csv_dft_type_survives_export(tmp_path, clean_casein_spectrum):
    # Synthetic fixture tests metadata transport; no real DFT data is used.
    source=clean_casein_spectrum.copy()
    source.metadata["type"]="dft_simulated"
    p=tmp_path / "ref.csv"
    save_spectrum_csv(source, p)
    ref=load_reference(str(p), preprocess=False)
    result=compare_to_reference(clean_casein_spectrum, ref)
    assert result["overall_quality"] == "fragment_reference"
    assert "confirm the crosslink" not in result["recommended_action"]


def test_html_escapes_labels_and_hides_private_path(tmp_path):
    spec=Spectrum([3000,2000,1000],[0.1,0.5,0.2],name="<script>alert(1)</script>",source_file=str(tmp_path / "private" / "sample.csv"))
    html=html_report(spec,ratios={"<img src=x onerror=alert(1)>":1.0})
    assert "<script>" not in html and "<img src=x" not in html
    assert "&lt;script&gt;" in html
    assert str(tmp_path) not in html
    assert str(tmp_path) not in text_report(spec)
    path=tmp_path / "export.csv"
    save_spectrum_csv(spec,path)
    assert str(tmp_path) not in path.read_text(encoding="utf-8")


def test_no_overlap_fails():
    with pytest.raises(ValueError, match="overlapping"):
        resample_to_common_axis(Spectrum([4000,3000],[0,1]),Spectrum([2000,1000],[0,1]))


def test_common_axis_never_extrapolates():
    a,b=resample_to_common_axis(Spectrum([10,3],[1,0]),Spectrum([9,2],[1,0]),resolution=4)
    assert a.wavenumber.min() >= 3
    assert a.wavenumber.max() <= 9


def test_short_and_irregular_smoothing_rejected():
    with pytest.raises(ValueError):
        smooth_savgol(Spectrum([3,2,1],[1,2,1]))
    with pytest.raises(ValueError, match="evenly"):
        smooth_savgol(Spectrum(np.arange(20)**2,np.ones(20)))


def test_constant_als_does_not_produce_nan():
    result=baseline_als(Spectrum(np.arange(50),np.zeros(50)))
    assert np.isfinite(result.absorbance).all()


@pytest.mark.parametrize("argv", [
    ["analyze","missing.csv","--deconv-export","x.csv"],
    ["analyze","missing.csv","--subtract-ref","x.csv"],
    ["analyze","missing.csv","--reference","ref.csv","--use-default-galalithe-ref"],
    ["gaussian","convert"],
    ["gaussian","convert","input.log"],
])
def test_cli_rejects_incomplete_commands(argv):
    with pytest.raises(SystemExit) as error:
        main(argv)
    assert error.value.code == 2


def test_requested_fit_failure_returns_error(tmp_path):
    p=tmp_path / "short.csv"
    p.write_text("1650,0.2\n1648,0.3\n1646,0.2",encoding="utf-8")
    with pytest.raises(SystemExit) as error:
        main(["analyze",str(p),"--no-preprocess","--deconv"])
    assert error.value.code == 2


def test_gaussian_missing_intensities_not_fabricated(tmp_path):
    p=tmp_path / "missing.log"
    p.write_text("Frequencies -- 1000 1500 1700",encoding="utf-8")
    with pytest.raises(ValueError,match="intensities"):
        gaussian_log_to_spectrum(p)


@pytest.mark.parametrize("kwargs", [{"resolution":0}, {"wn_min":4000,"wn_max":600}, {"noise_level":-1}, {"species":"unknown"}])
def test_invalid_simulation_parameters(kwargs):
    with pytest.raises(ValueError):
        simulate_casein(**kwargs)


def test_synthetic_crosslink_parameter_is_bounded():
    with pytest.raises(ValueError):
        simulate_galalithe(crosslink_degree=1.5)
