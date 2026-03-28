import ternviz
from ternviz.lib import _macos_vmd_env, resolve_vmd, vmd_script


def test_ternviz():
    assert ternviz.__version__


def test_macos_vmd_env_from_bundle_binary(monkeypatch, tmp_path):
    binary = tmp_path / "VMD2b1.app" / "Contents" / "vmd2b1" / "vmd_MACOSXARM64"
    binary.parent.mkdir(parents=True)
    binary.touch()
    (binary.parent / "stride_MACOSXARM64").touch()
    (binary.parent / "surf_MACOSXARM64").touch()
    (binary.parent / "tachyon_MACOSXARM64").touch()
    monkeypatch.setenv("PYTHONPATH", "/existing/python")

    env = _macos_vmd_env(binary)

    assert env["VMDDIR"] == str(binary.parent)
    assert env["TCL_LIBRARY"] == str(binary.parent / "scripts" / "tcl")
    assert env["PYTHONPATH"] == (
        f"{binary.parent / 'scripts' / 'python'}:/existing/python"
    )
    assert env["STRIDE_BIN"] == str(binary.parent / "stride_MACOSXARM64")
    assert env["SURF_BIN"] == str(binary.parent / "surf_MACOSXARM64")
    assert env["TACHYON_BIN"] == str(binary.parent / "tachyon_MACOSXARM64")


def test_resolve_vmd_accepts_app_bundle(tmp_path):
    app_path = tmp_path / "VMD2b1.app"
    binary = app_path / "Contents" / "vmd2b1" / "vmd_MACOSXARM64"
    binary.parent.mkdir(parents=True)
    binary.touch()

    resolved, env = resolve_vmd(str(app_path))

    assert resolved == str(binary)
    assert env["VMDDIR"] == str(binary.parent)


def test_resolve_vmd_accepts_existing_binary(tmp_path):
    binary = tmp_path / "vmd"
    binary.touch()

    resolved, env = resolve_vmd(str(binary))

    assert resolved == str(binary)
    assert env == {}


def test_vmd_script_can_strip_display_resize():
    script = vmd_script(
        800,
        "sample",
        script_name="render.vmd",
        strip_display_resize=True,
    )

    assert all("display resize" not in line for line in script)
