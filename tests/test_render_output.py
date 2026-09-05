"""Tests for the HTML renderer and ``write_outputs`` (WP2, part 2).

Golden tests: for each of the four fixture payloads, the rendered Markdown and HTML are compared
byte-for-byte against a committed file under ``tests/fixtures``. Regenerate the committed files by
running this module with the environment variable ``PBCHECK_UPDATE_GOLDEN=1`` set; a pytest CLI
option (``--update-golden``) cannot be registered for this instead, because ``tests/conftest.py``
belongs to another work package during this release and this module must not add one.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from pbcheck.render import render_html, render_markdown, write_outputs

FIXTURES_DIR = Path(__file__).parent / "fixtures"
FIXTURE_NAMES = ("complete", "naive_only", "design_only", "coarse_null")

_UPDATE_GOLDEN = os.environ.get("PBCHECK_UPDATE_GOLDEN") == "1"


def _load_fixture(name: str) -> dict:
    path = FIXTURES_DIR / f"audit_payload_{name}.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _check_golden(golden_path: Path, rendered: str) -> None:
    if _UPDATE_GOLDEN:
        with open(golden_path, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(rendered)
        return
    expected = golden_path.read_text(encoding="utf-8")
    assert rendered == expected, f"rendered output no longer matches {golden_path}"


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_golden_markdown(name: str) -> None:
    payload = _load_fixture(name)
    rendered = render_markdown(payload)
    _check_golden(FIXTURES_DIR / f"audit_payload_{name}.md", rendered)


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_golden_html(name: str) -> None:
    payload = _load_fixture(name)
    rendered = render_html(payload)
    _check_golden(FIXTURES_DIR / f"audit_payload_{name}.html", rendered)


def test_html_escapes_values() -> None:
    payload = _load_fixture("complete")
    rendered = render_html(payload)
    assert "<b>" not in rendered
    assert "&lt;b&gt;" in rendered


def test_html_has_no_script_or_external_asset() -> None:
    payload = _load_fixture("complete")
    rendered = render_html(payload)
    assert "<script" not in rendered.lower()
    assert "http://" not in rendered
    assert "https://" not in rendered
    assert "<link" not in rendered.lower()


def test_write_outputs_refuses_overwrite(tmp_path: Path) -> None:
    payload = _load_fixture("naive_only")
    write_outputs(payload, tmp_path)
    with pytest.raises(FileExistsError, match="pbcheck_audit.json"):
        write_outputs(payload, tmp_path)


def test_write_outputs_utf8(tmp_path: Path) -> None:
    payload = _load_fixture("design_only")
    out_dir = tmp_path / "étude-café"
    written = write_outputs(payload, out_dir)

    assert set(written) == {"json", "md", "html"}
    for fmt, path in written.items():
        assert path.parent == out_dir
        content = path.read_bytes()
        assert b"\r\n" not in content
        text = content.decode("utf-8")
        if fmt == "json":
            assert json.loads(text) == payload
        elif fmt == "md":
            assert text == render_markdown(payload)
        else:
            assert text == render_html(payload)
