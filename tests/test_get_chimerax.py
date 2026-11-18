"""The ChimeraX download helper, against a local fake of the UCSF download site."""

from __future__ import annotations

import hashlib
import http.server
import importlib.util
import json
import threading
import urllib.parse
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / ".agents" / "skills" / "ptmpsi" / "scripts" / "get_chimerax.py"
PAYLOAD = b"not really a disk image" * 1000
LINK = "9.9/mac_arm64/ChimeraX-9.9.dmg"


def _load():
    spec = importlib.util.spec_from_file_location("get_chimerax", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def site():
    state = {"posts": [], "sha": hashlib.sha256(PAYLOAD).hexdigest()}

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def _send(self, body, content_type):
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path.endswith("production.json"):
                build = {"name": "macOS (Apple Silicon)", "link": LINK, "size": "1 MiB", "sha256": state["sha"]}
                self._send(json.dumps(["9.9", "1 Jan 2027", {"mac_arm64": build}]).encode(), "application/json")
            elif "choice=Notified" in self.path:
                self._send(PAYLOAD, "application/octet-stream")
            else:
                gate = (f'<form method="POST"><input type="hidden" name="file" value="{LINK}"/>'
                        '<input type="submit" name="choice" value="Decline"/>'
                        '<input type="submit" name="choice" value="Accept"/></form>')
                self._send(gate.encode(), "text/html; charset=UTF-8")

        def do_POST(self):
            fields = urllib.parse.parse_qs(self.rfile.read(int(self.headers["Content-Length"])).decode())
            state["posts"].append(fields)
            if state.get("refresh"):
                refresh_url = f"/cgi-bin/secure/chimerax-get.py?ident=tok&file={urllib.parse.quote(LINK)}&choice=Notified"
                page = f'<meta http-equiv="Refresh" content="2;url={refresh_url}"/>'
                self._send(page.encode(), "text/html; charset=UTF-8")
            else:
                self._send(PAYLOAD, "application/octet-stream")

    server = http.server.HTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    state["url"] = f"http://127.0.0.1:{server.server_port}"
    yield state
    server.shutdown()


def test_without_the_flag_nothing_is_accepted_or_downloaded(site, tmp_path, capsys):
    module = _load()
    code = module.main(["--platform", "mac_arm64", "--base-url", site["url"], "-o", str(tmp_path)])
    assert code == 2
    assert site["posts"] == []
    assert list(tmp_path.iterdir()) == []
    assert "licence" in capsys.readouterr().out


def test_with_the_flag_it_accepts_downloads_and_verifies(site, tmp_path):
    module = _load()
    assert module.main(["--platform", "mac_arm64", "--accept-license", "--base-url", site["url"], "-o", str(tmp_path)]) == 0
    assert site["posts"] == [{"file": [LINK], "choice": ["Accept"]}]
    assert (tmp_path / "ChimeraX-9.9.dmg").read_bytes() == PAYLOAD


def test_meta_refresh_redirect_is_followed(site, tmp_path):
    site["refresh"] = True
    module = _load()
    assert module.main(["--platform", "mac_arm64", "--accept-license", "--base-url", site["url"], "-o", str(tmp_path)]) == 0
    assert (tmp_path / "ChimeraX-9.9.dmg").read_bytes() == PAYLOAD


def test_a_checksum_mismatch_deletes_the_download(site, tmp_path):
    site["sha"] = "0" * 64
    module = _load()
    with pytest.raises(SystemExit, match="mismatch"):
        module.main(["--platform", "mac_arm64", "--accept-license", "--base-url", site["url"], "-o", str(tmp_path)])
    assert list(tmp_path.iterdir()) == []


def test_list_and_unknown_platform(site, capsys):
    module = _load()
    assert module.main(["--list", "--base-url", site["url"]]) == 0
    assert "mac_arm64" in capsys.readouterr().out
    with pytest.raises(SystemExit, match="Unknown platform"):
        module.main(["--platform", "beos", "--base-url", site["url"]])
