#!/usr/bin/env python3
"""Download UCSF ChimeraX for this computer, after the user has agreed to its licence.

ChimeraX is free for academic, government, nonprofit and personal use under UCSF's non-commercial
licence (https://www.cgl.ucsf.edu/chimerax/docs/licensing.html). The download page asks the person
to press "Accept". This script presses it on their behalf, so it does nothing unless it is given
--accept-license, and the agent must have the user's yes before passing that flag.

    get_chimerax.py --list                       show the builds and the one chosen for this machine
    get_chimerax.py                              show the licence notice and stop (downloads nothing)
    get_chimerax.py --accept-license [-o DIR]    download, then check the SHA-256 the site publishes

It only downloads. It does not install anything.
"""

import argparse
import hashlib
import html
import json
import os
import platform
import re
import sys
import urllib.parse
import urllib.request

BASE = "https://www.cgl.ucsf.edu/chimerax"
LICENSE_URL = BASE + "/docs/licensing.html"


def default_platform():
    system, machine = platform.system(), platform.machine().lower()
    arm = machine in ("arm64", "aarch64")
    if system == "Darwin":
        return "mac_arm64" if arm else "macosx"
    if system == "Windows":
        return "windows"
    raise SystemExit("On Linux every ChimeraX package (.deb, .rpm, Flatpak) is built for a specific "
                     "system. Pick one with --platform (see --list). Installing it needs admin rights.")


def releases(base):
    with urllib.request.urlopen(f"{base}/data/release-info/production.json", timeout=30) as response:
        version, _date, builds = json.load(response)[:3]
    return version, builds


def download(base, link, destination, sha256):
    gate = f"{base}/cgi-bin/secure/chimerax-get.py"
    with urllib.request.urlopen(f"{gate}?file={urllib.parse.quote(link)}", timeout=30) as response:
        page = response.read().decode("utf-8", "replace")
    if not re.search(r'name="choice"\s+value="Accept"', page):
        raise SystemExit("The ChimeraX download page changed: no Accept button found. Download it in a browser.")
    hidden = dict(re.findall(r'<input type="hidden" name="(\w+)" value="([^"]*)"', page))
    body = urllib.parse.urlencode({**hidden, "choice": "Accept"}).encode()
    digest = hashlib.sha256()

    def _save(response):
        if "text/html" in response.headers.get("Content-Type", ""):
            raise SystemExit("ChimeraX returned a web page instead of the installer. Download it in a browser.")
        with open(destination, "wb") as out:
            for chunk in iter(lambda: response.read(1 << 20), b""):
                digest.update(chunk)
                out.write(chunk)

    with urllib.request.urlopen(urllib.request.Request(gate, data=body), timeout=60) as response:
        if "text/html" in response.headers.get("Content-Type", ""):
            refresh_page = response.read().decode("utf-8", "replace")
            match = re.search(r'<meta\s+http-equiv=["\']?Refresh["\']?\s+content=["\']?\d+;\s*url=([^"\'>]+)',
                              refresh_page, re.IGNORECASE)
            if not match:
                raise SystemExit("ChimeraX returned a web page instead of the installer. Download it in a browser.")
            next_url = urllib.parse.urljoin(gate, html.unescape(match.group(1)))
            with urllib.request.urlopen(next_url, timeout=120) as file_response:
                _save(file_response)
        else:
            _save(response)
    if digest.hexdigest() != sha256:
        os.remove(destination)
        raise SystemExit(f"SHA-256 mismatch for {destination}: got {digest.hexdigest()}, expected {sha256}. Deleted it.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--platform", help="build name from --list (default: chosen for this computer)")
    parser.add_argument("--list", action="store_true", help="list the builds and exit")
    parser.add_argument("--accept-license", action="store_true",
                        help="press Accept on the UCSF ChimeraX licence page on the user's behalf")
    parser.add_argument("-o", "--output-dir", default=".", help="where to save the installer")
    parser.add_argument("--base-url", default=BASE, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    version, builds = releases(args.base_url)
    if args.list:
        for name, build in builds.items():
            print(f"{name:14s} {build['size']:>11s}  {build['name']}")
        print(f"\nversion {version}; this computer: {default_platform_or_none()}")
        return 0

    name = args.platform or default_platform()
    if name not in builds:
        raise SystemExit(f"Unknown platform {name!r}. Choose from: {', '.join(builds)}")
    build = builds[name]

    if not args.accept_license:
        print(f"ChimeraX {version} ({build['name']}, {build['size']}) is released under UCSF's non-commercial "
              f"licence:\n  {LICENSE_URL}\n"
              "Downloading it means accepting that licence. Nothing was downloaded.\n"
              "If the user agrees, run again with --accept-license.")
        return 2

    os.makedirs(args.output_dir, exist_ok=True)
    destination = os.path.join(args.output_dir, os.path.basename(build["link"]))
    print(f"Accepting the UCSF ChimeraX licence ({LICENSE_URL}) on the user's behalf; downloading {destination}")
    download(args.base_url, build["link"], destination, build["sha256"])
    print(f"Downloaded and verified {destination}")
    return 0


def default_platform_or_none():
    try:
        return default_platform()
    except SystemExit:
        return "Linux: choose with --platform"


if __name__ == "__main__":
    sys.exit(main())
