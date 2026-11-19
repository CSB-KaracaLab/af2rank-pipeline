"""Download AF2Rank weights from DeepMind's dated AlphaFold archives.

Each archive is about 5.3 GB; this helper uses HTTP byte ranges to fetch only
the requested weight file and checks its SHA-256. Run as a script, it fetches
the default params_model_2_ptm.npz. No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import urllib.request
import zipfile
from pathlib import Path


def archive_url(release: str) -> str:
    """URL of a dated parameter archive: 2021-10-27 (AF2.1), 2022-03-02 (AF2.2), 2022-12-06 (AF2.3)."""
    return f"https://storage.googleapis.com/alphafold/alphafold_params_{release}.tar"


ARCHIVE_URL = archive_url("2022-12-06")
WEIGHT_NAME = "params_model_2_ptm.npz"
WEIGHT_SHA256 = "23645d9a82c4af2ed54cd48a7b3c1c2575dc6aa9fe931adb4d7203ca5f0dc398"


def get_range(start: int, end: int, url: str = ARCHIVE_URL):
    request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
    response = urllib.request.urlopen(request, timeout=60)
    expected = f"bytes {start}-{end}/"
    if response.status != 206 or not response.headers.get("Content-Range", "").startswith(expected):
        response.close()
        raise RuntimeError("The AlphaFold archive server did not honor the byte-range request")
    return response


def locate_weight(weight_name: str = WEIGHT_NAME, url: str = ARCHIVE_URL) -> tuple[int, int]:
    offset = 0
    for _ in range(100):
        with get_range(offset, offset + 511, url) as response:
            header = response.read()
        if len(header) != 512 or header[257:262] != b"ustar":
            raise RuntimeError(f"Invalid archive header at byte {offset}")
        # The 2021-10-27 archive stores members as ./params_*.npz.
        name = header[:100].split(b"\0", 1)[0].decode("ascii").removeprefix("./")
        size_field = header[124:136].split(b"\0", 1)[0].strip()
        size = int(size_field, 8)
        if name == weight_name:
            return offset + 512, size
        offset += 512 + ((size + 511) // 512) * 512
    raise RuntimeError(f"{weight_name} was not found in the official archive {url}")


def verify_weight(path: Path, sha256: str = WEIGHT_SHA256) -> None:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != sha256:
        raise RuntimeError(f"Unexpected checksum for {path}; remove the file and retry")
    with zipfile.ZipFile(path) as archive:
        corrupt_member = archive.testzip()
    if corrupt_member is not None:
        raise RuntimeError(f"Corrupt parameter archive member: {corrupt_member}")


def download_weight(params_dir: Path, weight_name: str = WEIGHT_NAME,
                    sha256: str = WEIGHT_SHA256, url: str = ARCHIVE_URL) -> Path:
    """Fetch one weight file into params_dir, or verify the copy already there."""
    params_dir = Path(params_dir).expanduser().resolve()
    params_dir.mkdir(parents=True, exist_ok=True)
    destination = params_dir / weight_name
    if destination.is_file():
        verify_weight(destination, sha256)
        print(f"Verified existing weight: {destination}")
        return destination

    start, size = locate_weight(weight_name, url)
    partial = destination.with_suffix(".npz.part")
    try:
        with get_range(start, start + size - 1, url) as response, partial.open("wb") as output:
            shutil.copyfileobj(response, output)
        if partial.stat().st_size != size:
            raise RuntimeError(f"Incomplete download: expected {size} bytes")
        verify_weight(partial, sha256)
        partial.replace(destination)
    finally:
        partial.unlink(missing_ok=True)
    print(f"Downloaded and verified: {destination}")
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("params_dir", type=Path, help="Destination params directory")
    args = parser.parse_args()
    download_weight(args.params_dir)


if __name__ == "__main__":
    main()
