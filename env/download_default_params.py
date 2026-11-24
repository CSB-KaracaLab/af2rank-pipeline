"""Download AF2Rank weights from DeepMind's dated AlphaFold archives.

Each archive is about 5.3 GB; this helper uses HTTP byte ranges to fetch only
the requested weight file and checks its SHA-256. Run as a script, it fetches
the default params_model_2_ptm.npz, or the weights named with --weight
(repeatable), e.g. --weight model_1_multimer_v2. No third-party Python packages
are required.
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

# The 17 template-capable weights: name -> (archive release, SHA-256 of params_<name>.npz).
WEIGHTS = {
    "model_1_ptm": ("2022-12-06", "5e564f79af5bcd54ccef6e2a6bb0ff01015d01650ebc41d4575e35f0de9ecc84"),
    "model_2_ptm": ("2022-12-06", WEIGHT_SHA256),
    "model_1_multimer": ("2021-10-27", "e12c732731fa2f4957979180c01f78be4bd20c64a615ae59902d811ccb6cfb17"),
    "model_2_multimer": ("2021-10-27", "a15017e8665da8811611381fb8c8e8e23d7244f7c72b25742b7674219e014ec0"),
    "model_3_multimer": ("2021-10-27", "3d9738bc14bfba44d4584d7afc37bcb52e6cb83412eecb23c9e2e45619bfec85"),
    "model_4_multimer": ("2021-10-27", "a6b00af878c13594023c7920005c5be7e470343df6745e702faec65a19d44a22"),
    "model_5_multimer": ("2021-10-27", "371dfdea9efdfa8b70d40a267d61286e6b1b85c822b7bfaba292c2e960f43b2d"),
    "model_1_multimer_v2": ("2022-03-02", "69e20de3e265c22d4b1616a210f0715cdc969b8c1f59d7ff882635f1e6c986d6"),
    "model_2_multimer_v2": ("2022-03-02", "7d78b62d75fed71b0688e4ea434fbb385fad29b2544dfbdc0457ed4a95990e27"),
    "model_3_multimer_v2": ("2022-03-02", "349e8b4ec8758f9348a05418e6c0105ccafd811641451baa2250ab0209e6bf32"),
    "model_4_multimer_v2": ("2022-03-02", "5f68bbf95ad633734b8ef8a3633d5d2b3d821bfe923bf1e72b8d06177ea7fb69"),
    "model_5_multimer_v2": ("2022-03-02", "73d4ed49917a5058e6d113ecf5af31799db7fd5a3c051539d292064d11e2e53b"),
    "model_1_multimer_v3": ("2022-12-06", "611da8fc7478928f68de12e8b226260ef1f4ce62bcc29b008572e52f4f212959"),
    "model_2_multimer_v3": ("2022-12-06", "51362b0844382ae0f5720c59b81dd13a43eea40fbf9995dd2573bdab88865378"),
    "model_3_multimer_v3": ("2022-12-06", "46d9bcad288edc7ad5a6362ee8e5f84307a69712e00e6c36b1ef9daf96ebc9ce"),
    "model_4_multimer_v3": ("2022-12-06", "59bdabd2d69c07fe26b37882d544acbd1b9f89f196828f4220da49e0610b572c"),
    "model_5_multimer_v3": ("2022-12-06", "917742be5a105d6b80f13f1f13f20459f27ec3fdcd34ea088b359f4502d6177f"),
}


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
    parser.add_argument("--weight", action="append", choices=sorted(WEIGHTS), metavar="NAME",
                        help="Weight to fetch, e.g. model_1_multimer_v2 (repeatable; default model_2_ptm)")
    args = parser.parse_args()
    for name in args.weight or ["model_2_ptm"]:
        release, sha256 = WEIGHTS[name]
        download_weight(args.params_dir, f"params_{name}.npz", sha256, archive_url(release))


if __name__ == "__main__":
    main()
