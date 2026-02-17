import os
import requests
from pathlib import Path


def human_size(size):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024:
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{size:.1f}PB"


def print_path(root: Path):
    if not root.exists():
        print(f"{root} does not exist")
        return

    print(root.name)
    _print_tree(root, prefix="")


def _print_tree(root: Path, prefix: str):
    # skip folder mmseqs2_tmp
    if root.name == 'mmseqs2_tmp':
        return
    entries = sorted(root.iterdir(), key=lambda p: (p.is_file(), p.name.lower()))
    for i, path in enumerate(entries):
        is_last = i == len(entries) - 1
        connector = "└── " if is_last else "├── "

        if path.is_file():
            size = human_size(path.stat().st_size)
            print(prefix + connector + f"{path.name} ({size})")
        else:
            print(prefix + connector + path.name)
            extension = "    " if is_last else "│   "
            _print_tree(path, prefix + extension)


def upload_blob_file(filepath, token, shock_url, hs_client):
    """Upload a file to Shock and get handle.

    Args:
        filepath: Path to file to upload

    Returns:
        Tuple of (shock_id, handle_id)
    """

    # Upload to Shock
    headers = {"Authorization": "OAuth " + token}

    # Get file size for Content-Length
    file_size = os.path.getsize(filepath)
    filename = os.path.basename(filepath)

    with open(filepath, "rb") as f:
        # Use multipart form with file and explicit content-type/size
        # The tuple format is (filename, fileobj, content_type, headers)
        files = {
            "upload": (
                filename,
                f,
                "application/octet-stream",
                {"Content-Length": str(file_size)},
            )
        }

        r = requests.post(
            shock_url + "/node",
            headers=headers,
            files=files,
            allow_redirects=True,
        )

        if not r.ok:
            error_msg = r.text
            try:
                error_data = r.json()
                error_msg = error_data.get("error", [r.text])[0]
            except:
                pass
            raise RuntimeError(f"Failed to upload file to Shock: {error_msg}")

        shock_node = r.json()["data"]
        shock_id = shock_node["id"]

    # Create handle
    hs = hs_client
    handle = hs.persist_handle(
        {
            "id": shock_id,
            "type": "shock",
            "url": shock_url,
        }
    )
    handle_id = handle

    print(f"File uploaded to Shock: {shock_id}, Handle: {handle_id}")
    return shock_id, handle_id