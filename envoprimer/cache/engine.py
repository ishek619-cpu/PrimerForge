"""
EnvoPrimer Persistent Cache Engine.

Caches every expensive stage of the pipeline.

Author: EnvoPrimer
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path


class CacheEngine:

    CACHE_VERSION = "1.0"

    def __init__(
        self,
        cache_dir: str = "cache",
    ):

        self.root = Path(cache_dir)

        self.root.mkdir(
            parents=True,
            exist_ok=True,
        )

    ####################################################################
    # Cache Key
    ####################################################################

    def key(
        self,
        species: str,
        marker: str,
    ) -> str:

        text = (
            f"{species}|{marker}|{self.CACHE_VERSION}"
        )

        return hashlib.sha256(
            text.encode()
        ).hexdigest()[:20]

    ####################################################################
    # Cache Directory
    ####################################################################

    def directory(
        self,
        species: str,
        marker: str,
    ) -> Path:

        folder = self.root / self.key(
            species,
            marker,
        )

        folder.mkdir(
            parents=True,
            exist_ok=True,
        )

        return folder

    ####################################################################
    # Stage Directories
    ####################################################################

    def downloads(
        self,
        species,
        marker,
    ):

        d = self.directory(
            species,
            marker,
        ) / "downloads"

        d.mkdir(
            parents=True,
            exist_ok=True,
        )

        return d

    def cleaned(
        self,
        species,
        marker,
    ):

        d = self.directory(
            species,
            marker,
        ) / "cleaned"

        d.mkdir(
            parents=True,
            exist_ok=True,
        )

        return d

    def alignments(
        self,
        species,
        marker,
    ):

        d = self.directory(
            species,
            marker,
        ) / "alignments"

        d.mkdir(
            parents=True,
            exist_ok=True,
        )

        return d

    def diagnostics(
        self,
        species,
        marker,
    ):

        d = self.directory(
            species,
            marker,
        ) / "diagnostics"

        d.mkdir(
            parents=True,
            exist_ok=True,
        )

        return d

    def primer3(
        self,
        species,
        marker,
    ):

        d = self.directory(
            species,
            marker,
        ) / "primer3"

        d.mkdir(
            parents=True,
            exist_ok=True,
        )

        return d

    def reports(
        self,
        species,
        marker,
    ):

        d = self.directory(
            species,
            marker,
        ) / "reports"

        d.mkdir(
            parents=True,
            exist_ok=True,
        )

        return d

    ####################################################################
    # Generic File Helpers
    ####################################################################

    def save_file(
        self,
        source,
        destination,
    ):

        source = Path(source)
        destination = Path(destination)

        destination.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        shutil.copy2(
            source,
            destination,
        )

    def load_file(
        self,
        path,
    ):

        path = Path(path)

        if path.exists():
            return path

        return None

    ####################################################################
    # JSON
    ####################################################################

    def save_json(
        self,
        file,
        data,
    ):

        file = Path(file)

        file.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with open(
            file,
            "w",
            encoding="utf-8",
        ) as handle:

            json.dump(
                data,
                handle,
                indent=4,
            )

    def load_json(
        self,
        file,
    ):

        file = Path(file)

        if not file.exists():

            return None

        with open(
            file,
            encoding="utf-8",
        ) as handle:

            return json.load(handle)

    ####################################################################
    # Pipeline Status
    ####################################################################

    def mark_complete(
        self,
        species,
        marker,
        metadata,
    ):

        self.save_json(

            self.directory(
                species,
                marker,
            ) / "complete.json",

            metadata,

        )

    def metadata(
        self,
        species,
        marker,
    ):

        return self.load_json(

            self.directory(
                species,
                marker,
            ) / "complete.json"

        )

    def exists(
        self,
        species,
        marker,
    ):

        return (

            self.directory(
                species,
                marker,
            ) / "complete.json"

        ).exists()

    ####################################################################
    # Clear Cache
    ####################################################################

    def clear(
        self,
        species,
        marker,
    ):

        folder = self.directory(
            species,
            marker,
        )

        if folder.exists():

            shutil.rmtree(folder)
