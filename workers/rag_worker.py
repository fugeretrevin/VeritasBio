"""
workers/rag_worker.py — DOMAIN: Science & Math Tools (Pillar 1: Biological Priors).

Purpose:
    Retrieves and parses verified biological metadata from UniProt (rest.uniprot.org)
    to provide the Strategic Planner (Meta Llama 3.1) with real functional constraints,
    chromophore/active-site coordinates, and mutagenesis context.

Invariants (ARCHITECTURE.md §3.1, §7; INNOVATION_ROADMAP Phase 2):
    1. Biological Grounding: Extracts real sequence, chromophore residues, and structural features.
    2. Zero-Cost / Local Resilience: Includes bundled offline cache (`workers/data/P42212.json`)
       for avGFP so offline testing and development never stall.
    3. Fail-Closed Parsing: Returns validated Pydantic models; rejects empty or malformed sequences.
"""

from __future__ import annotations

import json
import logging
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

DATA_DIR: Path = Path(__file__).resolve().parent / "data"
DEFAULT_ACCESSION: str = "P42212"  # avGFP (Aequorea victoria)
UNIPROT_REST_BASE: str = "https://rest.uniprot.org/uniprotkb"


class FeatureAnnotation(BaseModel):
    """Specific biological site or region annotation from UniProt."""
    feature_type: str = Field(..., description="e.g. Site, Region, Mutagenesis, Beta strand")
    description: str = Field(default="", description="Functional description or phenotypic effect")
    start: int = Field(..., ge=1, description="1-indexed start residue")
    end: int = Field(..., ge=1, description="1-indexed end residue")


class ProteinMetadata(BaseModel):
    """Validated biological metadata bundle used by the Strategic Planner."""
    accession: str
    entry_name: str
    protein_name: str
    organism: str
    sequence: str
    length: int
    chromophore_residues: List[int] = Field(default_factory=list)
    functional_features: List[FeatureAnnotation] = Field(default_factory=list)
    mutagenesis_hotspots: List[FeatureAnnotation] = Field(default_factory=list)

    def to_llm_summary(self) -> str:
        """
        Formats metadata into a concise, high-density scientific context
        block for the LLM prompt.
        """
        lines = [
            f"# Target Protein: {self.protein_name} ({self.entry_name})",
            f"- UniProt Accession: {self.accession}",
            f"- Organism: {self.organism}",
            f"- Sequence Length: {self.length} amino acids",
            f"- WT Sequence:\n{self.sequence}",
            "",
            "## Key Functional & Catalytic Sites:",
        ]

        if self.chromophore_residues:
            chrom_str = ", ".join(f"Residue {p} ({self.sequence[p-1]})" for p in self.chromophore_residues if p <= self.length)
            lines.append(f"- Core Chromophore Triad: {chrom_str}")

        lines.append("\n## Mutagenesis & Spectral Hotspots:")
        for mut in self.mutagenesis_hotspots[:10]:
            lines.append(f"- Pos {mut.start}-{mut.end}: {mut.description}")

        lines.append("\n## Key Structural Barrel Elements:")
        for feat in self.functional_features[:6]:
            lines.append(f"- {feat.feature_type} ({feat.start}-{feat.end}): {feat.description or 'Beta-canister scaffold'}")

        return "\n".join(lines)


def _load_cached_json(accession: str) -> Optional[Dict[str, Any]]:
    cache_path = DATA_DIR / f"{accession}.json"
    if cache_path.exists():
        try:
            return json.loads(cache_path.read_text(encoding="utf-8"))
        except Exception as exc:
            logger.warning("Failed to read local cache %s: %s", cache_path, exc)
    return None


def _save_cached_json(accession: str, data: Dict[str, Any]) -> None:
    try:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        cache_path = DATA_DIR / f"{accession}.json"
        cache_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        logger.debug("Saved UniProt cache to %s", cache_path)
    except Exception as exc:
        logger.warning("Failed to write UniProt cache for %s: %s", accession, exc)


def parse_uniprot_response(raw_data: Dict[str, Any]) -> ProteinMetadata:
    """Parses raw UniProt JSON payload into a structured ProteinMetadata model."""
    accession = raw_data.get("primaryAccession", "UNKNOWN")
    entry_name = raw_data.get("uniProtkbId", accession)
    
    # Description parsing
    desc_section = raw_data.get("proteinDescription", {})
    recommended_name = desc_section.get("recommendedName", {}).get("fullName", {}).get("value")
    protein_name = recommended_name or "Green fluorescent protein"

    organism = raw_data.get("organism", {}).get("scientificName", "Unknown organism")
    
    # Sequence parsing
    seq_data = raw_data.get("sequence", {})
    sequence = seq_data.get("value", "").strip().upper()
    length = seq_data.get("length", len(sequence))

    if not sequence:
        raise ValueError(f"No sequence found in UniProt record for {accession}")

    features = raw_data.get("features", [])
    functional_features: List[FeatureAnnotation] = []
    mutagenesis_hotspots: List[FeatureAnnotation] = []
    chromophore_residues: List[int] = []

    for f in features:
        ftype = f.get("type", "")
        fdesc = f.get("description", "")
        loc = f.get("location", {})
        start_val = loc.get("start", {}).get("value")
        end_val = loc.get("end", {}).get("value")

        if start_val is None or end_val is None:
            continue

        annot = FeatureAnnotation(
            feature_type=ftype,
            description=fdesc,
            start=int(start_val),
            end=int(end_val),
        )

        if "chromophore" in fdesc.lower() or "fluorophore" in fdesc.lower():
            for p in range(int(start_val), int(end_val) + 1):
                if p not in chromophore_residues:
                    chromophore_residues.append(p)

        if ftype == "Mutagenesis":
            mutagenesis_hotspots.append(annot)
        elif ftype in ("Site", "Region", "Active site", "Beta strand"):
            functional_features.append(annot)

    # avGFP default chromophore fallback (65, 66, 67) if not explicitly tagged
    if accession == "P42212" and not chromophore_residues:
        chromophore_residues = [65, 66, 67]

    return ProteinMetadata(
        accession=accession,
        entry_name=entry_name,
        protein_name=protein_name,
        organism=organism,
        sequence=sequence,
        length=length,
        chromophore_residues=sorted(chromophore_residues),
        functional_features=functional_features,
        mutagenesis_hotspots=mutagenesis_hotspots,
    )


def fetch_protein_metadata(
    accession: str = DEFAULT_ACCESSION,
    use_cache: bool = True,
    timeout: float = 8.0,
) -> ProteinMetadata:
    """
    Fetches protein metadata for the given UniProt accession.
    
    Falls back gracefully to the local offline cache if network is unavailable
    or request times out.
    """
    clean_acc = accession.strip().upper()

    # Step 1: If caching requested, check local disk first
    if use_cache:
        cached = _load_cached_json(clean_acc)
        if cached:
            logger.info("Loaded UniProt metadata for %s from local cache.", clean_acc)
            return parse_uniprot_response(cached)

    # Step 2: Hit UniProt REST API
    url = f"{UNIPROT_REST_BASE}/{clean_acc}.json"
    logger.info("Fetching UniProt metadata from %s...", url)
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "VeritasBio-CloudEconomicsEngine/3.1 (contact@veritasbio.ai)"},
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            _save_cached_json(clean_acc, data)
            return parse_uniprot_response(data)
    except Exception as exc:
        logger.warning("UniProt network request failed for %s (%s). Attempting local fallback...", clean_acc, exc)
        cached = _load_cached_json(clean_acc)
        if cached:
            logger.info("Using cached fallback for %s.", clean_acc)
            return parse_uniprot_response(cached)
        raise RuntimeError(f"Could not fetch protein metadata for {clean_acc} online or from local cache: {exc}") from exc
