"""Minimal client for HAPI (Heliophysics Application Programmer's Interface) servers.

Only the standard library and pandas are used. Requests use the HAPI 2
parameter names (id, time.min, time.max), which CDAWeb requires, and fall
back to the HAPI 3 names (dataset, start, stop) if a server rejects them.
"""

from __future__ import annotations

import io
import json
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

CDAWEB = "https://cdaweb.gsfc.nasa.gov/hapi"
_TIMEOUT = 180


class HapiError(RuntimeError):
    """The server answered with a HAPI status message instead of data."""


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "helios-shocks/0.1"})
    with urllib.request.urlopen(req, timeout=_TIMEOUT) as resp:
        return resp.read()


def _query(server: str, endpoint: str, v3: dict, v2: dict) -> bytes:
    # CDAWeb only understands the HAPI 2 names, and HAPI 3 servers must still accept them,
    # so those go first; the HAPI 3 names are the fallback.
    v3, v2 = v2, v3
    base = f"{server.rstrip('/')}/{endpoint}"
    try:
        return _get(f"{base}?{urllib.parse.urlencode(v3)}")
    except urllib.error.HTTPError as err:
        first = _status_message(err)
        if err.code not in (400, 404):
            raise HapiError(first) from err
    try:
        return _get(f"{base}?{urllib.parse.urlencode(v2)}")
    except urllib.error.HTTPError as err:
        raise HapiError(f"{first} | with HAPI 3 names: {_status_message(err)}") from err


def _status_message(err: urllib.error.HTTPError) -> str:
    """The HAPI status (code and message) carried in an error response body, if any."""
    try:
        body = err.read().decode("utf-8", "replace")
        status = json.loads(body).get("status", {})
        return f"HTTP {err.code}, HAPI {status.get('code')}: {status.get('message')}"
    except Exception:  # noqa: BLE001 - the body is best-effort context only
        return f"HTTP {err.code}"


def order_parameters(meta: dict, parameters: list[str]) -> list[str]:
    """Requested parameters in the dataset's own order.

    HAPI requires this: a request listing parameters in any other order is
    rejected with status 1411 (out-of-order or duplicate parameters).
    """
    known = [p["name"] for p in meta["parameters"][1:]]
    unknown = [p for p in parameters if p not in known]
    if unknown:
        raise HapiError(f"parameters not in dataset: {unknown}")
    wanted = set(parameters)
    return [name for name in known if name in wanted]


def catalog(server: str = CDAWEB) -> list[str]:
    """Dataset identifiers offered by the server."""
    return [d["id"] for d in json.loads(_get(f"{server.rstrip('/')}/catalog"))["catalog"]]


def info(dataset: str, server: str = CDAWEB) -> dict:
    """Metadata for one dataset: parameters, units, fill values, time coverage, cadence."""
    return json.loads(_query(server, "info", {"dataset": dataset}, {"id": dataset}))


def column_names(meta: dict, parameters: list[str] | None = None) -> list[str]:
    """Flat column names for a CSV response: vector parameters expand to name_0, name_1, ..."""
    names: list[str] = []
    for i, p in enumerate(meta["parameters"]):
        if i > 0 and parameters is not None and p["name"] not in parameters:
            continue
        size = int(np.prod(p.get("size", [1])))
        names.extend([p["name"]] if size == 1 else [f"{p['name']}_{k}" for k in range(size)])
    return names


def parse_csv(text: str, meta: dict, parameters: list[str] | None = None) -> pd.DataFrame:
    """Parse a HAPI CSV body into a DataFrame indexed by UTC time, with fill values as NaN."""
    names = column_names(meta, parameters)
    if not text.strip():
        return pd.DataFrame(columns=names[1:], index=pd.DatetimeIndex([], name="time", tz="UTC"))
    df = pd.read_csv(io.StringIO(text), header=None, names=names)
    df[names[0]] = pd.to_datetime(df[names[0]], utc=True, format="ISO8601")
    df = df.set_index(names[0]).rename_axis("time")
    for i, p in enumerate(meta["parameters"]):
        if i == 0 or (parameters is not None and p["name"] not in parameters):
            continue
        fill = p.get("fill")
        if fill in (None, "", "null"):
            continue
        cols = [c for c in df.columns if c == p["name"] or c.startswith(p["name"] + "_")]
        try:
            fill_value = float(fill)
        except (TypeError, ValueError):
            continue
        for c in cols:
            col = pd.to_numeric(df[c], errors="coerce")
            df[c] = col.mask(np.isclose(col, fill_value, rtol=1e-6, atol=0.0))
    return df


def fetch(dataset: str, parameters: list[str], start: str, stop: str, server: str = CDAWEB,
          meta: dict | None = None) -> pd.DataFrame:
    """Time series for the given parameters over [start, stop) (ISO 8601 UTC strings)."""
    meta = meta or info(dataset, server)
    parameters = order_parameters(meta, parameters)
    plist = ",".join(parameters)
    body = _query(
        server, "data",
        {"dataset": dataset, "parameters": plist, "start": start, "stop": stop, "format": "csv"},
        {"id": dataset, "parameters": plist, "time.min": start, "time.max": stop, "format": "csv"},
    )
    text = body.decode("utf-8")
    if text.lstrip().startswith("{"):  # a status object delivered with HTTP 200
        status = json.loads(text).get("status", {})
        raise HapiError(f"HAPI {status.get('code')}: {status.get('message')}")
    return parse_csv(text, meta, parameters)
